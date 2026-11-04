"""
source_ip_corr_MP.py
--------------------------
Instrumental polarisation correction pipeline with automatic calibrator
selection, written to run on multiple cores.

Workflow
--------
The pipeline finds the uncorrected science scans in ``red_uncorr`` and
loads calibrator scans from ``red_calibrators``. For each science scan
candidate calibrators that meet the desired criteria (e.g. within an
elevation difference of ``ELEV_DIFF`` degrees) are identified. The
pipeline is split into five stages which can be run independently via
the ``STAGES`` configuration variable:

1. correct     — run corr_instr_pol_1scan.piic for all scan/calibrator/array
                combinations in parallel (N_CORES workers). Each worker
                runs PIIC in an isolated tmp_i directory with its own
                $HOME/.gag/logs to avoid log-file collisions between
                concurrent PIIC processes.
2. select      — load corrected maps, compute Q/U residual metrics, and
                 identify the best calibrator per scan
3. refine      — iteratively refine the per-scan calibrator assignment by
                 minimising the stack standard deviation
4. write_lists — write LIST files for the best-corrected maps and rgw files
5. combine     — run ``combine_maps.piic`` to co-add the best corrected maps

Usage
-----
Set ``SOURCE``, ``REPO_ROOT`` and ``STAGES`` at the top of the file, then
run::

    python source_ip_corr_MP.py

To run only a subset of stages set::

    STAGES = {"correct"}

To run all stages set::

    STAGES = {"all"}

But check README as to whether you want to run with refine!!    
"""

import os
import sys
import subprocess
import random
import shutil
from pathlib import Path
from astropy.io import fits
from multiprocessing import Pool
import matplotlib
matplotlib.use("Agg")          
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
from tqdm import tqdm

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------

# Target
SOURCE = "???"

# Stages to run — subset of {"correct", "select", "refine", "write_lists",
# "combine"} or {"all"} to run the full pipeline
STAGES = {"correct", "select", "write_lists", "combine"}

# Save diagnostic figures
SAVE_FIGURES = True

# Maximum elevation difference (degrees) between calibrator and science scan
ELEV_DIFF = 10

# Flag scans with residuals above this percentage of Imax
RESIDUAL_THRESHOLD = 5

# Metric used to select and refine the best calibrator per scan
# Options: "stdv", "max", "sum" — evaluated over the central CROP_FRAC
# of each map
METRIC = "stdv"
CROP_FRAC = 0.5

# Tukey taper radii for the science map [arcsec] — inner, outer
TAPER_RADII = (250, 500)

# Diagnostic plot colour scale — units are percentage of Stokes/I (leakage)
VLIM = 2

# Set to a specific science scan ID to process only that scan
# (useful for testing)
SCIENCE_ID_FILTER = None

# Number of cores to run on
N_CORES = 50

# Leave as None to auto-detect the repo root. Set manually if PIIC/GILDAS
# gives "unbalanced quote count" errors from long paths - see README.
REPO_ROOT = None
# REPO_ROOT = Path.home() / "redNIKA2pol"


# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------
try:
    plt.style.use("standard")
except OSError:
    pass

# ---------------------------------------------------------------------------
# Derived directory paths
# ---------------------------------------------------------------------------
def find_repo_root(marker="setup.sh"):
    """Walk up from this file's location to find the repo root, identified
    by the presence of ``marker`` (default: setup.sh)."""
    path = Path(__file__).resolve()
    for parent in path.parents:
        if (parent / marker).exists():
            return parent
    raise FileNotFoundError(
        f"Could not find repo root (no {marker} found in any parent directory)"
    )

if REPO_ROOT is None:
    repo_root = find_repo_root()
else:
    repo_root = REPO_ROOT
    
base                = repo_root / "reductions" / SOURCE / "03_corr"
red_dir             = base / "red"
red_uncorr_dir      = base / "red_uncorr"
red_calibrators_dir = base / "red_calibrators"
fig_dir             = repo_root / "figures" / SOURCE 

# ---------------------------------------------------------------------------
# Pipeline command wrappers
# ---------------------------------------------------------------------------

def run_combination(maps_list, rgw_list, outpath):
    """Run ``combine_maps.piic`` to co-add a set of maps."""
    command = (
        f"source ~/.bashrc; gagpiic; piic -nl @ combine_maps "
        f"red red {maps_list} {rgw_list} {outpath} ; "
        f"rm -f PIIC*"
    )
    log_path = base / "combine.log"
    print(f"  Running combine_maps")

    with open(log_path, "w") as log:
        subprocess.run(
            command, shell=True, executable="/bin/bash",
            cwd=base,                                  # run from 03_corr/
            stdout=log, stderr=subprocess.STDOUT,
        )

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def load_iqu(path):
    """Load Stokes I, Q, U planes from a FITS file."""
    with fits.open(path, memmap=False) as hdul:
        I, Q, U = hdul[0].data
    return I, Q, U

def crop_center(data, crop_frac = CROP_FRAC):
    """Return the central rectangular crop of a 2D array."""
    ny, nx = data.shape
    half_x = int(nx * crop_frac) // 2
    half_y = int(ny * crop_frac) // 2
    return data[
        ny // 2 - half_y : ny // 2 + half_y,
        nx // 2 - half_x : nx // 2 + half_x,
    ]

def plot_map(ax, data, title, cbar_label):
    """Plot the central crop of a 2D map on ``ax``."""
    cropped = crop_center(data)
    im = ax.imshow(cropped, origin="lower", cmap="PiYG", vmin=-VLIM, vmax=VLIM)
    ax.tick_params(labelbottom=False, labelleft=False)
    ax.text(
        0.95, 0.1, 
        f"max = {np.nanmax(cropped):.1f}%", 
        transform=ax.transAxes, 
        ha='right')
    if title:
        ax.set_title(title)
    cbar = ax.get_figure().colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(cbar_label)

def highlight_axes(axes, color, linewidth = 3):
    """Draw a coloured border around each Axes in ``axes``."""
    for ax in axes:
        for spine in ax.spines.values():
            spine.set_edgecolor(color)
            spine.set_linewidth(linewidth)

def find_science_scans(directory):
    """Find science scans in ``directory``.

    Expects files named ``...-1-{scanID}_...0s.fits`` (array 1); the matching
    array-3 file is assumed to have ``-3-`` in place of ``-1-``.
    """
    science_scans = {}
    for path in sorted(directory.glob("*-1-*0s.fits")):
        parts = path.name.split("-1-")
        science_id = parts[1].split("_")[0]
        with fits.open(path) as hdul:
            header = hdul[0].header
        science_scans[science_id] = {
            "path_a1": path,
            "path_a3": Path(str(path).replace("-1-", "-3-")),
            "elev":   header.get("ELEND", np.nan)
        }
    print(f"Found {len(science_scans)} science scans in {directory}")
    return science_scans

def build_calib_catalogue(directory):
    """Build catalogue of calibrator scans from ``directory``.

    Expects files named ``{OBJECT}_{scanID}_a1_...n10.fits``
    (e.g. ``Neptune_20250214s130_a1_...``); the matching array-3 file is
    assumed to have ``_a3_`` in place of ``_a1_``.
    """
    catalogue = {}
    skipped = []
    for path in sorted(directory.glob("*_a1_*n10.fits")):
        calib_id = path.name.split("_")[1]
        path_a3 = Path(str(path).replace("_a1_", "_a3_"))
        if not path_a3.exists():
            skipped.append(calib_id)
            continue
        with fits.open(path) as hdul:
            header = hdul[0].header
        catalogue[calib_id] = {
            "path_a1": path,
            "path_a3": path_a3,
            "object": header.get("OBJECT", ""),
            "elev":   header.get("ELEND", np.nan),
        }
    if skipped:
        print(
            f"WARNING: skipping {len(skipped)} calibrator(s) missing "
            f"their a3 file: {skipped}"
        )
    print(
        f"Calibrator catalogue: {len(catalogue)} entries found in {directory}"
    )
    return catalogue

def find_candidate_calibrators(
    science_scans, 
    calib_catalogue, 
    elev_diff=ELEV_DIFF
):
    """Assign candidate calibrators to each science scan.

    A calibrator is a candidate if its elevation is within ``elev_diff``
    degrees of the science scan. The entry ``'Uncorrected'`` is
    always prepended as a baseline option.
    """
    for science_id, science_meta in science_scans.items():
        science_meta['candidate_calibrators'] = ['Uncorrected'] + [
            calib_id 
            for calib_id, calib_meta in calib_catalogue.items()
            if abs(calib_meta["elev"] - science_meta["elev"]) <= elev_diff
        ]
        print(
            f"  {science_id}: elev={science_meta['elev']:.2f}° "
            f"-> {len(science_meta['candidate_calibrators'])-1} "
            f"matching calibrators"
        )
    return science_scans

def compute_stack_std(science_scans_corr, assignment, exclude=frozenset()):
    """Compute the total stack standard deviation for corrected Q and U.

    For each science scan (excluding those in ``exclude``) the Q and U maps
    are normalised by their Stokes I peak, cropped to the central region,
    and stacked. The function returns the sum of pixel-wise standard
    deviations across the four array/Stokes combinations.
    """
    qu_stacks = {"Q_a1": [], "U_a1": [], "Q_a3": [], "U_a3": []}

    for science_id, calib_id in assignment.items():
        if science_id in exclude:
            continue
        I1, Q1, U1 = science_scans_corr[science_id][calib_id]["iqu_a1"]
        I3, Q3, U3 = science_scans_corr[science_id][calib_id]["iqu_a3"]
        qu_stacks["Q_a1"].append(crop_center(Q1 / np.nanmax(I1)))
        qu_stacks["U_a1"].append(crop_center(U1 / np.nanmax(I1)))
        qu_stacks["Q_a3"].append(crop_center(Q3 / np.nanmax(I3)))
        qu_stacks["U_a3"].append(crop_center(U3 / np.nanmax(I3)))

    return sum(
        np.nansum(np.nanstd(np.array(stack), axis=0))
        for stack in qu_stacks.values()
    )

# ---------------------------------------------------------------------------
# Internal parallelisation functions
# ---------------------------------------------------------------------------

REAL_HOME = Path(os.environ["HOME"])
KEEP_TMP = False

def _build_correction_tasks(science_scans, calib_catalogue):
    """Return list of tasks for workers to perform individually. 
    """
    tasks = []
    for science_id, science_meta in science_scans.items():
        for calib_id in science_meta["candidate_calibrators"]:
            if calib_id == "Uncorrected":
                continue
            calib_meta = calib_catalogue[calib_id]
            for array in ("a1", "a3"):
                science_p = science_meta[f"path_{array}"]
                calib_p = calib_meta[f"path_{array}"]
                final_p = red_dir / f"{science_p.stem}_{calib_id}_IP_corrected.fits"
                if final_p.exists():
                    continue
                tasks.append((science_p, science_id, calib_p, calib_id, final_p))
    return [(i, *task) for i, task in enumerate(tasks)]  # index each task

def _make_worker_home(tmp_dir):
    """Build a per-worker $HOME that mirrors the real one, but with an
    isolated .gag/logs directory to avoid piic/cube's log-rename collisions
    when multiple workers run concurrently. Allows for mulitple reductions
    at the same time :). 
    """
    worker_home = tmp_dir / "home"
    worker_gag = worker_home / ".gag"
    worker_gag.mkdir(parents=True, exist_ok=True)

    real_gag = REAL_HOME / ".gag"
    for item in real_gag.iterdir():
        if item.name == "logs":
            continue  # skip — give this its own fresh directory below
        link = worker_gag / item.name
        if not link.exists() and not link.is_symlink():
            os.symlink(item.resolve(), link)

    (worker_gag / "logs").mkdir(exist_ok=True)
    return worker_home

def _run_one_correction(task):
    """Run corr_instr_pol_1scan.piic for one (science, calibrator, array) task.
    Runs PIIC inside an isolated tmp_i directory (own cwd and own
    $HOME/.gag/logs) to avoid collisions between concurrent workers.
    """
    i, science_p, science_id, calib_p, calib_id, final_p = task
    tmp_dir = base / f"tmp_{i}" 
    tmp_dir.mkdir(parents=True, exist_ok=True)
    tmp_red = tmp_dir / "red"
    tmp_red.mkdir(parents=True, exist_ok=True)

    for src in (science_p, calib_p):
        dst = tmp_dir / "red" / src.name
        if not dst.exists() and not dst.is_symlink():
            os.symlink(src.resolve(), dst)

    worker_home = _make_worker_home(tmp_dir)
    env = os.environ.copy()
    env["HOME"] = str(worker_home)

    command = (
        f"source {REAL_HOME}/.bashrc; gagpiic; "
        f"piic -nl @ {base}/corr_instr_pol_1scan "
        f"{science_p.name} {calib_p.name} {calib_id} "
        f"{TAPER_RADII[0]} {TAPER_RADII[1]} ;"
        f"rm PIIC*"
    )

    log_path = tmp_dir / "piic_stdout.log2"

    success = False

    try:
        with open(log_path, "w") as log:
            subprocess.run(
                command, shell=True, executable="/bin/bash",
                cwd=tmp_dir, env=env, check=True,
                stdout=log, stderr=subprocess.STDOUT,   
            )

        out = tmp_red / final_p.name

        if out.exists():
            success = True
            out.replace(red_dir / final_p.name)

    except subprocess.CalledProcessError:
        pass

    finally:
        if success and not KEEP_TMP:
            shutil.rmtree(tmp_dir, ignore_errors=True)

    return {"science_id": science_id, "science": science_p.stem,
            "calib": calib_id, "success": success, "log": str(log_path)}

# ---------------------------------------------------------------------------
# Internal figure functions
# ---------------------------------------------------------------------------   

def _figures_comparison(
    science_scans,
    science_scans_corr,
    calib_catalogue,
    assignment,
    tag,
    bad_scans=frozenset(),
    changed_from=None,
):
    """Save one figure per Stokes parameter (Q and U) comparing the
    uncorrected and corrected maps for every science scan.

    Rows are science scans; columns are
    [A1 uncorrected, A1 corrected, A3 uncorrected, A3 corrected].
    The corrected map is taken from ``assignment[science_id]``.
    """
    changed_from = changed_from or {}
    science_ids = list(science_scans.keys())
    n_rows = len(science_ids)

    # Agg can't save images taller than 2^16 px, so lower the dpi when
    # there are many scans
    dpi = min(300, int(60000 / (4 * n_rows)))

    for stokes, idx in (("Q", 1), ("U", 2)):
        col_titles = [
            f"A1 {stokes} uncorrected", f"A1 {stokes} corrected",
            f"A3 {stokes} uncorrected", f"A3 {stokes} corrected",
        ]
        cbar_label = rf"${stokes}/I_\mathrm{{max}}$ (%)"

        fig, axes = plt.subplots(
            n_rows, 4, figsize=(16, 4 * n_rows), squeeze=False
        )

        for row, science_id in enumerate(science_ids):
            calib_id = assignment[science_id]
            uncorr = science_scans_corr[science_id]["Uncorrected"]
            corr = science_scans_corr[science_id][calib_id]

            # Column order: A1 uncorr, A1 corr, A3 uncorr, A3 corr
            panels = []
            for array in ("a1", "a3"):
                for entry in (uncorr, corr):
                    iqu = entry[f"iqu_{array}"]
                    panels.append(iqu[idx] / np.nanmax(iqu[0]) * 100)

            for col, panel in enumerate(panels):
                plot_map(
                    axes[row, col],
                    panel,
                    title=col_titles[col] if row == 0 else None,
                    cbar_label=cbar_label,
                )

            # Row label + annotations
            axes[row, 0].set_ylabel(science_id, fontsize=14)
            science_elev = science_scans[science_id]["elev"]
            if calib_id == "Uncorrected":
                corr_str = "calib: none"
            else:
                calib_elev = calib_catalogue[calib_id]["elev"]
                corr_str = f"calib: {calib_id}\nelev = {calib_elev:.1f}°"
            if science_id in changed_from:
                corr_str += f"\nwas: {changed_from[science_id]}"

            for col in (0, 2):
                axes[row, col].text(
                    0.05, 0.95, f"elev = {science_elev:.1f}°",
                    transform=axes[row, col].transAxes, va="top",
                )
            for col in (1, 3):
                axes[row, col].text(
                    0.05, 0.95, corr_str,
                    transform=axes[row, col].transAxes, va="top",
                )

            if science_id in bad_scans:
                highlight_axes(axes[row], color="tab:red")
            elif science_id in changed_from:
                highlight_axes(axes[row], color="tab:orange")

        fig.suptitle(
            f"{SOURCE}: Stokes {stokes}, uncorrected vs corrected "
            f"({tag}, {METRIC})",
            fontsize=16,
            y=1.0 + (0.3 / (4 * n_rows)),
        )
        fig.tight_layout()
        out_path = (
            fig_dir / f"{SOURCE}_{stokes}_uncorr_vs_corr_{tag}_{METRIC}.png"
        )
        fig.savefig(out_path, bbox_inches="tight", dpi=dpi)
        plt.close(fig)
        print(f"  Saved: {out_path}")

# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def stage_correct(science_scans, calib_catalogue):
    print("\n=== Stage 1: running IP corrections (parallel) ===")
    tasks = _build_correction_tasks(science_scans, calib_catalogue)
    if not tasks:
        print("  Nothing to do: all corrections already exist.")
        return

    totals = Counter(t[2] for t in tasks) 
    done = Counter()
    failed = Counter()
    scan_bars = {}
    failures = []

    print(f"  {len(tasks)} corrections over {len(totals)} science IDs "
      f"({2 * len(totals)} A1/A3 maps)")

    sys.stdout.reconfigure(line_buffering=True)
    interactive = sys.stderr.isatty()

    with Pool(N_CORES) as pool, \
         tqdm(total=len(tasks), desc="All scans", unit="run",
              position=0, smoothing=0,
              disable=not interactive) as overall:

        for result in pool.imap_unordered(_run_one_correction, tasks):
            sid = result["science_id"]

            if sid not in scan_bars:
                scan_bars[sid] = tqdm(
                    total=totals[sid], desc=f"  {sid}", unit="run",
                    leave=False, disable=not interactive,
                )
            bar = scan_bars[sid]

            bar.update(1)
            overall.update(1)
            done[sid] += 1

            if not result["success"]:
                failed[sid] += 1
                failures.append(result)
                bar.set_postfix(failed=failed[sid])

            # scan finished: remove its bar and print one summary line
            if done[sid] == totals[sid]:
                bar.close()
                del scan_bars[sid]
                status = "all ok" if failed[sid] == 0 else f"{failed[sid]} FAILED"
                overall.write(f"  {sid}: {totals[sid]} corrections, {status}")

    print(f"  Done: {len(tasks) - len(failures)} succeeded, {len(failures)} failed")
    if failures:
        fail_path = base / "failed_corrections.txt"
        with open(fail_path, "w") as f:
            for r in failures:
                f.write(f"{r['science']}\t{r['calib']}\t{r['log']}\n")
        print(f"  Failed list written to {fail_path}")
    
def stage_select(science_scans, calib_catalogue):
    """Stage 2 — load corrected maps and select the best calibrator per
    scan.

    For each science scan and each candidate calibrator, loads the
    IP-corrected maps, computes a residual metric (``METRIC``, summed over
    Q and U for both arrays), and computes the peak absolute residual.
    The calibrator with the lowest metric is selected as ``'best_calib'``.

    Scans whose *selected* best-calibrator correction still has a
    ``max_residual`` exceeding ``RESIDUAL_THRESHOLD`` are flagged as
    ``bad_scans`` and a warning is printed. This is a diagnostic flag only:
    - In `stage_refine`, bad scans are excluded from the stack std
      calculation and their calibrator assignment is left unchanged.
    - They are NOT excluded from `stage_write_lists` or `stage_combine` —
      they are still written to the LIST files and included in the final
      combined map. Review flagged scans manually; there is currently no
      option to drop them from the final product automatically.
    """
    print("\n=== Stage 2: selecting best calibrators ===")

    metric_funcs = {
        "stdv": np.nanstd,
        "max": lambda x: np.nanmax(np.abs(x)),
        "sum": lambda x: np.nansum(np.abs(x)),
    }

    if METRIC not in metric_funcs:
        raise ValueError(f'METRIC must be one of {list(metric_funcs.keys())}')

    science_scans_corr = {}

    for science_id, science_meta in science_scans.items():
        science_scans_corr[science_id] = {}
        
        for calib_id in science_meta['candidate_calibrators']:

            if calib_id == 'Uncorrected':
                path_a1 = science_meta['path_a1']
                path_a3 = science_meta['path_a3']
            else:
                corr_name_a1 = (
                    science_meta['path_a1'].stem
                    + f"_{calib_id}_IP_corrected.fits"
                )
                corr_name_a3 = (
                    science_meta['path_a3'].stem
                    + f"_{calib_id}_IP_corrected.fits"
                )
                path_a1 = red_dir / corr_name_a1
                path_a3 = red_dir / corr_name_a3

            iqu_a1 = load_iqu(path_a1)
            iqu_a3 = load_iqu(path_a3)

            I1, Q1, U1 = iqu_a1
            I3, Q3, U3 = iqu_a3

            # Cropped normalised panels
            norm_crops = [
                crop_center(Q1 / np.nanmax(I1) * 100),
                crop_center(U1 / np.nanmax(I1) * 100),
                crop_center(Q3 / np.nanmax(I3) * 100),
                crop_center(U3 / np.nanmax(I3) * 100),
            ]

            metric_fn = metric_funcs[METRIC]
            metric = sum(metric_fn(nc) for nc in norm_crops)
            max_residual = max(np.nanmax(np.abs(nc)) for nc in norm_crops)

            science_scans_corr[science_id][calib_id] = {
                'path_a1': path_a1,
                'path_a3': path_a3,
                'iqu_a1': iqu_a1,
                'iqu_a3': iqu_a3,
                'metric': metric,
                'max_residual': max_residual
            }

        # Select calibrator with lowest metric
        best_calib = min(
            science_scans_corr[science_id],
            key=lambda c: science_scans_corr[science_id][c]['metric']
        )

        best = science_scans_corr[science_id][best_calib]
        print(
            f"  {science_id}: best calibrator = {best_calib} "
            f"({METRIC} = {best['metric']:.2f}, "
            f"max residual = {best['max_residual']:.2f}%)"
        )
        science_scans[science_id]['best_calib'] = best_calib

    # Identifying 'bad scans' - where remaining IP signal is above threshold
    bad_scans = {
        sid for sid, meta in science_scans.items()
        if science_scans_corr[sid][meta['best_calib']]['max_residual']
        > RESIDUAL_THRESHOLD
    }
    if bad_scans:
        print(
            f"WARNING: {len(bad_scans)} scans flagged for "
            f"high residual: {bad_scans}"
        )

    if SAVE_FIGURES:
        assignment = {
            sid: meta["best_calib"] for sid, meta in science_scans.items()
        }
        _figures_comparison(
            science_scans, science_scans_corr, calib_catalogue,
            assignment, tag="select", bad_scans=bad_scans,
        )

    return science_scans, science_scans_corr, bad_scans


def stage_refine(
    science_scans,
    science_scans_corr,
    calib_catalogue,
    bad_scans,
):
    """Stage 3 — refine the per-scan calibrator assignment by iterative
    optimisation.

    Starting from the per-scan best calibrators identified in
    :func:`stage_select`, this stage iteratively reassigns calibrators to
    minimise :func:`compute_stack_std` — the sum of pixel-wise standard
    deviations of normalised Q and U maps across all scans.

    At each iteration the science scans are visited in a random order. For
    each scan every candidate calibrator is trialled in turn; the assignment
    is updated if it reduces the global stack standard deviation. Iteration
    continues until a full pass produces no changes (convergence). Scans in
    ``bad_scans`` are excluded from the std calculation and their calibrator
    assignment is left unchanged.

    The final assignment is written back into ``science_scans`` as
    ``'best_calib'``. A summary of which scans changed calibrator is
    printed, and diagnostic figures are saved if ``SAVE_FIGURES`` is True.
    """

    # Initialise assignment from stage_select
    assignment = {
        science_id: science_meta["best_calib"]
        for science_id, science_meta in science_scans.items()
    }

    assignment_initial = assignment.copy()
    std_initial = compute_stack_std(
                    science_scans_corr,
                    assignment,
                    exclude=bad_scans,
                )
    print(f"  Initial std: {std_initial:.3f}")

    iteration = 0
    science_ids = list(science_scans.keys())

    while True:
        iteration += 1
        changed = False
        random.shuffle(science_ids)

        for science_id in science_ids:

            if science_id in bad_scans:
                continue

            science_meta = science_scans[science_id]
            original = assignment[science_id]
            best_calib = original
            best_std = compute_stack_std(
                science_scans_corr,
                assignment,
                exclude=bad_scans,
            )

            for calib_id in science_meta["candidate_calibrators"]:

                # Skip if already using this calibrator
                if calib_id == original:
                    continue

                # Temporarily swap this scan to calib_id
                assignment[science_id] = calib_id
                std = compute_stack_std(
                    science_scans_corr,
                    assignment,
                    exclude=bad_scans,
                )

                # Update if total std of stack decreased
                if std < best_std:
                    best_std = std
                    best_calib = calib_id
                
            assignment[science_id] = best_calib

            if best_calib != original:
                changed = True

        print(
            f"  Iteration {iteration}: std = {best_std:.3f}, "
            f"changed = {changed}"
        )
        if not changed:
            print("  Converged.")
            break

    # Write results back into science_scans
    for science_id, calib_id in assignment.items():
        science_scans[science_id]["best_calib"] = calib_id
        print(f"  {science_id}: best calibrator = {calib_id}")

    assignment_final = assignment.copy()

    science_ids = list(science_scans.keys())
    changed_ids = {
        sid for sid in science_ids
        if assignment_initial[sid] != assignment_final[sid]
    }
    
    print(f"  {len(changed_ids)} scans changed calibrator during refinement")
    for changed_id in changed_ids:
        old = assignment_initial[changed_id]
        new = assignment_final[changed_id]
        print(f"{changed_id}: {old}->{new}")

    if SAVE_FIGURES:
        changed_from = {
            sid: assignment_initial[sid] for sid in changed_ids
        }
        _figures_comparison(
            science_scans, science_scans_corr, calib_catalogue,
            assignment_final, tag="refine", bad_scans=bad_scans,
            changed_from=changed_from,
        )

    return science_scans

def stage_write_lists(science_scans, science_scans_corr):
    """Stage 4 — write LIST files of best-corrected maps and rgw weight
    files.
 
    For each science scan the best-calibrator corrected maps (arrays A1 and
    A3) are written to ``all_corr_maps.LIST``. The corresponding rgw weight
    files share the uncorrected map name with ``rgw`` appended
    (``..._0s.fits`` -> ``..._0srgw.fits``); they are located in
    ``red_uncorr_dir``, symlinked into ``red_dir``, and written to
    ``all_rgw.LIST``. Both files are written to ``base``.
    """
    print("\n=== Stage 4: writing LIST files ===")
 
    maps_path = base / "all_corr_maps.LIST"
    rgw_path  = base / "all_rgw.LIST"
 
    with open(maps_path, "w") as f_maps, open(rgw_path, "w") as f_rgw:
        for science_id, science_meta in science_scans.items():
            best_calib = science_meta['best_calib']
            for array in ("a1", "a3"):
                corr_path = science_scans_corr[science_id][best_calib][f"path_{array}"]
                rgw_name  = f"{science_meta[f'path_{array}'].stem}rgw.fits"
 
                f_maps.write(f"{corr_path.name}\n")
                f_rgw.write(f"{rgw_name}\n")
 
                dst = red_dir / rgw_name
                if not dst.exists() and not dst.is_symlink():
                    os.symlink(red_uncorr_dir / rgw_name, dst)
 
    print(f"  Saved map list : {maps_path}")
    print(f"  Saved rgw list : {rgw_path}")
 
def stage_combine():
    """Stage 5 — co-add the best IP-corrected maps using
    ``combine_maps.piic``.
 
    Reads the LIST files written by :func:`stage_write_lists` and passes
    them to :func:`run_combination`. The combined output is written as
    ``{SOURCE}_IP_corrected.fits``.
    """
    print("\n=== Stage 5: combining best-calibrated maps ===")
 
    maps_path = "all_corr_maps.LIST"
    rgw_path  = "all_rgw.LIST"
    combined_output = f"{SOURCE}_IP_corrected.fits"
    out = red_dir / combined_output
    log_path = base / "combine.log"
 
    out.unlink(missing_ok=True)   # remove any old result first
 
    run_combination(maps_path, rgw_path, combined_output)
 
    if out.exists():
        log_path.unlink(missing_ok=True)
        print(f"  Combined map written to {out}")
    else:
        print(f"  WARNING: {out} not found, log kept at {log_path}")
     
# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Find scans, build the calibrator catalogue, and run selected stages."""
    print("Source:", SOURCE)
    active = {
        "correct",
        "select",
        "refine",
        "write_lists",
        "combine",
    } if "all" in STAGES else STAGES
    print("Active stages:", active)

    # -- Science and calibrator scan finding --------------------------------
    science_scans = find_science_scans(red_uncorr_dir)
    calib_catalogue = build_calib_catalogue(red_calibrators_dir)
    science_scans = find_candidate_calibrators(science_scans, calib_catalogue)

    for meta in science_scans.values():
        for array in ("a1", "a3"):
            src = meta[f"path_{array}"]
            dst = red_dir / src.name
            if not dst.exists() and not dst.is_symlink():
                os.symlink(src.resolve(), dst)
    
    # (for testing)
    if SCIENCE_ID_FILTER is not None:
        science_scans = {SCIENCE_ID_FILTER: science_scans[SCIENCE_ID_FILTER]}
        print(f"Filtering to single science scan: {SCIENCE_ID_FILTER}")

    # -- Stage 1: correct 
    if "correct" in active:
        stage_correct(science_scans, calib_catalogue)

    # -- Stage 2: select best calibrator
    if "select" in active:
        science_scans, science_scans_corr, bad_scans = stage_select(
            science_scans, calib_catalogue
        )

    # -- Stage 3: refine through minimising std
    if "refine" in active:
        if "select" not in active:
            raise ValueError("'refine' requires 'select'")
        else:
            science_scans = stage_refine(
                science_scans,
                science_scans_corr,
                calib_catalogue,
                bad_scans,
            )

    # -- Stage 4: write lists of best corrected scans
    if "write_lists" in active:
        if "select" not in active:
            raise ValueError("'write_lists' requires at least 'select'")
        else:
            stage_write_lists(science_scans, science_scans_corr)

    # -- Stage 5: combine
    if "combine" in active:
        stage_combine()

    print("\nDone.")

if __name__ == "__main__":
    main()
