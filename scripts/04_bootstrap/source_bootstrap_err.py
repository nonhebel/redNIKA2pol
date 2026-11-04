"""
source_bootstrap_err.py
--------------------------
Bootstrap uncertainty across science maps.

Overview
--------
Starting from a set of individually reduced and corrected cubes
(Stokes I, Q, U), this script:

  1. If BSM > 1, smooths each scan to the target BSM using PIIC and combines
     the full list of smoothed scans into red/{SOURCE}_IP_corrected_bsm{BSM}.fits.
     This smoothed science map is built the same way as the bootstrap
     realisations (smooth each scan, then combine), so it is consistent with
     the bootstrap covariance map. At BSM = 1 the native-resolution science
     map comes from 03_corr, and no combination is done here.
  2. Draws N_BOOTSTRAP bootstrap resamples of the original scan list and
     combines each resample into a single map via PIIC. Each bootstrap
     realisation uses a fixed random seed equal to its index, so results are
     fully reproducible.
  3. Computes a per-pixel 3×3 covariance matrix over Stokes I, Q, U from
     the full set of bootstrap realisations.

Stages
------
bootstrap - Draw resamples and produce combined maps (parallelised).
covariance - Compute the per-pixel covariance cube from all bootstrap maps.

Usage
-----
Set configuration at the top of the script and then run with

  python source_bootstrap_err.py

Pass STAGES = {"all"} to run everything, or a subset such as
{"bootstrap"} to run a single stage.

"""

import os
import subprocess
from pathlib import Path
import shutil
from multiprocessing import Pool
import tqdm
import numpy as np
from astropy.io import fits

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------
SOURCE          = "???"                                                             # Name of target
BSM             = 1                                                                 # Desired smoothing (BSM in PIIC, 1 for native resolution)
N_BOOTSTRAP     = 5000                                                              # Number of bootstrap realisations
N_CORES         = 50                                                                # Number of cores to run on
STAGES          = {"all"}                                                           # Pipeline stages to execute

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

repo_root = find_repo_root()

base                = repo_root / "reductions" / SOURCE / "04_bootstrap"
red_dir             = base / "red"
combined_dir        = base / f"combined_bsm{BSM}"
combined_dir.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# PIIC command wrappers
# ---------------------------------------------------------------------------

def run_bsm(in_dir, out_dir, map_path, bsm_val):
    """
    Execute the bsm.piic script to smooth a map to BSM value.
    """
    command = (
        f"source ~/.bashrc; gagpiic; piic -nl @ bsm "
        f"{in_dir} {out_dir} {map_path.name} {bsm_val}; "
        f"rm PIIC*"
    )
    subprocess.run(
        command, shell=True, executable="/bin/bash", cwd=base,
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

def run_combination(work_dir, in_dir, out_dir, maps_list, rgw_list, outpath):
    """
    Execute the combine_maps.piic script to combine a set of maps in list.
    """
    script  = os.path.relpath(base / "combine_maps", work_dir)
    in_rel  = os.path.relpath(in_dir, work_dir)
    out_rel = os.path.relpath(out_dir, work_dir)

    command = (
        f"source ~/.bashrc; gagpiic; piic -nl @ {script} "
        f"{in_rel} {out_rel} {Path(maps_list).name} {Path(rgw_list).name} {Path(outpath).name} ; "
        f"rm PIIC*"
    )
    subprocess.run(
        command, shell=True, executable="/bin/bash", cwd=work_dir,
        stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )

# ---------------------------------------------------------------------------
# Bootstrapping procedure
# ---------------------------------------------------------------------------

def run_bootstrap(args):
    """
    Run a single bootstrap realisation.

    Returns True on success and False on failure, so that one failed
    realisation does not stop the others. The run_dir of a failed
    realisation is kept for inspection.
    """

    i, base, maps_list, rgw_list, combined_dir = args

    # Skip realisations that are already completed to allow safe restart
    final_combined = combined_dir / f"combined_{i}.fits"
    if final_combined.exists():
        return True

    try:
        # Separate working directory per realisation avoids cross-contamination
        run_dir = base / f"run_{i}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # Draw bootstrap indices
        n   = len(maps_list)
        rng = np.random.default_rng(seed=i)  # random seed so fully reproducible per realisation
        idx = rng.choice(n, size=n, replace=True)

        maps_list_i = run_dir / 'maps.LIST'
        rgw_list_i = run_dir / 'rgw.LIST'
        combined_i = run_dir / 'combined.fits'

        # Write lists of bootstrapped maps
        np.savetxt(maps_list_i, maps_list[idx], fmt="%s")
        np.savetxt(rgw_list_i, rgw_list[idx], fmt="%s")

        # Combine bootstrapped maps (written to run_dir first, then moved, so a
        # killed job never leaves a partial file in combined_dir)
        run_combination(run_dir, red_dir, run_dir, maps_list_i, rgw_list_i, combined_i)

        shutil.move(combined_i, final_combined)
        shutil.rmtree(run_dir)

    except Exception:
        return False

    return True

# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def stage_bootstrap(maps_list, rgw_list, base, combined_dir, n_bootstrap, n_cores):
    """
    Stage 1 - generate bootstrap realisations in parallel.
    """
    print("\n=== Stage 1: bootstrapping over individual corrected maps ===")

    # Only submit realisations that are not yet done
    todo = [
        i for i in range(n_bootstrap)
        if not (combined_dir / f"combined_{i}.fits").exists()
    ]
    n_done = n_bootstrap - len(todo)
    print(f"  {n_done} / {n_bootstrap} realisations already done, {len(todo)} to run")
    if not todo:
        return

    args = [
        (i, base, maps_list, rgw_list, combined_dir)
        for i in todo
    ]

    n_failed = 0
    with Pool(processes=n_cores) as p:
        for ok in tqdm.tqdm(
            p.imap_unordered(run_bootstrap, args),
            total=len(todo),
            desc="Bootstrap realisations",
            unit="map",
        ):
            n_failed += not ok

    n_ok = len(todo) - n_failed
    print(f"  Completed {n_ok} / {len(todo)} this run ({n_done + n_ok} / {n_bootstrap} in total)")
    if n_failed:
        print(f"  {n_failed} realisation(s) failed - rerun the script to retry them")

def stage_covariance(combined_dir, base, bsm):
    """
    Stage 2 - compute a per-pixel (3x3) covariance matrix over Stokes (I, Q, U)
    from all bootstrap realisations.
    """
    print("\n=== Stage 2: calculating pixel-by-pixel covariance map ===")

    boot_maps = []
    for path in combined_dir.glob("*.fits"):
        with fits.open(path, memmap=False) as hdul:
            data = hdul[0].data
        boot_maps.append(data)

    boot_maps = np.array(boot_maps)

    ny, nx = boot_maps.shape[2], boot_maps.shape[3]
    cov_map = np.zeros([ny, nx, 3, 3])

    for y in range(ny):
        for x in range(nx):
            samples = boot_maps[:, :, y, x]
            mask = np.all(np.isfinite(samples), axis=1)
            if mask.sum() > 1:
                cov_map[y, x] = np.cov(samples[mask].T)
            else:
                cov_map[y, x] = np.nan

    cov_hdu = fits.PrimaryHDU(cov_map)
    cov_hdu.writeto(base / f"cov_map_bsm{bsm}.fits", overwrite=True)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    """Smooth scans and combine them if need be, then run selected stages."""
    print("Source:", SOURCE)
    active = {"bootstrap", "covariance"} if "all" in STAGES else STAGES
    print("Active stages:", active)

    maps_list = np.loadtxt(base / "all_maps.LIST", dtype=str)
    rgw_list = np.loadtxt(base / "all_rgw.LIST", dtype=str)

    if BSM != 1:
        # Smoothing individual scans
        bsm_maps_list = base / f"all_maps_bsm{BSM}.LIST"
        bsm_rgw_list = base / f"all_rgw_bsm{BSM}.LIST"

        print(f"  Smoothing scans to BSM {BSM}...")

        with open(bsm_maps_list, "w") as f_maps, open(bsm_rgw_list, "w") as f_rgw:

            for map_path, rgw_path in zip(maps_list, rgw_list):

                bsm_map = map_path.replace(".fits", f"_bsm{BSM}.fits")
                bsm_rgw = rgw_path.replace(".fits", f"_bsm{BSM}.fits")

                # Checking if bsm files already exist
                bsm_map_out = red_dir / bsm_map
                if not bsm_map_out.exists():
                    run_bsm('red', 'red', red_dir / map_path, BSM)
                bsm_rgw_out = red_dir / bsm_rgw
                if not bsm_rgw_out.exists():
                    run_bsm('red', 'red', red_dir / rgw_path, BSM)

                f_maps.write(f"{bsm_map}\n")
                f_rgw.write(f"{bsm_rgw}\n")

        print(f"  Saved map list : {bsm_maps_list}")
        print(f"  Saved rgw list : {bsm_rgw_list}")

        maps_list = np.loadtxt(bsm_maps_list, dtype=str)
        rgw_list  = np.loadtxt(bsm_rgw_list,  dtype=str)

        # Combining the full list of smoothed scans into the smoothed science
        final_out = red_dir / f"{SOURCE}_IP_corrected_bsm{BSM}.fits"
        if final_out.exists():
            print(f"  {final_out.name} already exists, skipping combination")
        else:
            run_combination(base, red_dir, red_dir, bsm_maps_list, bsm_rgw_list, final_out)
            if final_out.exists():
                print(f"  Saved smoothed science map : {final_out}")
            else:
                print(f"  Combination failed: {final_out.name} was not created")

    # Stage 1: bootstrapping
    if "bootstrap" in active:
        stage_bootstrap(maps_list, rgw_list, base, combined_dir, N_BOOTSTRAP, N_CORES)

    # Stage 2: calculate covariance map
    if "covariance" in active:
        stage_covariance(combined_dir, base, BSM)

    print("\nDone.")

if __name__ == "__main__":
    main()