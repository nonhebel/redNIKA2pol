"""
source_create_products.py
--------------------------
Assembles final data products for a source by collecting reduction
outputs, smoothing for display, and computing error maps.

Workflow
--------
1. Symlinks the final-iteration uncorrected maps (from ``01_general/red``)
   and the IP-corrected combined map (from ``03_corr/red``), plus their
   associated rgw weight and RMS polygons (.pol), into
   ``05_products/``.
2. Smooths each map (corrected, uncorrected, rgw) to each value in
   ``BSM_VALS`` using ``bsm.piic``, for vector-plot-ready versions.
3. Computes RMS maps for the corrected/uncorrected maps at native and
   smoothed resolutions using ``creaRMS_IQU.piic``.
4. Overwrites the Q and U RMS planes of the IP-corrected RMS maps with
   standard deviations derived from the bootstrap covariance maps
   (``04_bootstrap/cov_map_bsm{BSM}.fits``), since bootstrapping gives a
   more reliable noise estimate than the pipeline RMS calculation for
   these terms.

Usage
-----
Set ``SOURCE`` and ``BSM_VALS`` at the top of the file, then run::

    python source_create_products.py
"""

import os
import subprocess
from pathlib import Path
import glob
from astropy.io import fits
import numpy as np

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------
SOURCE    = "???"
BSM_VALS = [5]     # don't need to put in 1 here for unsmoothed, just for additional smoothed maps

# ---------------------------------------------------------------------------
# Derived paths
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

base       = repo_root / "reductions" / SOURCE / "05_products"
uncorr_dir = repo_root / "reductions" / SOURCE / "01_general" / "red"
corr_dir   = repo_root / "reductions" / SOURCE / "03_corr" / "red"
bootstrap_dir = repo_root / "reductions" / SOURCE / "04_bootstrap" 

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
    print("RUN:", command)
    subprocess.run(command, shell=True, executable="/bin/bash", check=True)

def run_create_rms(map_path, rgw_path, I_pol_path, Q_pol_path, U_pol_path):
    """
    Execute the creaRMS_IQU.piic script to write the RMS map for given science 
    map.
    """
    command = (
        f"source ~/.bashrc; gagpiic; piic -nl @ creaRMS_IQU "
        f"{map_path.name} {rgw_path.name} {I_pol_path.name} {Q_pol_path.name} {U_pol_path.name} ; "
        f"rm PIIC*"
    )
    print("RUN:", command)
    subprocess.run(command, shell=True, executable="/bin/bash", check=True)

# ---------------------------------------------------------------------------
# Updating Q and U RMS terms with bootstrapped values
# ---------------------------------------------------------------------------

def overwrite_QU_rms(rms_path, cov_map_path):
    """
    Replace Q and U RMS planes in an IP-corrected RMS FITS file with
    the bootstrap-derived standard deviations.
    """
    with fits.open(cov_map_path) as hdul:
        cov_map = hdul[0].data  # shape: (ny, nx, 3, 3)

    sigma_Q = np.sqrt(cov_map[:, :, 1, 1])
    sigma_U = np.sqrt(cov_map[:, :, 2, 2])

    with fits.open(rms_path, mode='update') as hdul:
        hdul[0].data[1] = sigma_Q
        hdul[0].data[2] = sigma_U
        hdul.flush()

    print(f"Updated Q/U RMS planes in: {rms_path.name}")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

# Finding final iteration of reduction
tot_iter = int(glob.glob("*MP*n*.fits", root_dir=uncorr_dir)[0].split("i")[1].split("n")[0])
print(f"Using iteration: {tot_iter}")

# Symlink outputs into products directory
files = [
    (uncorr_dir, f"*MP*n{tot_iter}.fits",        f"{SOURCE}_uncorrected.fits"),
    (uncorr_dir, f"*MP*n{tot_iter}rgw.fits",     f"{SOURCE}_rgw.fits"),
    (uncorr_dir, f"*MP*n{tot_iter}_I.pol",       f"{SOURCE}_I.pol"),
    (uncorr_dir, f"*MP*n{tot_iter}_Q.pol",       f"{SOURCE}_Q.pol"),
    (uncorr_dir, f"*MP*n{tot_iter}_U.pol",       f"{SOURCE}_U.pol"),
    (corr_dir,   f"{SOURCE}_IP_corrected.fits",  f"{SOURCE}_IP_corrected.fits"),
]

for src_dir, src_name, dst_name in files:
    src = next(src_dir.glob(src_name))
    dst = base / dst_name
    if not dst.exists() and not dst.is_symlink():
        os.symlink(src, dst)

# Smooth maps for vector plotting
for bsm in BSM_VALS:
    for map_path in [
        base / f"{SOURCE}_IP_corrected.fits",
        base / f"{SOURCE}_uncorrected.fits",
        base / f"{SOURCE}_rgw.fits",
    ]:
        run_bsm('.', '.', map_path, bsm)

# Creating RMS maps
for map_path, rgw_path in [
    (base / f"{SOURCE}_IP_corrected.fits",  base / f"{SOURCE}_rgw.fits"),
    (base / f"{SOURCE}_uncorrected.fits",   base / f"{SOURCE}_rgw.fits"),
    *[(base / f"{SOURCE}_IP_corrected_bsm{bsm}.fits", base / f"{SOURCE}_rgw_bsm{bsm}.fits") for bsm in BSM_VALS],
    *[(base / f"{SOURCE}_uncorrected_bsm{bsm}.fits",  base / f"{SOURCE}_rgw_bsm{bsm}.fits") for bsm in BSM_VALS],
]:
    run_create_rms(
        map_path, rgw_path,
        base / f"{SOURCE}_I.pol",
        base / f"{SOURCE}_Q.pol",
        base / f"{SOURCE}_U.pol",
    )

# Replacing the Q and U entries for the IP corrected case with results from bootstrapping
overwrite_QU_rms(
    rms_path     = base / f"{SOURCE}_IP_corrected_rms.fits",
    cov_map_path = bootstrap_dir / "cov_map_bsm1.fits"
)
for bsm in BSM_VALS:
    overwrite_QU_rms(
        rms_path     = base / f"{SOURCE}_IP_corrected_bsm{bsm}_rms.fits",
        cov_map_path = bootstrap_dir / f"cov_map_bsm{bsm}.fits"
)
