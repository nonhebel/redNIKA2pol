"""
source_create_products.py
--------------------------
Assembles final data products for a source at a single smoothing (BSM) by
collecting reduction outputs, smoothing for display, and computing error maps.
 
Workflow
--------
1. Symlinks the final-iteration uncorrected maps (from ``01_general/red``)
   and the IP-corrected combined map (from ``03_corr/red``), plus their
   associated rgw weight and RMS polygons (.pol), into ``05_products/``.
2. If BSM > 1:
     - symlinks the smoothed IP-corrected map
       (``04_bootstrap/{SOURCE}_IP_corrected_bsm{BSM}.fits``), made by
       ``source_bootstrap_err.py`` by smoothing each scan and then combining,
       so that it is consistent with the bootstrap covariance map;
3. Computes RMS maps for the corrected and uncorrected maps at this BSM
   using ``creaRMS_IQU.piic``.
4. Overwrites the Q and U RMS planes of the IP-corrected RMS map with
   standard deviations derived from the bootstrap covariance map
   (``04_bootstrap/cov_map_bsm{BSM}.fits``), since bootstrapping gives a
   more reliable noise estimate than the pipeline RMS calculation for
   these terms.
 
Requires ``source_bootstrap_err.py`` to have been run with the same BSM first.
 
Usage
-----
Set ``SOURCE`` and ``BSM`` at the top of the file, then run::
 
    python source_create_products.py
 
Run once per BSM value (e.g. BSM = 1, BSM = 5).
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
SOURCE    = "???"       # Name of target
BSM       = 1           # Desired smoothing (BSM in PIIC, 1 for native resolution)
 
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

# File-name suffix for this smoothing: none at native resolution
sfx = "" if BSM == 1 else f"_bsm{BSM}"

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

print("Source:", SOURCE)
print("BSM:", BSM)

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

# Smoothed maps for vector plotting
if BSM != 1:
    # IP-corrected
    src = corr_dir / f"{SOURCE}_IP_corrected{sfx}.fits"
    dst = base / src.name
    if not dst.exists() and not dst.is_symlink():
        os.symlink(src, dst)

    # Uncorrected map and rgw
    for map_path in [
        base / f"{SOURCE}_uncorrected.fits",
        base / f"{SOURCE}_rgw.fits",
    ]:
        run_bsm('.', '.', map_path, BSM)

# Creating RMS maps at this BSM
for map_path in [
    base / f"{SOURCE}_IP_corrected{sfx}.fits",
    base / f"{SOURCE}_uncorrected{sfx}.fits",
]:
    run_create_rms(
        map_path,
        base / f"{SOURCE}_rgw{sfx}.fits",
        base / f"{SOURCE}_I.pol",
        base / f"{SOURCE}_Q.pol",
        base / f"{SOURCE}_U.pol",
    )

# Replacing the Q and U entries for the IP corrected case with results from bootstrapping
overwrite_QU_rms(
    rms_path     = base / f"{SOURCE}_IP_corrected{sfx}_rms.fits",
    cov_map_path = bootstrap_dir / f"cov_map_bsm{BSM}.fits",
)
 
print("\nDone.")