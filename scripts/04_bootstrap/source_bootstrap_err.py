"""
bootstrap_parllel.py
=====================
Bootstrap uncertainty of across science maps. 

Overview
--------
Starting from a set of individually reduced and corrected cubes 
(Stokes I, Q, U), this pipeline:

  1. Optionally smooths each scan to a target BSM using PIIC. 
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

  python bootstrap_parallel.py 

Pass STAGES = {"all"}`` to run everything, or a subset such as
{"bootstrap"} to run a single stage. 

"""

import sys
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
SOURCE          = ???                                                               # Name of target
REPO_ROOT       = ???                                                               # Location of repo
BSM             = 3                                                                 # Desired smoothing (BSM in PIIC)
N_BOOTSTRAP     = 5000                                                              # Number of bootstrap realisations
N_CORES         = 50                                                                # Number of cores to run on
STAGES          = {"all"}                                                           # Pipeline stages to execute

# ---------------------------------------------------------------------------
# Derived directory paths
# ---------------------------------------------------------------------------

base                = REPO_ROOT / "reductions" / SOURCE / "04_bootstrap"
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
    print("RUN:", command)
    subprocess.run(command, shell=True, executable="/bin/bash", cwd=base, check=True)

def run_combination(fits_dir, maps_list, rgw_list, outpath): 
    """
    Execute the combine_maps.piic script to combine a set of maps in list. 
    """
    command = (
        f"source ~/.bashrc; gagpiic; piic -nl @ ../combine_maps "
        # f"{str(outDir)} {str(inDir)} {maps_list.name} {rgw_list.name} {outpath.name}; "
        f"../red . {maps_list.name} {rgw_list.name} {outpath.name} ; "
        f"rm PIIC*"
    )
    print("RUN:", command)
    subprocess.run(command, shell=True, executable="/bin/bash", cwd=fits_dir, check=True)

# ---------------------------------------------------------------------------
# Bootstrapping procedure
# ---------------------------------------------------------------------------

def run_bootstrap(args):
    """
    Run a single bootstrap realisation.
    """
    
    i, base, maps_list, rgw_list, combined_dir = args

    # Skip realisations that are already completed to allow safe restart
    final_combined = combined_dir / f"combined_{i}.fits"
    if final_combined.exists():
        print(f"[{i}] already exists, skipping")
        return i
    
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

    # Combine bootstrapped maps
    run_combination(run_dir, maps_list_i, rgw_list_i, combined_i)

    shutil.move(combined_i, final_combined)
    shutil.rmtree(run_dir)

    return i

# ---------------------------------------------------------------------------
# Pipeline stages
# ---------------------------------------------------------------------------

def stage_bootstrap(maps_list, rgw_list, base, combined_dir, n_bootstrap, n_cores):
    """
    Stage 1 - generate bootstrap realisations in parallel.
    """
    print("\n=== Stage 1: bootstrapping over individual corrected maps ===")
    
    args = [
        (i, base, maps_list, rgw_list, combined_dir)
        for i in range(n_bootstrap)
    ]

    with Pool(processes=n_cores) as p:
        for _ in tqdm.tqdm(
            p.imap_unordered(run_bootstrap, args), 
            total=n_bootstrap, 
            desc="Bootstrap realisations"
        ):
            pass

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
    """Smooth scans if need be and run selected stages."""
    print("Source:", SOURCE)
    active = {"bootstrap", "covariance"} if "all" in STAGES else STAGES
    print("Active stages:", active)

    maps_list = np.loadtxt(base / "all_maps.LIST", dtype=str)
    rgw_list = np.loadtxt(base / "all_rgw.LIST", dtype=str)

    # Performing smoothing of individual scans
    if BSM != 1:
        bsm_maps_list = base / f"all_maps_bsm{BSM}.LIST"
        bsm_rgw_list = base / f"all_rgw_bsm{BSM}.LIST"

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

    # Stage 1: bootstrapping
    if "bootstrap" in active:
        stage_bootstrap(maps_list, rgw_list, base, combined_dir, N_BOOTSTRAP, N_CORES)

    # Stage 2: calculate covariance map
    if "covariance" in active:
        stage_covariance(combined_dir, base, BSM)

    print("\nDone.")

if __name__ == "__main__":
    main()
