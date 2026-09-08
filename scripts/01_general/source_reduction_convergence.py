"""
source_reduction_convergence.py
--------------------------
Diagnostic plots for checking convergence of maps as a function of iteration
number.

Workflow
--------
Loads all per-iteration I/Q/U cubes for a given source from
``01_general/red``, and for each iteration:

- computes the RMS of Stokes I, Q and U within PIIC-defined RMS polygon
- computes the total Stokes I and polarised intensity over the full map

Two figures are produced:

1. RMS convergence — RMS I, Q, U vs. iteration number
2. Flux convergence — total I and total Ipol vs. iteration number

Each panel annotates the fractional change between the last two
iterations as a quick convergence check. Figures are saved to
``figures/{SOURCE}/{SOURCE}_rms_convergence.pdf`` and
``..._flux_convergence.pdf``.

Usage
-----
Set ``SOURCE`` at the top of the file, then run::

    python source_reduction_convergence.py
"""

from pathlib import Path
import glob

import numpy as np
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u

from regions import RectangleSkyRegion

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------
SOURCE    = "???"
# ---------------------------------------------------------------------------
try:
    plt.style.use("standard")
except OSError:
    pass

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

red_dir = repo_root / "reductions" / SOURCE / "01_general" / "red"
fig_dir = repo_root  / "figures" / SOURCE 
fig_dir.mkdir(parents=True, exist_ok=True)

# Infer total number of iterations from filename
tot_iter = int(glob.glob("*MP*n*.fits", root_dir=red_dir)[0].split("i")[1].split("n")[0])

class Cube:
    def __init__(self, path):
        self.path   = Path(path)
        self.hdu    = fits.open(self.path)[0]
        self.header = self.hdu.header
        self.wcs    = WCS(self.header).celestial
        self.I_map, self.Q_map, self.U_map = self.hdu.data
        self.Ipol_map = np.sqrt(self.Q_map**2 + self.U_map**2)

    def calc_rms(self, stokes):
        pol_path = self.path.with_name(self.path.stem + f"_{stokes}.pol")
        lb, rb, rt, lt = np.loadtxt(pol_path)
        corners   = np.array([lb, rb, rt, lt])
        dx        = np.abs(lb[0] - rb[0])
        dy        = np.abs(rt[1] - rb[1])
        cx_offset = corners[:, 0].mean()
        cy_offset = corners[:, 1].mean()
        cx        = self.header["CRVAL1"] + (cx_offset / 3600) / np.cos(np.deg2rad(self.header["CRVAL2"]))
        cy        = self.header["CRVAL2"] + (cy_offset / 3600)
        centre    = SkyCoord(cx, cy, unit="deg", frame="fk5", equinox="J2000")
        region_p  = RectangleSkyRegion(centre, dx * u.arcsec, dy * u.arcsec).to_pixel(self.wcs)
        data      = getattr(self, f"{stokes}_map")
        mask      = region_p.to_mask().to_image(data.shape).astype(bool)
        setattr(self, f"rms_{stokes}", np.nanstd(data[mask]))

metrics = {k: np.zeros(tot_iter) for k in ("rms_I", "rms_Q", "rms_U", "tot_I", "tot_Ipol")}

for i in range(tot_iter):
    path = red_dir / glob.glob(f"*MP*i{tot_iter}n{i}.fits", root_dir=red_dir)[0]
    cube = Cube(path)
    for stokes in ("I", "Q", "U"):
        cube.calc_rms(stokes)
    metrics["rms_I"][i]    = cube.rms_I
    metrics["rms_Q"][i]    = cube.rms_Q
    metrics["rms_U"][i]    = cube.rms_U
    metrics["tot_I"][i]    = np.nansum(cube.I_map)
    metrics["tot_Ipol"][i] = np.nansum(cube.Ipol_map)

iters = np.arange(tot_iter)

# Figure 1: RMS convergence
fig, axes = plt.subplots(1, 3, figsize=(16, 4))
rms_panels = [
    (axes[0], "rms_I", "RMS I (mJy/beam)"),
    (axes[1], "rms_Q", "RMS Q (mJy/beam)"),
    (axes[2], "rms_U", "RMS U (mJy/beam)"),
]
for ax, key, ylabel in rms_panels:
    ax.scatter(iters, metrics[key], s=20, marker="x")
    ax.set_xlabel("Iterations")
    ax.set_ylabel(ylabel)
    diff = (metrics[key][-1] - metrics[key][-2]) / metrics[key][-1] * 100
    ax.annotate(f"Δ = {diff:.2f}%", xy=(0.1, 0.9), xycoords="axes fraction")

plt.tight_layout()
fig.savefig(fig_dir / f"{SOURCE}_rms_convergence.png", bbox_inches="tight", dpi=300)
plt.close()

# Figure 2: Intensity
fig, axes = plt.subplots(1, 2, figsize=(8, 4))
tot_panels = [
    (axes[0], "tot_I",    "Total I (mJy/beam)"),
    (axes[1], "tot_Ipol", "Total Ipol (mJy/beam)"),
]

for ax, key, ylabel in tot_panels:
    ax.scatter(iters, metrics[key], s=20, marker="x")
    ax.set_xlabel("Iterations")
    ax.set_ylabel(ylabel)
    diff = (metrics[key][-1] - metrics[key][-2]) / metrics[key][-1] * 100
    ax.annotate(f"Δ = {diff:.2f}%", xy=(0.1, 0.9), xycoords="axes fraction")

plt.tight_layout()
fig.savefig(fig_dir / f"{SOURCE}_flux_convergence.png", bbox_inches="tight", dpi=300)
plt.close()

print(f"Saved {fig_dir}/{SOURCE}_rms_convergence.png and {fig_dir}/{SOURCE}_flux_convergence.png")
