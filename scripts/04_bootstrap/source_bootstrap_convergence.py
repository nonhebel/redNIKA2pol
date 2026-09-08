"""
source_bootstrap_convergence.py
--------------------------
Diagnostic plot for checking convergence of bootstrap-resampled
polarisation maps.

Workflow
--------
Loads all ``combined_*.fits`` bootstrap realisation maps from
``04_bootstrap/combined_bsm{BSM}`` for a given source, then tracks the
running standard deviation of Stokes Q and U at two representative pixels
(map centre and an outer pixel) as a function of the number of bootstrap
realisations included. This is used to check whether enough bootstrap
realisations have been generated for the noise estimates to stabilise.

The resulting convergence plot (STD vs. number of realisations, for both
pixels and both Stokes parameters) is saved to
``figures/{SOURCE}/{SOURCE}_bootstrap_convergence_bsm{BSM}.pdf``.

Usage
-----
Set ``SOURCE`` and ``BSM`` at the top of the file, then run::

    python source bootstrap_convergence.py
"""

import glob
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from astropy.io import fits

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------
SOURCE          = "W3IRS4"                    # Name of target
BSM             = 5                           # Desired smoothing (BSM in PIIC)

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
combined_dir        = base / f"combined_bsm{BSM}"
fig_dir             = repo_root / "figures" / SOURCE

# ---------------------------------------------------------------------------
# Matplotlib style
# ---------------------------------------------------------------------------
try:
    plt.style.use("standard")
except OSError:
    pass

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
print("Running...")

fits_names = glob.glob("combined_*.fits", root_dir=combined_dir)
num_maps = max(int(name.split("_")[1].split(".fits")[0]) for name in fits_names)

boot_maps = []
for path in combined_dir.glob("combined_*.fits"):
    with fits.open(path, memmap=False) as hdul:
        data = hdul[0].data
    boot_maps.append(data)

boot_maps = np.array(boot_maps)  # (n_boot, 3, ny, nx)

ny, nx = boot_maps.shape[-2:]

# centre pixel and an outer pixel to compare convergence at
pixels = {
    "centre": (ny // 2, nx // 2),
    "outer":  (ny // 4, nx // 4),
}

indices = np.linspace(1, num_maps, 100).astype(int)

fig, axes = plt.subplots(1, 2, figsize=(8, 4), sharey=False)

for ax, (label, (y_pix, x_pix)) in zip(axes, pixels.items()):
    cov_i = np.zeros([len(indices), 3, 3])
    for i, index in enumerate(indices):
        samples = boot_maps[:index, :, y_pix, x_pix]
        mask = np.all(np.isfinite(samples), axis=1)
        if mask.sum() > 1:
            cov_i[i] = np.cov(samples[mask].T)

    ax.plot(indices, np.sqrt(cov_i[:, 1, 1]), color='tab:red', label='Stokes Q')
    ax.plot(indices, np.sqrt(cov_i[:, 2, 2]), color='tab:blue', label='Stokes U')
    ax.set_xlabel("Number of bootstrap realisations")
    ax.set_ylabel("STD at pixel")
    ax.set_title(f"{label} pixel (y={y_pix}, x={x_pix})")
    ax.legend()

fig.tight_layout()

out_path = fig_dir / f"{SOURCE}_bootstrap_convergence_bsm{BSM}.png"
fig.savefig(out_path, bbox_inches="tight", dpi=300)

print(f"Saved figured to {out_path}")