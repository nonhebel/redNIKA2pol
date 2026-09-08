"""
source_uncorrected_v_corrected.py
-------------------------------------
Produces a side-by-side comparison of magnetic field vector maps
(uncorrected vs. IP-corrected) overlaid on Stokes I intensity, for
visually assessing the effect of instrumental polarisation correction.

Workflow
--------
Loads the native-resolution and BSM-smoothed I/Q/U maps and their
associated RMS maps (from ``05_products``) for both the uncorrected and
IP-corrected reductions. For each, computes debiased polarised
intensity, polarisation fraction/angle, and magnetic field vectors
(polarisation angle rotated by 90°), masking pixels below the
``SIGMA``-level polarisation significance thresholds.

Produces a two-panel figure: Stokes I shown as a log-scaled image with
SNR contours, overlaid with magnetic field vectors drawn from the
BSM-smoothed cube and projected onto the native-resolution WCS. Panels
are zoomed to the region defined by ``ZOOM_X``/``ZOOM_Y`` and labelled
"Uncorrected" / "IP corrected".

Output is saved to
``figures/{SOURCE}/{SOURCE}_uncorr_corr_vector_map_bsm{BSM_VAL}.png``.

Usage
-----
Set ``SOURCE``, ``BSM_VAL``, and plot parameters (contour levels,
significance, zoom region, vector scale) at the top of the file, then
run::

    python plot_vector_map.py
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.wcs import WCS

# ---------------------------------------------------------------------------
# Top-level configuration
# ---------------------------------------------------------------------------
SOURCE       = "W3IRS4"
BSM_VAL      = 5
CMAP         = "inferno"
VMIN         = 5
VMAX         = 5000
SNR_LEVELS   = [5, 20, 50, 100, 500, 800, 1000, 2000, 5000, 10000]
SIGMA        = 3            # polarisation significance
ZOOM_X       = (0.2, 0.75)  # fractional x limits
ZOOM_Y       = (0.2, 0.75)  # fractional y limits
VECTOR_SCALE = 40
  
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

base   = repo_root / "reductions" / SOURCE / "05_products"
fig_dir = repo_root / "figures" / SOURCE 

bsm_suffix = "" if BSM_VAL == 1 else f"_bsm{BSM_VAL}"

paths = {
    'uncorr':          base / f"{SOURCE}_uncorrected.fits",
    'uncorr_rms':      base / f"{SOURCE}_uncorrected_rms.fits",
    'uncorr_bsm':      base / f"{SOURCE}_uncorrected{bsm_suffix}.fits",
    'uncorr_bsm_rms':  base / f"{SOURCE}_uncorrected{bsm_suffix}_rms.fits",
    'corr':            base / f"{SOURCE}_IP_corrected.fits",
    'corr_rms':        base / f"{SOURCE}_IP_corrected_rms.fits",
    'corr_bsm':        base / f"{SOURCE}_IP_corrected{bsm_suffix}.fits",
    'corr_bsm_rms':    base / f"{SOURCE}_IP_corrected{bsm_suffix}_rms.fits",
    }

# ---------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------
try:
    plt.style.use("standard")
except OSError:
    pass

cmap = plt.colormaps[CMAP].copy()
cmap.set_bad(cmap(0))

# ---------------------------------------------------------------------------
# Class to store cube info
# ---------------------------------------------------------------------------

class Cube:
    
    def __init__(self, path, rms_path):

        self.path = path
            
        # open main cube
        with fits.open(self.path) as hdul:
            hdu = hdul[0]
            self.header = hdu.header
            self.wcs = WCS(self.header).celestial
            self.I_map, self.Q_map, self.U_map = hdu.data

        # open rms cube
        with fits.open(rms_path) as hdul:
                hdu = hdul[0]
                self.I_error_map, self.Q_error_map, self.U_error_map = hdu.data

    def calc_polarisation(self, sigma=None, mask=None):

        if (sigma is None) == (mask is None):
            raise ValueError("Provide exactly one of 'sigma' or 'mask'")

        # Debiased polarised intensity
        self.Ipol_error_map = np.sqrt(self.Q_map**2*self.Q_error_map**2 + self.U_map**2*self.U_error_map**2) \
                        / np.sqrt(self.Q_map**2 + self.U_map**2)
        self.Ipol_map = np.sqrt(self.Q_map**2 + self.U_map**2 - self.Ipol_error_map**2)

        # Polarisation fraction 
        self.pol_frac_map = self.Ipol_map/self.I_map
        self.pol_frac_error_map = self.pol_frac_map * np.sqrt((self.Ipol_error_map/self.Ipol_map)**2 \
                            + (self.I_error_map/self.I_map)**2) 

        # Polarisation angle
        self.pol_angle_map = 0.5*np.arctan2(self.U_map, self.Q_map)
        self.pol_angle_error_map = 0.5*1/(self.Q_map**2 + self.U_map**2)*np.sqrt(self.Q_map**2*self.U_error_map**2 \
                                        + self.U_map**2*self.Q_error_map**2)
        
        # Polarisation vectors 
        x_pol = -np.sin(self.pol_angle_map)
        y_pol = np.cos(self.pol_angle_map)

        # Rotating by 90 degrees for magnetic field vectors
        self.x_mag = y_pol
        self.y_mag = -x_pol 

        # Masking by Ipol rms
        if mask is None:
            self.mask = (self.Ipol_map > sigma*self.Ipol_error_map) & (self.I_map > 5*self.I_error_map)
        else:
            self.mask = mask

        # Masking polarisation properties
        self.pol_frac_map_masked = np.ma.masked_array(self.pol_frac_map, mask=~self.mask)
        self.pol_angle_map_masked = np.ma.masked_array(self.pol_angle_map, mask=~self.mask)
        self.x_mag_masked = np.ma.masked_array(self.x_mag, mask=~self.mask)
        self.y_mag_masked = np.ma.masked_array(self.y_mag, mask=~self.mask)

# ---------------------------------------------------------------------------
# Plot function
# ---------------------------------------------------------------------------

def plot_vector_map(fig, ax, cube, bsm_cube):
    
    # Stokes I intensity
    im = ax.imshow(cube.I_map, norm=LogNorm(vmin=VMIN, vmax=VMAX), cmap=cmap)
    cbar = fig.colorbar(im, ax=ax, location='top', pad=0.0)
    cbar.set_label('mJy/beam', fontsize=14)
    
    # SNR contours
    I_snr_map = cube.I_map / cube.I_error_map
    ny, nx = I_snr_map.shape
    ax.contour(
        np.arange(nx), np.arange(ny), I_snr_map,
        levels=SNR_LEVELS,
        alpha=0.5, zorder=10, cmap='Accent', linewidths=0.8,
    )

    # Vectors from BSM cube projected onto unsmoothed cube
    ny_bsm, nx_bsm = bsm_cube.I_map.shape
    y_bsm, x_bsm = np.mgrid[0:ny_bsm, 0:nx_bsm]
    x_unsm, y_unsm = cube.wcs.world_to_pixel(bsm_cube.wcs.pixel_to_world(x_bsm, y_bsm))

    # Magnetic field vectors
    ax.quiver(
        x_unsm, y_unsm,
        bsm_cube.x_mag_masked, bsm_cube.y_mag_masked,
        headwidth=0, headlength=0, headaxislength=0,
        color='black', zorder=10, scale=VECTOR_SCALE, pivot='middle',
    )

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    cmap = plt.cm.get_cmap(CMAP).copy()
    cmap.set_bad(cmap(0))

    # Load cubes
    print("Loading cubes...")
    uncorrCube     = Cube(paths['uncorr'],     paths['uncorr_rms'])
    uncorrCube_bsm = Cube(paths['uncorr_bsm'], paths['uncorr_bsm_rms'])
    corrCube       = Cube(paths['corr'],       paths['corr_rms'])
    corrCube_bsm   = Cube(paths['corr_bsm'],   paths['corr_bsm_rms'])

    # Compute polarisation
    print(f"Computing polarisation quantities (sigma={SIGMA})...")
    for cube in (uncorrCube, uncorrCube_bsm, corrCube, corrCube_bsm):
        cube.calc_polarisation(sigma=SIGMA)

    # Plot
    print("Plotting...")
    fig, (ax_uncorr, ax_corr) = plt.subplots(
        1, 2, figsize=(15, 10),
        subplot_kw={'projection': uncorrCube.wcs},
    )

    plot_vector_map(fig, ax_uncorr, uncorrCube,  uncorrCube_bsm)
    plot_vector_map(fig, ax_corr,   corrCube,    corrCube_bsm)

    ax_uncorr.annotate('Uncorrected', (0.05,0.9), xycoords='axes fraction', 
                       ha='left', zorder=20, fontsize=20, color='white')
    ax_corr.annotate('IP corrected', (0.05,0.9), xycoords='axes fraction', 
                     ha='left', zorder=20, fontsize=20, color='white')

    ny_full, nx_full = uncorrCube.I_map.shape
    for ax in (ax_uncorr, ax_corr):
        ax.set_xlim(nx_full * ZOOM_X[0],  nx_full * ZOOM_X[1])
        ax.set_ylim(ny_full * ZOOM_Y[0],  ny_full * ZOOM_Y[1])
        ax.set_xlabel('RA', fontsize=14)
        ax.set_ylabel('Dec', fontsize=14)
        ax.annotate(SOURCE, (0.95, 0.9), xycoords='axes fraction',
                    ha='right', fontsize=20, color='white', zorder=20)

    outpath = fig_dir / f"{SOURCE}_uncorr_corr_vector_map_bsm{BSM_VAL}.png"
    plt.savefig(outpath, bbox_inches='tight', dpi=300)
    print(f"Saved: {outpath}")

if __name__ == '__main__':
    main()
