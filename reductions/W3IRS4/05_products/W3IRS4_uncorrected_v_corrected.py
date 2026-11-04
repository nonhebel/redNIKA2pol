"""
source_uncorrected_v_corrected.py
-------------------------------------
Produces a side-by-side comparison of magnetic field vector maps
(uncorrected vs. IP-corrected) overlaid on Stokes I intensity, for
visually assessing the effect of instrumental polarisation correction.

Workflow
--------
Loads I/Q/U maps and their RMS maps (from ``05_products``) for both the
uncorrected and IP-corrected reductions, at two independently chosen
smoothing levels:

- ``BSM_BG``  : smoothing of the background, i.e. the Stokes I image and
  SNR contours. This also sets the pixel grid/WCS of the plot.
- ``BSM_VEC`` : smoothing of the magnetic field vectors. The polarisation
  quantities, significance masking and vector angles all come from this
  cube. Vectors are then placed onto the background pixel grid with WCS.

A value of 1 means native resolution (no smoothing).

For the vector cube, computes debiased polarised intensity, polarisation
fraction/angle and magnetic field vectors (polarisation angle rotated by
90°). Pixels below the ``SIGMA``-level polarisation significance
threshold are masked.

Produces a two-panel figure: Stokes I shown as a log-scaled image with
SNR contours, overlaid with the magnetic field vectors. Panels are
zoomed to the region defined by ``ZOOM_X``/``ZOOM_Y`` and labelled
"Uncorrected" / "IP corrected".

The Stokes I colour scale is set by the ``VMIN_PCT``/``VMAX_PCT``
percentiles. These are computed from the finite, positive background
pixels inside the zoom region of both panels together, so the two
panels share one colour scale.

Output is saved to
``figures/{SOURCE}/{SOURCE}_uncorr_corr_vector_map_bg{BSM_BG}_vec{BSM_VEC}.png``.

Usage
-----
Set ``SOURCE``, ``BSM_BG``, ``BSM_VEC`` and plot parameters (colour
percentiles, contour levels, significance, zoom region, vector scale) at
the top of the file, then run::

    python source_uncorrected_v_corrected.py
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
BSM_BG       = 1            # smoothing for Stokes I background + SNR contours (1 = native)
BSM_VEC      = 3            # smoothing for B-field vectors (1 = native)
CMAP         = "inferno"
VMIN_PCT     = 1            # lower percentile for Stokes I colour scale
VMAX_PCT     = 99.5         # upper percentile for Stokes I colour scale
SNR_LEVELS   = [5, 20, 50, 100, 500, 800, 1000, 2000, 5000, 10000]
SIGMA        = 3            # polarisation significance threshold (applied to vectors)
ZOOM_X       = (0.2, 0.75)  # fractional x limits (of background grid)
ZOOM_Y       = (0.2, 0.75)  # fractional y limits (of background grid)
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

base    = repo_root / "reductions" / SOURCE / "05_products"
fig_dir = repo_root / "figures" / SOURCE


def cube_paths(kind, bsm):
    """Return (cube_path, rms_path) for ``kind`` ('uncorrected' or
    'IP_corrected') at smoothing ``bsm`` (1 = native)."""
    suffix = "" if bsm == 1 else f"_bsm{bsm}"
    stem = f"{SOURCE}_{kind}{suffix}"
    return base / f"{stem}.fits", base / f"{stem}_rms.fits"

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
# Helpers
# ---------------------------------------------------------------------------

def percentile_limits(cubes, x_slice, y_slice, lo=VMIN_PCT, hi=VMAX_PCT):
    """
    Colour limits for Stokes I from percentiles of the zoomed region.

    Pools pixels from all ``cubes`` so the panels share one colour scale.
    Only finite, positive pixels are used, since LogNorm needs vmin > 0.
    """
    vals = np.concatenate([c.I_map[y_slice, x_slice].ravel() for c in cubes])
    vals = vals[np.isfinite(vals) & (vals > 0)]
    return np.percentile(vals, [lo, hi])

# ---------------------------------------------------------------------------
# Plot function
# ---------------------------------------------------------------------------

def plot_vector_map(fig, ax, bg_cube, vec_cube, vmin, vmax):
    """
    Plot Stokes I with SNR contours from ``bg_cube``, overlaid with
    magnetic field vectors from ``vec_cube``.

    Parameters
    ----------
    bg_cube : Cube
        Background cube (smoothing BSM_BG). Used for the Stokes I image,
        the SNR contours and the plot's pixel grid.
    vec_cube : Cube
        Vector cube (smoothing BSM_VEC). Used only for the B-field vectors
        and their significance masking. Must have had calc_polarisation()
        called.
    vmin, vmax : float
        Colour limits for the Stokes I image (LogNorm).
    """

    # Stokes I intensity
    im = ax.imshow(bg_cube.I_map, norm=LogNorm(vmin=vmin, vmax=vmax), cmap=cmap)
    cbar = fig.colorbar(im, ax=ax, location='top', pad=0.0)
    cbar.set_label('mJy/beam', fontsize=14)

    # SNR contours
    I_snr_map = bg_cube.I_map / bg_cube.I_error_map
    ny, nx = I_snr_map.shape
    ax.contour(
        np.arange(nx), np.arange(ny), I_snr_map,
        levels=SNR_LEVELS,
        alpha=0.5, zorder=10, cmap='Accent', linewidths=0.8,
    )

    # Vector pixel positions projected onto the background grid
    ny_vec, nx_vec = vec_cube.I_map.shape
    y_vec, x_vec = np.mgrid[0:ny_vec, 0:nx_vec]
    x_bg, y_bg = bg_cube.wcs.world_to_pixel(vec_cube.wcs.pixel_to_world(x_vec, y_vec))

    # Magnetic field vectors
    ax.quiver(
        x_bg, y_bg,
        vec_cube.x_mag_masked, vec_cube.y_mag_masked,
        headwidth=0, headlength=0, headaxislength=0,
        color='black', zorder=10, scale=VECTOR_SCALE, pivot='middle',
    )

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():

    # Load cubes
    print(f"Loading cubes (background BSM{BSM_BG}, vectors BSM{BSM_VEC})...")
    uncorr_bg  = Cube(*cube_paths("uncorrected",  BSM_BG))
    uncorr_vec = Cube(*cube_paths("uncorrected",  BSM_VEC))
    corr_bg    = Cube(*cube_paths("IP_corrected", BSM_BG))
    corr_vec   = Cube(*cube_paths("IP_corrected", BSM_VEC))

    # Polarisation is only needed for the vector cubes
    print(f"Computing polarisation for vector cubes (sigma={SIGMA})...")
    for cube in (uncorr_vec, corr_vec):
        cube.calc_polarisation(sigma=SIGMA)

    # Zoom region (in background pixels)
    ny_full, nx_full = uncorr_bg.I_map.shape
    x0, x1 = int(nx_full * ZOOM_X[0]), int(nx_full * ZOOM_X[1])
    y0, y1 = int(ny_full * ZOOM_Y[0]), int(ny_full * ZOOM_Y[1])

    # Shared colour limits from percentiles of the zoomed background
    vmin, vmax = percentile_limits((uncorr_bg, corr_bg), slice(x0, x1), slice(y0, y1))
    print(f"Stokes I colour limits: {vmin:.3g} – {vmax:.3g} "
          f"({VMIN_PCT}th–{VMAX_PCT}th percentile)")

    # Plot
    print("Plotting...")
    fig, (ax_uncorr, ax_corr) = plt.subplots(
        1, 2, figsize=(15, 10),
        subplot_kw={'projection': uncorr_bg.wcs},
    )

    plot_vector_map(fig, ax_uncorr, uncorr_bg, uncorr_vec, vmin, vmax)
    plot_vector_map(fig, ax_corr,   corr_bg,   corr_vec,   vmin, vmax)

    ax_uncorr.annotate('Uncorrected', (0.05, 0.9), xycoords='axes fraction',
                       ha='left', zorder=20, fontsize=20, color='white')
    ax_corr.annotate('IP corrected', (0.05, 0.9), xycoords='axes fraction',
                     ha='left', zorder=20, fontsize=20, color='white')

    for ax in (ax_uncorr, ax_corr):
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        ax.set_xlabel('RA', fontsize=14)
        ax.set_ylabel('Dec', fontsize=14)
        ax.annotate(SOURCE, (0.95, 0.9), xycoords='axes fraction',
                    ha='right', fontsize=20, color='white', zorder=20)

    fig_dir.mkdir(parents=True, exist_ok=True)
    outpath = fig_dir / f"{SOURCE}_uncorr_corr_vector_map_bg{BSM_BG}_vec{BSM_VEC}.png"
    plt.savefig(outpath, bbox_inches='tight', dpi=300)
    print(f"Saved: {outpath}")

if __name__ == '__main__':
    main()
