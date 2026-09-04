# Reduction of NIKA2pol observations

Welcome to the GitHub page outlining the data reduction pipeline used by the DYNAMAG large programme for the reduction of NIKA2 polarimetry data, including the correction of instrumental polarisation. Scripts are shared here to help others similarly reducing NIKA2pol data - if you use them please acknowledge this page! While the scripts provided here are intended to be general, in places choices have been made specific to the dataset in question that may not be the best fit for all data.

The scripts provided here are a mix of Python and PIIC scripts. The Pointing and Imaging in Continuum (PIIC) software is the data reduction pipeline provided by IRAM. For information on PIIC, please see the very helpful handbook (https://www.iram.fr/~g\
ildas/dist/piic.pdf), which outlines the installation and usage of PIIC. The pipeline introduced here makes use of many of the core PIIC scripts, with additional scripts to ease especially the process of correcting for instrumental polarisation. The main novelties are source_ip_corr_MP.py and source_bootstrap_err.py scripts.  

The pipeline is not perfect and likely improvements could be made - feel free to get in touch with suggestions or questions. 

## Initial setup

To get started, first clone this repository with:

git clone https://github.com/<username>/redNIKA2pol.git

After cloning, the repo has the following structure:

      redNIKA2pol/
      ├── README.md
      ├── setup.sh
      ├── raw/
      │   ├── science/
      │   └── calibrators/
      ├── scripts/
      │   ├── 01_general/
      │   ├── 02_individual/
      │   ├── 03_corr/
      │   ├── 04_bootstrap/
      │   ├── 05_products/
      ├── reductions/
      │   └── <target>/
      │       ├── 01_general/
      │       ├── 02_individual/
      │       ├── 03_corr/
      │       ├── 04_bootstrap/
      │       └── 05_products/
      └── figures/

Place your raw NIKA2 scan files under `raw/science`. Download the calibrators directory from the PIIC webpage and place under `raw/calibrators`.

The repo is currently set up for the target W3IRS4, included as an example. To generate	the same file tree for another target, open up the setup.sh script. Enter your	target name and	run the	script with ./setup.sh. This will create the same structure as above for your given target, with the necessary scripts for each step of the pipeline where they need to be. Scripts should be execute within the directories in which they are found. 

## Reducing your data

Now that everything is in the correct place, we can carry out the data reduction. This procedure is separated into 5 different steps, each with a separate direction in recuctions/<targets>.

### 01_general

The first step in the reduction process produces a science map of the target from the raw detector timelines uncorrected for the effect of instrumental polarisation (IP). The reduction is an iterative process that works repeatedly to remove sky signal and instrument instabilities from the source signal. The procecure followed is outlined in the PIIC handbook and makes use of the scripts set up for multi-core processing, source_a13_MP.piic and source_setMPpar.piic. A list of the scans to be reduced must be made, with the file name source_a13.LIST, which is then separated into multiple lists, one for each core (see PIIC 5.2). Within the source_a13_MP.piic script, edit the source name, the number of procedures (cores to run on), and the number of iterations. Initially set the number of interations to 0. Within the source_setMPpar.piic, you can select the smSNRpar and blOrderOrig parameters. These will depend on the specfics of your target (see PIIC 3.2.6 and 4.14.4) - here we set smSNRpar=0 and blOrderOrig=4. Also within source_setMPpar.piic, source and base range polygons can be set, as well as a zeroing radius (see PIIC 6.5.1 for information on how to define the polygons). For the initial zero-iteration reduction, leave the polygons commented out. To run the reduction in batch mode, use the command:

nohup piic @ source_a13_MP.piic > source_a13_MP.log2 2>&1 &

The output science map will give a first idea of the source structure, from which to define the source and base range polygons (source.pol and base.pol). The scripts should then be edited to run the reduction again, this time with the polygon definitions and with many iterations. The number of iterations needed for convergence will change with source. For the DYNAMAG reductions, around 50 iterations are needed. Convergence of the reductions can be checked using the check_reduction_convergence.py script, which produces plots of the changing flux and noise with increasing iterations.

### 02_individual

The second step then produces reduced maps for each individual science scan, using the final source definition (sbSource) produced in 01_general. This is carried out using the script source_a13_save_ind.piic. The reduction should be run with the same parameters as above, i.e. the same smSNRpar and blOrderOrig, as well as the same source and base range polygons, and zeroing radius. By saving individual maps, we can then individually correct them for IP which will vary from scan-to-scan. As above, run the script using:

nohup piic @ source_a13_save_ind.piic > source_a13_save_ind.log2 2>&1 &

### 03_corr

The third step aims to correct the individual maps produced in 02_individual for the contribution to Stokes Q and U from leakage from Stokes I. We call this instrumental polarisation, or IP. The tricky part about the IP present in NIKA2pol observations is that it is sensitive to observing conditions, including the elevation of the target, offset from optimal focus, and weather conditions. A blanket correction to all scans therefore cannot be made, and instead each scan must be corrected individually. 

One way to determine the IP present in observations is to observe an unpolarised target, i.e. a planet, where all polarised signal can be attributed to leakage from Stokes I. The standard NIKA2pol observing procedure aims to observe a planet before and after science observations, which can then be used to correct for the IP by deconvolving this signal from the science map. There are, however, number of problems with this method. Often the elevation of the planet will differ from that of the science scan, leading to a poor characterisation of the IP pattern given the sensitivity of IP to elevation. Elevation differences > 10 degrees can lead to large residuals when the IP correction is performed, leading to artefacts in the final science maps. An even larger issue, however, is when no planets at all available to observe at the time observations. In this case, calibration measurements are instead taken on QSOs, which have a time-varying degree of polarisation and therefore need additional calibrator measurements to remove the intrinsic polarisation of the QSO. 

Here we don't limit ourselves to calibration measurements taken before and after the science scans, and instead match each science scan to an appropriate planet calibration scan via different means. This is carried out by the source_ip_corr_MP.py.

#### source_ip_corr_MP.py

The aim of source_ip_corr_MP.py is to correct each science scan with the 'best' calibrator scan from a catalogue of such scans (scans placed in `raw/calibrators`). This catalogue was produced by reducing all available planet calibration scans (polarised calib_1scans) from the nikas-24, nikaw-24, nikas-25, nikaw-25 and 084-25 projects. The catalogue will be expanded to include planet observations from other projects in the future. The observations cover a large range of observing conditions and thus IP patterns. 

source_ip_corr_MP.py works by finding candidate calibration scans from the catalogue of scans, for each science scan. It does this by filtering by elevation, selecting only calibrators that are within a user-set ELEV_THRESH from the science target. It then finds which calibrator produces the 'best' IP correction for each science scan and produces the final IP corrected science map. It does this through a number of stages:

- Stage 1) - CORRECT: Each science scan is corrected for IP with every identified candidate calibrator using corr_instr_pol.piic. Note that the version of corr_instr_pol.piic used here is *not* identical as the standard script provided in PIIC (i.e. in piic/pro); the script has been modified to work on individual scans rather than scan lists. corr_instr_pol.piic corrects for the IP of the science scan via deconvolving the candidate calibrator signal from the science signal, with appropriate scaling by the Stokes I component (see PIIC handbook). The procedure is carried out for both array 1 and array 3 individually - it is important that this is the case due to focus offsets between the arrays, leading to different IP patterns.
- Stage 2) - SELECT: The 'best' calibrator is selected for each science scan. Unlike matching the elevation of the science and calibration measurements, matching the remaining conditions which may influence the IP (e.g. the offset from the optimal focus at the time of observations, variations in the temperature of the telescope etc.) is non-trivial. Instead, the best calibrator is decided *a posteriori*, via a metric calculated for the IP-corrected science maps created using each candidate calibrator. The metric is chosen to be the standard deviation across the inner region of the corrected science map. The reasoning behind this is that the IP signal in Stokes Q and U is known to be cloverleaf pattern with positive and negative lobes, and thus the better the IP correction, the lower the standard deviation should be. The central region of the map is focussed on such that noise at the map-edges doesn't have an impact. Moreoever, the IP scales with Stokes I, and the centre of the map, where the bright Stokes I source is present, will be most affected by the IP correction. The metric is calculated for each candidate calibrator, with the best calibrator that which minimising it. A figure is produced of all of the possible IP corrected maps, with the best map highlighted.
- Stage 3) - REFINE: An *optional* additional step to attempt to improve the correction is to refine the calibrator selection by reducing the 'stack' standard deviation. The idea behind this step is that the best calibrators should lead to the lowest standard deviation *between* IP-corrected maps - as while the IP signal will vary from scan-to-scan, the true science signal will not, beyond the usual noise. To minimise the 'stack' standard deviation - the pixel-wise standard deviation of Q and U maps normalised by I in arrays 1 and 3 across all scans - a coordinate descent optimisation is performed: Starting from the initial best calibrator assignments from Stage 2), for each science scan, each candidate calibrator is trialled in turn, while holding the remainder of the corrected science maps the same. The best calibrator is updated if the stack standard deviation decreases. This procedure is performed iteratively until the standard deviation no longer decreases, i.e. convergence has been achieved. Two figures are produced; a figure showing the best corrected maps after Step 2), and a figure showing the best corrected maps after Step 3). Comparing the two figures, the variation between the science scans should be lower. *However*, whether this really leads to a 'better' IP correction is not guaranteed - just because the scans are now consistent, does not mean that are consistently correct. From testing, the optional REFINE step tends to swap a large amount of the best calibrators, with the resulting corrected scans indeed more similar, but often with suspicously strong polarisation and patterns reminiscent of IP. Conclusion: use with caution and check the diagnostic plots.
- Stage 4) - WRITE_LISTS: The best corrected maps are written to all_corr_maps.LIST and the associated rgw weights maps are written to all_rgw.LIST.
- Stage 5) - COMBINE: combine_maps.piic is used to co-add the corrected maps, weighted appropriately by the corresponding rgw maps, to produce a final IP corrected science map of the target.

### 04_bootstrap

The fourth step aims to estimate the uncertainty in Stokes Q and U, taking into account additional uncertainty introduced by imperfect IP correction. This is done by bootstrapping over the final IP-corrected individual scans, using the source_boostrap_err.py script. Bootstrapping is a means of estimating statistics of a dataset via randomly sampling the dataset with replacement: source_boostrap_err.py randomly samples with replacement the list of corrected scans produced in the 03_corr step, drawing a user-set N_BOOTSTRAP samples. The per-pixel 3x3 covariance matrix over Stokes I, Q and U is then calculated over the full set of bootstrap realisations. source_bootstrap_err.py is parallelised and is run over N_CORES cores. 

As noted above, a 'better' IP correction will lead to reduced variation in the Stokes Q and U signal between scans, tending to the limit of simply the noise in the case IP has perfectly been removed / is not present. The covariance matrix at each point in the map therefore gives us an idea of additional error, in excess of the typical RMS noise, due to residual signal from the IP. 

source_bootstrap_err.py can produce bootstrap covariance maps at the native map resolution or for spatially smoothed maps. bsm.piic applies a BSM x BSM pixel box average to the original cube, producing smooth maps with pixels BSM x BSM times larger than in the original.  Setting BSM > 1 causes the pipeline to first smooth each individual scan before bootstrapping over the smoothed scans. Using a BSM > 1 may be desirable in cases of low signal-to-noise. 

### 05_products

The final step produces science-ready data products. create_products.py symlinks the uncorrected and IP-corrected maps, along with the associated RGW map and the Stokes I, Q, U RMS polygons output by PIIC, into the 05_products directory. Spatially smoothed versions of the maps are produced for the BSM values set in BSM_VALS. RMS noise maps for each science map are then generated via creaRMS_IQU.piic, which rescales the RGW weight map by the RMS measured within the polygons.

For the IP-corrected maps, the PIIC-derived RMS maps are compared against the bootstrap covariance maps from 04_bootstrap. Away from the map centre — where Stokes I emission is weak and IP leakage is negligible — the two agree closely. Toward the centre, the bootstrap uncertainty in Stokes Q and U exceeds the PIIC noise estimate, reflecting residual IP that the correction has not fully removed. The Stokes Q and U planes of the RMS maps for the IP-corrected case are therefore replaced with the bootstrap-derived standard deviations, giving a more conservative and physically motivated noise estimate in regions of significant IP leakage.



 


