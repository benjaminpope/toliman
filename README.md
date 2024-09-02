# toliman
Code related to the TOLIMAN astrometric mission concept

## Files

These are pretty much untouched from Ben's original work
* `bandwidth.ipynb` investigating wavelength stability of a star
* `poppy_toliman_test.ipynb` investigating PSF for Toliman rosette using `Poppy`
* `synthetic.ipynb` investigating wavelength stability for various band pass filters
* `synthetic-mcmc.ipynb` finding parameters from synthetic data using MCMC
* `poppy_toliman_doubly_diffractive.ipynb` fitting using least squares and Richardson–Lucy deconvolution (NOTE: some bugs - may need to revert to an earlier version to run properly)
* `diffraction.tex` Ben's write-up of the above investigations

These are new files added by Bryn
* `jupyter_setup.md` Bryn's notes on setting up everything to run Ben's existing Toliman notebooks
* `toliman_image_simulation.ipynb` Model images with noise (in progress, based upon Ben's code from above)
* `toliman_optics.py` Code factored out from Ben's notebooks for the rosette aperture
