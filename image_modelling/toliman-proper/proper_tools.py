import numpy as np
import proper

def normalise_sampling(wavefronts, samplings, common_sampling, npsf):
    """Resample each wavefront to a common grid
    
    Parameters
    ----------
    wavefronts : list of PROPER WaveFront class objects
        The wavefronts to be resampled
        
    samplings: list of floats 
        wavefront samplings for each wavefront, in metres
    
    common_sampling : float
        Target wavefront sampling in metres
        
    npsf : int
        Conic constant, <-1 for hyperbolas, -1 for parabolas, between -1 and 0 
        for ellipses, 0 for spheres, and greater than 0 for oblate ellipsoids

        
    Returns
    -------
    out : numpy ndarray
        Returns stack of wavefronts with common dimensions and sampling.
    """
    n = len(wavefronts)
    out = np.zeros([n, npsf, npsf], dtype = np.float64)
    # Resample and weight
    for i in range(n):
        wf = wavefronts[i] # np.abs(wavefronts[i])**2
        mag = samplings[i] / common_sampling
        out[i,:,:] = proper.prop_magnify(wf, mag, npsf, CONSERVE = True)
    return out

def combine_psfs(psfs, weights):
    """Combine stack of PSFs into a single 2D array 
    
    Parameters
    ----------
    psfs : 3D numpy array
        The PSFs to be combined, as a stack of 2D arrays of same dimensions
        
    weights : float
        Target wavefront sampling in meters
        
    Returns
    -------
    out : numpy ndarray
        Returns single 2D image array
    """
    out = psfs[0,:,:] * weights[0]    
    for i in range(1,len(weights)):
        out += psfs[i,:,:] * weights[i]
    return out

from numpy.fft import fft2, ifft2

# Looks like proper.prop_pixellate is broken, so paste directly and hack it up
def fix_prop_pixellate(image_in, sampling_in, sampling_out, n_out = 0):
    """Integrate a sampled PSF onto detector pixels. 
    
    This routine takes as input a sampled PSF and integrates it over pixels of a 
    specified size. This is done by convolving the Fourier transform of the input 
    PSF with a sinc function representing the transfer function of an idealized 
    square pixel and transforming back. This result then represents the PSF 
    integrated onto detector-sized pixels with the same sampling as the PSF. 
    The result is interpolated to get detector-sized pixels at detector-pixel 
    spacing.
    
    Parameters
    ----------
    image_in : numpy ndarray
        2D floating image containing PSF
        
    sampling_in : float
        Sampling of image_in in meters/pixel
        
    sampling_out : float
        Size(=sampling) of detector pixels
        
    n_out : int
        Output image dimension (n_out by n_out)
        
    Returns
    -------
    new : numpy ndarray
        Returns image integrated over square detector pixels.
    """
    n_in = image_in.shape[0]
    
    n_in_half = int(n_in/2)
    
    # Compute pixel transfer function (MTF)
    psize = 0.5 * (sampling_out / sampling_in)
    constant = psize * np.pi
    mag = sampling_in / sampling_out
    
    arr = np.arange(n_in, dtype = np.float64) - n_in/2.
    # BJ
#    x = np.roll(arr, -n_in/2, 0) / (n_in/2.)
    x = np.roll(arr, -n_in_half, 0) / (n_in/2.)
    t = x * constant
    y = np.zeros(n_in, dtype = np.float64)
    y[1:] = np.sin(t[1:]) / t[1:]
    y[0] = 1.
    
    #pixel_mtf = np.tile(np.vstack(y), (1, n_in)) * y
    pixel_mtf = np.dot(y[:,np.newaxis], y[np.newaxis,:])
    
    # Convolve image with detector pixel
    image_mtf = fft2(np.roll(np.roll(image_in, -n_in_half, 0), -n_in_half, 1)) / image_in.size
    image_mtf *= pixel_mtf
    
    convolved_image = np.abs(ifft2(image_mtf) * image_mtf.size)
    image_mtf = 0
    convolved_image = np.roll(np.roll(convolved_image/mag**2, n_in_half, 0), n_in_half, 1)
    
    # Image is integrated over pixels but has original sampling; now, resample
    # pixel sampling
    if n_out != 0:
        n_out = int(np.fix(n_in * mag))
        
    new = proper.prop_magnify(convolved_image, mag, n_out)
    
    return new

def form_detector_image(prescription, sources, gridsize, detector_pitch, npixels, multi=True):
    source_psfs = []
    common_sampling = detector_pitch/2. # for Nyquist 
    npsf = npixels*2
    for source in sources:
        settings = source['settings']
        wavelengths = source['wavelengths']
        wl_weights = source['weights']

        if multi is True:
            (wavefronts, samplings) = proper.prop_run_multi(prescription, wavelengths, gridsize = gridsize, QUIET=True, PRINT_INTENSITY=False, PASSVALUE=settings)
            # prop_run_multi returns complex arrays, even when PSFs are intensity, so make real with abs
            psfs = normalise_sampling(np.abs(wavefronts), samplings, common_sampling, npsf)
        else:
            wavefronts = []
            samplings = []
            for wl in wavelengths:
                (wavefront, sampling) = proper.prop_run(prescription, wl, gridsize = gridsize, QUIET=True, PRINT_INTENSITY=False, PASSVALUE=settings)
                wavefronts.append(wavefront)
                samplings.append(sampling)
            psfs = normalise_sampling(wavefronts, samplings, common_sampling, npsf)
            
        source_psfs.append(combine_psfs(psfs, wl_weights))

    psf_all = combine_psfs(np.stack(source_psfs), [1. for i in range(len(source_psfs))])
    return fix_prop_pixellate(psf_all, common_sampling, detector_pitch)

