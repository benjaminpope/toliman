import math
import numpy as np

def gen_phasemap(phase_func, ngrid, sampling, use_cached=True, save_cached=True):
    """Generate the phase map for a phase pupil
    
    Parameters
    ----------
    phase_func : function
        Function to calculate the phase shift at a given polar co-ordinate 
                
    ngrid : integer
        Size of (square) grid for wavefront
    
    sampling : float
        Sampling distance for grid, in metres
    
        
    Returns
    -------
    phase_map: nparray
        Phase map of optical path difference for each position of wavefront
    """

    cached_name = '{}_{}_{}.npy'.format(phase_func.__name__, ngrid, sampling)
    cached_exists = False
    if use_cached is True:
        try:
#            print("Using cached file {}".format(cached_name))
            phase_map = np.load(cached_name)
            cached_exists = True
        except IOError:
            print("Couldn't load file {}".format(cached_name))
            cached_exists = False
    
    if cached_exists is False:
        phase_map = np.zeros([ngrid, ngrid], dtype = np.float64)
        c = ngrid/2.
        for i in range(ngrid):
            for j in range(ngrid):
                x = i - c
                y = j - c
                phi = math.atan2(y, x)
                r = sampling*math.hypot(x,y)
                phase_map[i][j] = phase_func(r, phi)
            
    if cached_exists is False and save_cached is True:
#        print("Caching file {}".format(cached_name))
        np.save(cached_name, phase_map)
        
    return phase_map