import math
import numpy as np

def gen_opdmap(opd_func, ngrid, sampling, use_cached=True, save_cached=True):
    """Generate the OPD map for a phase pupil
    
    Parameters
    ----------
    opd_func : function
        Function to calculate the OPD introduced at a given polar co-ordinate 
                
    ngrid : integer
        Size of (square) grid for wavefront
    
    sampling : float
        Sampling distance for grid, in metres
    
        
    Returns
    -------
    opd_map: nparray
        OPD map of optical path difference for each position of wavefront
    """

    cached_name = '{}_{}_{}.npy'.format(opd_func.__name__, ngrid, sampling)
    cached_exists = False
    if use_cached is True:
        try:
#            print("Using cached file {}".format(cached_name))
            opd_func = np.load(cached_name)
            cached_exists = True
        except IOError:
            print("Couldn't load file {}".format(cached_name))
            cached_exists = False
    
    if cached_exists is False:
        opd_map = np.zeros([ngrid, ngrid], dtype = np.float64)
        c = ngrid/2.
        for i in range(ngrid):
            for j in range(ngrid):
                x = i - c
                y = j - c
                phi = math.atan2(y, x)
                r = sampling*math.hypot(x,y)
                opd_map[i][j] = opd_func(r, phi)
            
    if cached_exists is False and save_cached is True:
#        print("Caching file {}".format(cached_name))
        np.save(cached_name, opd_map)
        
    return opd_map