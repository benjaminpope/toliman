import numpy as np
import proper
import glob, os

def gen_cached_name(label, wfo):
    prefix = 'cached'
    # Most of the cacheable steps depend upon wfo for:
    ngrid = proper.prop_get_gridsize(wfo)
    beamradius = proper.prop_get_beamradius(wfo)
    sampling = proper.prop_get_sampling(wfo)
    return '{}_{}_{}_{}_{}.npy'.format(prefix, label, ngrid, sampling, beamradius)

def clear_all_cached():
    prefix = 'cached'
    for f in glob.glob(prefix+"P*.npy"):
        print("Removing cached file {}".format(f))
        os.remove(f)
    
def load_cached_grid(cached_name):
    try:
        grid = np.load(cached_name)
        print("Found cached file {}".format(cached_name))       
    except IOError:
        print("Couldn't load file {}".format(cached_name))
        grid = None 
    return grid

def save_cached_grid(cached_name, grid):
    print("Caching file {}".format(cached_name))
    np.save(cached_name, grid)

# TODO: Use once I've worked out how to do currying in Python
#def load_cacheable_grid(label, wfo, func, use_caching=True):
#   cachename = gen_cached_name(label, wfo)
#   if use_caching:
#       opd = load_cached_grid(cachename)
#   if opd is None:
#       opd = func
#       if use_caching:
#           save_cached_grid(cachename, opd)