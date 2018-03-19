import proper
import math
import numpy as np
from prop_tilt import prop_tilt

def prescription_quad(wavelength, gridsize, PASSVALUE = {}):
    # Assign parameters from PASSVALUE struct or use defaults
    diam           = PASSVALUE.get('diam',0.3)                    # telescope diameter in meters
    m1_fl          = PASSVALUE.get('m1_fl',0.5717255)             # primary focal length (m)
    beam_ratio     = PASSVALUE.get('beam_ratio',0.2)              # initial beam width/grid width
    tilt_x         = PASSVALUE.get('tilt_x',0.)                   # Tilt angle along x (arc seconds)
    tilt_y         = PASSVALUE.get('tilt_y',0.)                   # Tilt angle along y (arc seconds)
    noabs          = PASSVALUE.get('noabs',False)                 # Output complex amplitude?
    """Prescription for a single quad lens system        
    """
    
    # Define the wavefront
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)
    
    # Point off-axis
    prop_tilt(wfo, tilt_x, tilt_y)
    
    # Input aperture
    proper.prop_circular_aperture(wfo, diam/2.)
        
    # Define entrance
    proper.prop_define_entrance(wfo)

    # Primary mirror (treat as quadratic lens)
    proper.prop_lens(wfo, m1_fl, "primary")
    if 'phase_func' in PASSVALUE:
        phase_func = PASSVALUE['phase_func']
        ngrid = proper.prop_get_gridsize(wfo)
        sampling = proper.prop_get_sampling(wfo)
        phase_map = np.zeros([ngrid, ngrid], dtype = np.float64)
        c = ngrid/2.
        opd = wavelength/4.
        for i in range(ngrid):
            for j in range(ngrid):
                x = i - c
                y = j - c
                phi = math.atan2(y, x)
                r = sampling*math.hypot(x,y)
                phase_map[i][j] = phase_func(r, phi, opd)
        
        proper.prop_add_phase(wfo, phase_map)

    # Focus
    proper.prop_propagate(wfo, m1_fl, "focus", TO_PLANE=True)

    # End
    (wfo, sampling) = proper.prop_end(wfo)

    return (wfo, sampling)
