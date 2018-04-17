import proper
import math
import numpy as np
import matplotlib.pyplot as plt
from prop_tilt import prop_tilt
from gen_phasemap import gen_phasemap

def prescription_quad(wavelength, gridsize, PASSVALUE = {}):
    # Assign parameters from PASSVALUE struct or use defaults
    diam           = PASSVALUE.get('diam',0.3)                    # telescope diameter in meters
    m1_fl          = PASSVALUE.get('m1_fl',0.5717255)             # primary focal length (m)
    beam_ratio     = PASSVALUE.get('beam_ratio',0.2)              # initial beam width/grid width
    tilt_x         = PASSVALUE.get('tilt_x',0.)                   # Tilt angle along x (arc seconds)
    tilt_y         = PASSVALUE.get('tilt_y',0.)                   # Tilt angle along y (arc seconds)
    noabs          = PASSVALUE.get('noabs',False)                 # Output complex amplitude?
    m1_hole_rad    = PASSVALUE.get('m1_hole_rad',None)            # Inner hole diameter
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
    if 'phase_func' in PASSVALUE:
        phase_func = PASSVALUE['phase_func']
        ngrid = proper.prop_get_gridsize(wfo)
        sampling = proper.prop_get_sampling(wfo)
        phase_map = gen_phasemap(phase_func, ngrid, sampling)
        proper.prop_add_phase(wfo, phase_map)
    if m1_hole_rad is not None:
        proper.prop_circular_obscuration(wfo, m1_hole_rad)
    wf = proper.prop_get_wavefront(wfo)
    #plt.imshow(np.angle(wf))
    #plt.show()
    #plt.imshow(np.abs(wf))
    #plt.show()
    #plt.imshow(np.log10(np.abs(np.fft.fftshift(np.fft.fft2(wf)))))
    #plt.show()
    proper.prop_lens(wfo, m1_fl, "primary")
    
    # Focus
    proper.prop_propagate(wfo, m1_fl, "focus") #, TO_PLANE=True)

    # End
    (wfo, sampling) = proper.prop_end(wfo)

    return (wfo, sampling)
