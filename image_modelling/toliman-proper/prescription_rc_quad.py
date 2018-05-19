import proper
import math
import numpy as np
from prop_conic import prop_conic
from prop_tilt import prop_tilt
from gen_opdmap import gen_opdmap
from build_prop_circular_aperture import build_prop_circular_aperture
from build_prop_circular_obscuration import build_prop_circular_obscuration
from build_prop_rectangular_obscuration import build_prop_rectangular_obscuration
from build_phase_map import build_phase_map

def prescription_rc_quad(wavelength, gridsize, PASSVALUE = {}):
    # Assign parameters from PASSVALUE struct or use defaults
    diam           = PASSVALUE.get('diam',0.3)                    # telescope diameter in meters
    m1_fl          = PASSVALUE.get('m1_fl',0.5717255)             # primary focal length (m)
    m1_hole_rad    = PASSVALUE.get('m1_hole_rad',0.035)           # Radius of hole in primary (m)
    m1_m2_sep      = PASSVALUE.get('m1_m2_sep',0.549337630333726) # primary to secondary separation (m)
    m2_fl          = PASSVALUE.get('m2_fl',-0.023378959)          # secondary focal length (m)
    bfl            = PASSVALUE.get('bfl',0.528110658881)          # nominal distance from secondary to focus (m)
    beam_ratio     = PASSVALUE.get('beam_ratio',0.2)              # initial beam width/grid width
    m2_rad         = PASSVALUE.get('m2_rad',0.059)                # Secondary half-diameter (m)
    m2_strut_width = PASSVALUE.get('m2_strut_width',0.01)         # Width of struts supporting M2 (m)
    m2_supports    = PASSVALUE.get('m2_supports',5)               # Number of support structs (assumed equally spaced)
    tilt_x         = PASSVALUE.get('tilt_x',0.)                   # Tilt angle along x (arc seconds)
    tilt_y         = PASSVALUE.get('tilt_y',0.)                   # Tilt angle along y (arc seconds)
    noabs          = PASSVALUE.get('noabs',False)                 # Output complex amplitude?
    use_caching    = PASSVALUE.get('use_caching',False)           # Use cached files if available?
    # Can also specify a phase_func function with signature phase_func(r, phi)
    
    # Define the wavefront
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)

# Disable state saving as by default this saves state even when not used.
#    if proper.prop_is_statesaved(wfo) == False:
    # Point off-axis
    prop_tilt(wfo, tilt_x, tilt_y)

    # Input aperture
    wfo.wfarr *= build_prop_circular_aperture(wfo, diam/2)
    # NOTE: could prop_propagate() here if some baffling included

    # Secondary and structs obscuration
    wfo.wfarr *= build_prop_circular_obscuration(wfo, m2_rad) # secondary mirror obscuration
    # Spider struts/vanes, arranged evenly radiating out from secondary
    strut_length = diam/2 - m2_rad
    strut_step = 360/m2_supports
    strut_centre = m2_rad + strut_length/2
    for i in range(0, m2_supports):
        angle = i*strut_step
        radians = math.radians(angle) 
        xoff = math.cos(radians)*strut_centre
        yoff = math.sin(radians)*strut_centre
        wfo.wfarr *= build_prop_rectangular_obscuration(wfo, m2_strut_width, 									strut_length,
                                            xoff, yoff,
                                            ROTATION = angle + 90)

    # Normalize wavefront
    proper.prop_define_entrance(wfo)

    proper.prop_propagate(wfo, m1_m2_sep, "primary")
    # Primary mirror
    if 'phase_func' in PASSVALUE:
        phase_func = PASSVALUE['phase_func']
        ngrid = proper.prop_get_gridsize(wfo)
        sampling = proper.prop_get_sampling(wfo)
        opd_map = gen_opdmap(phase_func, ngrid, sampling, use_cached=use_caching, save_cached=use_caching)
        wfo.wfarr *= build_phase_map(wfo, opd_map)
    if 'm1_conic' in PASSVALUE:
        prop_conic(wfo, m1_fl, PASSVALUE['m1_conic'], "conic primary")
    else:
        proper.prop_lens(wfo, m1_fl, "primary")
    wfo.wfarr *= build_prop_circular_obscuration(wfo, m1_hole_rad)

    # Secondary mirror
    proper.prop_propagate(wfo, m1_m2_sep, "secondary")
    if 'phase_func_sec' in PASSVALUE:
        phase_func = PASSVALUE['phase_func_sec']
        ngrid = proper.prop_get_gridsize(wfo)
        sampling = proper.prop_get_sampling(wfo)
        opd_map = gen_opdmap(phase_func, ngrid, sampling, use_cached=use_caching, save_cached=use_caching)
        wfo.wfarr *= build_phase_map(wfo, opd_map)
    if 'm1_conic' in PASSVALUE:
        prop_conic(wfo, m2_fl, PASSVALUE['m2_conic'], "conic secondary")
    else:
        proper.prop_lens(wfo, m2_fl, "secondary")
    wfo.wfarr *= build_prop_circular_aperture(wfo, m2_rad)    

#    proper.prop_state(wfo)

    # Hole through primary
    if m1_m2_sep<bfl:
        proper.prop_propagate(wfo, m1_m2_sep, "M1 hole")
        wfo.wfarr *= build_prop_circular_aperture(wfo, m1_hole_rad)    


    # Focus - bfl can be varied between runs
    if m1_m2_sep<bfl:
        proper.prop_propagate(wfo, bfl-m1_m2_sep, "focus", TO_PLANE=True)
    else:
        proper.prop_propagate(wfo, bfl, "focus", TO_PLANE=True)

    # End
    return proper.prop_end(wfo, NOABS = noabs)
