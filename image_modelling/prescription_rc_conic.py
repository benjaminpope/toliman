import proper
import math
from prop_conic import prop_conic
from prop_tilt import prop_tilt

def prescription_rc_conic(wavelength, gridsize, 
                          PASSVALUE = {'diam': 0.3, 'm1_fl': 0.5717255, 'm1_m2_sep': 0.549337630333726, 
                                       'm2_fl': -0.023378959,  'bfl': 0.528110658881, 
                                       'beam_ratio': 0.2, 'm2_rad': 0.059,
                                       'm2_strut_width': 0.01,
                                       'm2_supports': 5,
                                       'm1_conic': -1.0,
                                       'm2_conic': -1.0,
                                       'tilt_x': 0.0,
                                       'tilt_y': 0.0
                                      }):
    diam           = PASSVALUE['diam']           # telescope diameter in meters
    m1_fl          = PASSVALUE['m1_fl']          # primary focal length (m)
    m1_m2_sep      = PASSVALUE['m1_m2_sep']      # primary to secondary separation (m)
    m2_fl          = PASSVALUE['m2_fl']          # secondary focal length (m)
    bfl            = PASSVALUE['bfl']            # nominal distance from secondary to focus (m)
    beam_ratio     = PASSVALUE['beam_ratio']     # initial beam width/grid width
    m2_rad         = PASSVALUE['m2_rad']         # Secondary half-diameter (m)
    m2_strut_width = PASSVALUE['m2_strut_width'] # Width of struts supporting M2 (m)
    m2_supports    = PASSVALUE['m2_supports']    # Number of support structs (assumed equally spaced)
    m1_conic       = PASSVALUE['m1_conic']       # Conic constant for M1
    m2_conic       = PASSVALUE['m2_conic']       # Conic constant for M2
    
    tilt_x         = PASSVALUE['tilt_x']         # Tilt angle along x (arc seconds)
    tilt_y         = PASSVALUE['tilt_y']         # Tilt angle along y (arc seconds)

    # Define the wavefront
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)

    # Point off-axis
    prop_tilt(wfo, tilt_x, tilt_y)
    
    
    # Input aperture
    proper.prop_circular_aperture(wfo, diam/2)
    # NOTE: could prop_propagate() here if some baffling included
    # Secondary and structs obscuration
    proper.prop_circular_obscuration(wfo, m2_rad)                      # secondary mirror obscuration
    # Spider struts/vanes, arranged evenly radiating out from secondary
    strut_length = diam/2 - m2_rad
    strut_step = 360/m2_supports
    strut_centre = m2_rad + strut_length/2
    for i in range(0, m2_supports):
        angle = i*strut_step
        radians = math.radians(angle) 
        xoff = math.cos(radians)*strut_centre
        yoff = math.sin(radians)*strut_centre
        proper.prop_rectangular_obscuration(wfo, m2_strut_width, strut_length,
                                            xoff, yoff,
                                            ROTATION = angle + 90)
        
    # Define entrance
    proper.prop_define_entrance(wfo)

    # Primary mirror (treat as conic lens)
    prop_conic(wfo, m1_fl, m1_conic, "primary")

    # Propagate the wavefront
    proper.prop_propagate(wfo, m1_m2_sep, "secondary")

    # Secondary mirror (another conic lens)
    prop_conic(wfo, m2_fl, m2_conic, "secondary")
    
    # NOTE: hole through primary?

    # Focus
    proper.prop_propagate(wfo, bfl, "focus", TO_PLANE=True)

    # End
    (wfo, sampling) = proper.prop_end(wfo)

    return (wfo, sampling)
