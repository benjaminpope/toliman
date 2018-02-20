import proper
import math

def prescription_rc_quad(wavelength, gridsize, 
                          PASSVALUE = {'diam': 0.3, 'm1_fl': 0.5717255, 'm1_m2_sep': 0.549337630333726, 
                                       'm2_fl': -0.023378959,  'bfl': 0.528110658881, 
                                       'beam_ratio': 0.2, 'm2_rad': 0.059,
                                       'm2_strut_width': 0.01,
                                       'm2_supports': 5
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
    
    # Define the wavefront
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)
    
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

    # Primary mirror (treat as quadratic lens)
    proper.prop_lens(wfo, m1_fl, "primary")

    # Propagate the wavefront
    proper.prop_propagate(wfo, m1_m2_sep, "secondary")

    # Secondary mirror (another quadratic lens)
    proper.prop_lens(wfo, m2_fl, "secondary")
    
    # NOTE: hole through primary?

    # Focus
    proper.prop_propagate(wfo, bfl, "focus", TO_PLANE=True)

    # End
    (wfo, sampling) = proper.prop_end(wfo)

    return (wfo, sampling)
