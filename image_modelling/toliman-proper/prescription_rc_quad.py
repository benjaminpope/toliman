import proper
import math
from prop_conic import prop_conic
from prop_tilt import prop_tilt

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
    noabs          = PASSVALUE.get('noabs',False)                 # Tilt angle along y (arc seconds)
    
    
    # Define the wavefront
    wfo = proper.prop_begin(diam, wavelength, gridsize, beam_ratio)

    if proper.prop_is_statesaved(wfo) == False:
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

        # Normalize wavefront
        proper.prop_define_entrance(wfo)

        # Primary mirror
        if 'm1_conic' in PASSVALUE:
            prop_conic(wfo, m1_fl, PASSVALUE['m1_conic'], "conic primary")
        else:
            proper.prop_lens(wfo, m1_fl, "primary")
        proper.prop_circular_obscuration(wfo, m1_hole_rad)


        # Propagate the wavefront
        proper.prop_propagate(wfo, m1_m2_sep, "secondary")

        # Secondary mirror
        if 'm1_conic' in PASSVALUE:
            prop_conic(wfo, m2_fl, PASSVALUE['m2_conic'], "conic secondary")
        else:
            proper.prop_lens(wfo, m2_fl, "secondary")
        proper.prop_circular_aperture(wfo, m2_rad)    

        # Hole through primary
        proper.prop_propagate(wfo, m1_m2_sep, "M1 hole")
        proper.prop_circular_aperture(wfo, m1_hole_rad)    

    proper.prop_state(wfo)

    # Focus - bfl can be varied between runs
    proper.prop_propagate(wfo, bfl-m1_m2_sep, "focus", TO_PLANE=True)

    # End
    return proper.prop_end(wfo, NOABS = noabs)
