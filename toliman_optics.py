import poppy
from poppy import utils
import astropy.units as u
from poppy import AnalyticOpticalElement
from poppy.poppy_core import OpticalElement, Wavefront, PlaneType, _PUPIL, _IMAGE, _RADIANStoARCSEC
import numpy as np

class EllipticalAperture(AnalyticOpticalElement):
    """ Defines an ideal circular pupil aperture
    Parameters
    ----------
    name : string
        Descriptive name
    radius : float
        Radius of the pupil, in meters. Default is 1.0
    pad_factor : float, optional
        Amount to oversize the wavefront array relative to this pupil.
        This is in practice not very useful, but it provides a straightforward way
        of verifying during code testing that the amount of padding (or size of the circle)
        does not make any numerical difference in the final result.
    """

    @utils.quantity_input(semimajor=u.meter)
    def __init__(self, name=None, semimajor=1.0*u.meter, ecc=0.5, pad_factor=1.0, planetype=PlaneType.unspecified, 
                 **kwargs):

        if name is None:
            name = "Ellipse, a={}, e={}".format(semimajor,ecc)
        super(EllipticalAperture, self).__init__(name=name, planetype=planetype, **kwargs)
        self.semimajor=semimajor
        self.ecc=ecc
        # for creating input wavefronts - let's pad a bit:
        self.pupil_diam = pad_factor * 2 * self.semimajor

    def get_transmission(self, wave):
        """ Compute the transmission inside/outside of the aperture.
        """
        if not isinstance(wave, Wavefront):  # pragma: no cover
            raise ValueError("EllipticalAperture get_transmission must be called with a Wavefront "
                             "to define the spacing")
        assert (wave.planetype != _IMAGE)

        y, x = self.get_coordinates(wave)
        semimajor = self.semimajor.to(u.meter).value
        semiminor = semimajor*(1.-self.ecc**2.)
        r = np.sqrt((x/semimajor) ** 2 + (y/semiminor) ** 2)
        del x
        del y

        w_outside = np.where(r > 1.)
        del r
        self.transmission = np.ones(wave.shape)
        self.transmission[w_outside] = 0
        return self.transmission
    
    def get_opd(self, wave):
        transmission = self.get_transmission(wave)
        phase = 325e-9
        return phase * transmission

class CompositeAnalyticOptic(poppy.CompoundAnalyticOptic):

    # CompoundAnalyticOptic starts off with ones and multiplies by each element. This
    # version starts off with zero and adds
    def get_transmission(self,wave):
        trans = np.zeros(wave.shape, dtype=np.float)
        for optic in self.opticslist:
            trans += optic.get_transmission(wave)
        return trans

    # CompoundAnalyticOptic starts off with ones and multiplies each optic's phasor
    # This version starts with zeros and adds each optic's scaled phasor
    def get_phasor(self, wave):
        phasor = np.zeros(wave.shape, dtype=np.complex)
        for optic in self.opticslist:
            nextphasor = optic.get_phasor(wave)/float(len(self.opticslist))
            phasor += nextphasor
        return phasor

def rosette(symm, radius, ecc, semimajor):
    """Create a custom TOLIMAN rosette
    
    Defines both an elliptical aperture, and a composite aperture which sums individual subapertures.
    
    Keyword arguments:
    symm -- Degree of rotational symmetry
    radius -- radial position of elliptical apertures (in metres)
    ecc -- eccentricity of each aperture ellipse
    semimajor -- length of each aperture ellipse's semimajor (in metres)
    """
    aps = []

    for j in np.arange(symm):
        pa = 2.*np.pi/symm*(j+1.)
        xc, yc = radius*np.cos(pa), radius*np.sin(pa)
        aps.append(EllipticalAperture(semimajor=semimajor,ecc=ecc,shift_x=xc,shift_y=yc,rotation=180./np.pi*pa,
                                      pad_factor=6))

    for j in np.arange(symm):
        pa = 2.*np.pi/symm*(j+1.) + np.pi/symm
        xc, yc = radius*np.cos(pa), radius*np.sin(pa)
        aps.append(EllipticalAperture(semimajor=semimajor,ecc=ecc,shift_x=xc,shift_y=yc,rotation=180./np.pi*pa +90.,
                                      pad_factor=6))
    return CompositeAnalyticOptic(aps)

def phase_rosette(edge, symm, radius, ecc, semimajor, phase):
    """Create a custom TOLIMAN rosette
    
    Defines both an elliptical aperture, and a composite aperture which sums individual subapertures.
    
    Keyword arguments:
    edge -- radius of edge of element
    symm -- Degree of rotational symmetry
    radius -- radial position of elliptical apertures (in metres)
    ecc -- eccentricity of each aperture ellipse
    semimajor -- length of each aperture ellipse's semimajor (in metres)
    phase -- phase shift of rosette
    """
    aps = []
    aps.append(poppy.CircularAperture(edge))

    for j in np.arange(symm):
        pa = 2.*np.pi/symm*(j+1.)
        xc, yc = radius*np.cos(pa), radius*np.sin(pa)
        aps.append(EllipticalAperture(semimajor=semimajor,ecc=ecc,shift_x=xc,shift_y=yc,rotation=180./np.pi*pa,
                                      pad_factor=6))

    for j in np.arange(symm):
        pa = 2.*np.pi/symm*(j+1.) + np.pi/symm
        xc, yc = radius*np.cos(pa), radius*np.sin(pa)
        aps.append(EllipticalAperture(semimajor=semimajor,ecc=ecc,shift_x=xc,shift_y=yc,rotation=180./np.pi*pa +90.,
                                      pad_factor=6))
#    return CompositeAnalyticOptic(aps)
    return poppy.CompoundAnalyticOptic(aps)

# poppy.fresnel.ConicLens is broken (as of 0.6.1) so this is a local version
# that hopefully fixes it.
class ConicLens(poppy.CircularAperture):
    @u.quantity_input(f_lens=u.m, radius=u.m)
    def __init__(self,
                 f_lens=1.0 * u.m,
                 K=1.0,
                 radius=1.0 * u.m,
                 planetype=PlaneType.unspecified,
                 name="Conic lens",
                **kwargs):
        """Conic Lens/Mirror
        Parabolic, elliptical, hyperbolic, or spherical powered optic.
        Parameters
        ----------------
        f_lens : astropy.quantities.Quantity of dimension length
            Focal length of the optic
        K : float
            Conic constant
        radius: astropy.quantities.Quantity of dimension length
            Radius of the clear aperture of the optic as seen on axis.
        name : string
            Descriptive name
        planetype : poppy.PlaneType, optional
            Optional optical plane type specifier
        """
        CircularAperture.__init__(self, name=name, radius=radius.to(u.m).value, planetype=planetype, **kwargs)
        self.f_lens = f_lens
        self.K = K
