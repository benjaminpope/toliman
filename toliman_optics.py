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
            raise ValueError("CircularAperture get_transmission must be called with a Wavefront "
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

class CompositeAnalyticOptic(AnalyticOpticalElement):
    """ Define a compound analytic optical element made up of the combination
    of two or more individual optical elements.
    This is just a convenience routine for semantic organization of optics.
    It can be useful to keep the list of optical planes cleaner, but
    you can certainly just add a whole bunch of planes all in a row without
    using this class to group them.
    All optics should be of the same plane type (pupil or image); propagation between
    different optics contained inside one compound is not supported.
    Parameters
    ----------
    opticslist : list
        A list of AnalyticOpticalElements to be merged together.
    """

    def _validate_only_analytic_optics(self, optics_list):
        for optic in optics_list:
            if isinstance(optic, AnalyticOpticalElement):
                continue  # analytic elements are allowed
            elif isinstance(optic, InverseTransmission):
                if isinstance(optic.uninverted_optic, AnalyticOpticalElement):
                    continue  # inverted elements are allowed, as long as they're analytic elements
                else:
                    return False  # inverted non-analytic elements aren't allowed, skip the rest
            else:
                return False  # no other types allowed, skip the rest of the list
        return True

    def __init__(self, opticslist=None, name="unnamed", verbose=True, **kwargs):
        if opticslist is None:
            raise ValueError("Missing required opticslist argument to CompoundAnalyticOptic")
        AnalyticOpticalElement.__init__(self, name=name, verbose=verbose, **kwargs)

        self.opticslist = []
        self._default_display_size = 3*u.arcsec
        self.planetype = None

        for optic in opticslist:
            if not self._validate_only_analytic_optics(opticslist):
                raise ValueError("Supplied optics list to CompoundAnalyticOptic can "
                                 "only contain AnalyticOptics")
            else:
                # if we are adding the first optic in the list, check what type of optical plane
                # it has
                # for subsequent optics, validate they have the same type
                if len(self.opticslist) == 0:
                    self.planetype = optic.planetype
                elif (self.planetype != optic.planetype and self.planetype != PlaneType.unspecified and
                        optic.planetype != PlaneType.unspecified):
                    raise ValueError("Cannot mix image plane and pupil plane optics in "
                                     "the same CompoundAnalyticOptic")

                self.opticslist.append(optic)
                if hasattr(optic, '_default_display_size'):
                    self._default_display_size = max(self._default_display_size,
                                                     optic._default_display_size)
                if hasattr(optic,'pupil_diam'):
                    if not hasattr(self,'pupil_diam'):
                        self.pupil_diam = optic.pupil_diam
                    else:
                        self.pupil_diam = max(self.pupil_diam, optic.pupil_diam)

        if self.planetype == _PUPIL:
            if all([hasattr(o, 'pupil_diam') for o in self.opticslist]):
                self.pupil_diam = np.asarray([o.pupil_diam.to(u.meter).value for o in self.opticslist]).max() * u.meter

    def get_transmission(self,wave):
        trans = np.zeros(wave.shape, dtype=np.float)
        for optic in self.opticslist:
            trans += optic.get_transmission(wave)
        return trans

    def get_opd(self,wave):
        opd = np.zeros(wave.shape, dtype=np.float)
        for optic in self.opticslist:
            opd += optic.get_opd(wave)
        return opd

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