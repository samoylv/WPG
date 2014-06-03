from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptL as Lens
from wpg.srwlib import srwl
import wpg.srwlib


class WPGOpticalElement(object):

    """Base class for optical elements"""

    def __init__(self):
        pass


class Empty(WPGOpticalElement):

    """Optical element: Empty
    This is empty propagator used for sampling and zooming wavefront
    """

    def __init__(self):
        super(Empty, self).__init__()

    def __str__(self):
        return ''

    def propagate(self, wfr, propagation_parameters):
        """
        Propagate wavefront through empty propagator,
        used for sampling and resizing wavefront
        """
        beamline = wpg.srwlib.SRWLOptC([], propagation_parameters)
        srwl.PropagElecField(wfr._srwl_wf, beamline)


class Use_PP(object):

    """Short version of propagation parameters"""

    def __init__(self, auto_resize_before=None, auto_resize_after=None,
                 releative_precision=None, semi_analytical_treatment=None,
                 fft_resizing=None,  zoom=None, zoom_h=None, zoom_v=None,
                 sampling=None, sampling_h=None, sampling_v=None, srw_pp=None
                 ):
        """
        :params srw_pp: propagation parameters in srw style
        """

        super(Use_PP, self).__init__()
        # default values fo propagation parameters
        # Wavefront Propagation Parameters:
        #[0]:  Auto-Resize (1) or not (0) Before propagation
        #[1]:  Auto-Resize (1) or not (0) After propagation
        #[2]:  Relative Precision for propagation with Auto-Resizing (1. is nominal)
        #[3]:  Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation
        #[4]:  Do any Resizing on Fourier side, using FFT, (1) or not (0)
        #[5]:  Horizontal Range modification factor at Resizing (1. means no modification)
        #[6]:  Horizontal Resolution modification factor at Resizing
        #[7]:  Vertical Range modification factor at Resizing
        #[8]:  Vertical Resolution modification factor at Resizing
        #[9]:  Type of wavefront Shift before Resizing (not yet implemented)
        #[10]: New Horizontal wavefront Center position after Shift (not yet implemented)
        #[11]: New Vertical wavefront Center position after Shift (not yet implemented)
        #         [ 0] [1] [2]  [3] [4] [5]  [6]  [7]  [8]  [9] [10] [11]
        self.pp = [0,  0, 1.0,  0,  0, 1.0, 1.0, 1.0, 1.0,  0,  0,   0]

        if not srw_pp is None:
            if isinstance(srw_pp, list) and len(srw_pp) == 12:
                self.pp = srw_pp

        if not auto_resize_before is None:
            self.auto_resize_before = auto_resize_before

        if not auto_resize_after is None:
            self.auto_resize_after = auto_resize_after

        if not releative_precision is None:
            self.releative_precision = releative_precision

        if not semi_analytical_treatment is None:
            self.semi_analytical_treatment = semi_analytical_treatment

        if not fft_resizing is None:
            self.fft_resizing = fft_resizing

        if not zoom is None:
            self.zoom = zoom

        if not zoom_h is None:
            self.zoom_h = zoom_h

        if not zoom_v is None:
            self.zoom_v = zoom_v

        if not sampling is None:
            self.sampling = sampling

        if not sampling_h is None:
            self.sampling_h = sampling_h

        if not sampling_v is None:
            self.sampling_v = sampling_v

    @property
    def auto_resize_before(self):
        """Auto-Resize (1) or not (0) Before propagation"""
        return self.pp[0]

    @auto_resize_before.setter
    def auto_resize_before(self, value):
        if value in [0, 1]:
            self.pp[0] = value

    @property
    def auto_resize_after(self):
        """Auto-Resize (1) or not (0) After propagation"""
        return self.pp[1]

    @auto_resize_after.setter
    def auto_resize_after(self, value):
        if value in [0, 1]:
            self.pp[1] = value

    @property
    def releative_precision(self):
        """Relative Precision for propagation with Auto-Resizing (1. is nominal)"""
        return self.pp[2]

    @releative_precision.setter
    def releative_precision(self, value):
        self.pp[2] = float(value)

    @property
    def semi_analytical_treatment(self):
        """Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation"""
        return self.pp[3]

    @semi_analytical_treatment.setter
    def semi_analytical_treatment(self, value):
        if value in [0, 1]:
            self.pp[3] = value

    @property
    def fft_resizing(self):
        """Do any Resizing on Fourier side, using FFT, (1) or not (0)"""
        return self.pp[4]

    @fft_resizing.setter
    def fft_resizing(self, value):
        if value in [0, 1]:
            self.pp[4] = value

    @property
    def zoom_h(self):
        """Horizontal Range modification factor at Resizing (1. means no modification)"""
        return self.pp[5]

    @zoom_h.setter
    def zoom_h(self, value):
        self.pp[5] = float(value)

    @property
    def sampling_h(self):
        """Horizontal Resolution modification factor at Resizing"""
        return self.pp[6]

    @sampling_h.setter
    def sampling_h(self, value):
        self.pp[6] = float(value)

    @property
    def zoom_v(self):
        """Vertical Range modification factor at Resizing"""
        return self.pp[7]

    @zoom_v.setter
    def zoom_v(self, value):
        self.pp[7] = float(value)

    @property
    def sampling_v(self):
        """Vertical Resolution modification factor at Resizing"""
        return self.pp[8]

    @sampling_v.setter
    def sampling_v(self, value):
        self.pp[8] = float(value)

    @property
    def zoom(self):
        """Range modification factor at Resizing (1. means no modification)"""
        if self.zoom_h == self.zoom_v:
            return self.zoom_h
        else:
            return "zoom_h != zoom_v"

    @zoom.setter
    def zoom(self, value):
        z = float(value)
        self.zoom_h = z
        self.zoom_v = z

    @property
    def sampling(self):
        """Resolution modification factor at Resizing"""
        if self.sampling_h == self.sampling_v:
            return self.sampling_h
        else:
            return "sampling_h != sampling_v"

    @sampling.setter
    def sampling(self, value):
        s = float(value)
        self.sampling_h = s
        self.sampling_v = s

    def get_srw_pp(self):
        return self.pp

    def __str__(self):
        return '\n'.join([
            "zoom_h = {}".format(self.zoom_h),
            "zoom_v = {}".format(self.zoom_v),
            "sampling_h = {}".format(self.sampling_h),
            "sampling_v = {}".format(self.sampling_v),
            "semi_analytical_treatment = {}".format(
                self.semi_analytical_treatment),
            "auto_resize_before = {}".format(self.auto_resize_before),
            "auto_resize_after = {}".format(self.auto_resize_after),
            "releative_precision = {}".format(self.releative_precision),
            "fft_resizing = {}".format(self.fft_resizing),
            '\n'
        ])