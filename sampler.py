import numpy as np
import ctypes
import platform

class Sampler(object):
    c_sampler_path = ''

    def __init__(self, **kwargs):
        self.base_vec_1 = kwargs.get('base_vec_1', None)
        self.base_vec_2 = kwargs.get('base_vec_1', None)
        self.angle = kwargs.get('angle', None)
        self.offset = kwargs.get('offset', None)
        self.image = kwargs.get('image', None)
        self.mask = kwargs.get('mask', None)
        self.average_unitcell = None
        self.average_unitcell_shape = kwargs.get('average_unitcell_shape', None)
        self._c_sampler = None

    def calculate_base_from_fft(self, fft_vec_1, fft_vec_2):
        pass

    def load_c_sampler(self):
        # for Windows
        if platform.system() == "Windows":
            self._c_sampler = ctypes.WinDLL(Sampler.c_sampler_path)
        # for Linux
        elif platform.system() == "Linux":
            self._c_sampler = ctypes.cdll.LoadLibrary(Sampler.c_sampler_path)
        else:
            raise RuntimeError('Cannot detect operating system, wil now stop.')

    def sample_image(self):
        """
        Calls the C Code to sample an image.
        Assumed C function:
            void sampleImage(float* average_uc, int* image, byte* mask, float angle, float offset, float base_vec_1_a,
                             float base_vec_1_b, float base_vec_2_a, float base_vec_2_b)
        """
        assert self._c_sampler is not None, 'Sampler not loaded, please call "load_c_sampler" first.'
        assert self.image is not None, 'No image to process'
        assert self.average_unitcell_shape is not None, 'Need to know the shape of the resulting unitcell'
        if self.mask is None:
            self.mask = np.ones(self.image.shape)
        assert self.image.shape == self.mask.shape, 'Image and mask have to have the same shape'
        
        
        c_int16_p = ctypes.POINTER(ctypes.c_uint16)
        image_p = self.image.ctypes.data_as(c_int16_p)
        
        c_int8_p = ctypes.POINTER(ctypes.c_uint8)
        mask_p = self.mask.ctypes.data_as(c_int8_p)
        
        c_float_p = ctypes.POINTER(ctypes.c_float)
        self.average_unitcell = np.zeros(tuple(self.average_unitcell_shape))
        average_uc_p = self.average_unitcell.ctypes.data_as(c_float_p)
        
        angle = ctypes.c_float()
        angle.value = self.angle
        
        offset = ctypes.c_float()
        offset.value = self.offset
        
        base_vec_1_a = ctypes.c_float()
        base_vec_1_a.value = self.base_vec_1[1]
        
        base_vec_1_b = ctypes.c_float()
        base_vec_1_b.value = self.base_vec_1[0]
        
        base_vec_2_a = ctypes.c_float()
        base_vec_2_a.value = self.base_vec_2[1]
        
        base_vec_2_b = ctypes.c_float()
        base_vec_2_b.value = self.base_vec_2[0]
        
        self._c_sampler.sampleImage(average_uc_p, image_p, mask_p, angle, offset, base_vec_1_a, base_vec_1_b,
                                    base_vec_2_a, base_vec_2_b)