import numpy as np
import ctypes
import platform

class Sampler(object):
    c_sampler_path = '/home/mittelberger2/git/sampler/sampleUnitCell/sampleUnitCell.so'

    def __init__(self, **kwargs):
        self.base_vec_1 = kwargs.get('base_vec_1', None)
        self.base_vec_2 = kwargs.get('base_vec_1', None)
        self.sample_rate = kwargs.get('sample_rate', None)
        self.offset = kwargs.get('offset', None)
        self.image = kwargs.get('image', None)
        self.mask = kwargs.get('mask', None)
        self.average_unitcell = None
        self.average_unitcell_shape = kwargs.get('average_unitcell_shape', None)
        self._c_sampler = None

    def calculate_base_from_fft(self, fft_vec_1, fft_vec_2):
        a1_s = np.array((fft_vec_1[1], fft_vec_1[0], 0))
        a2_s = np.array((fft_vec_2[1], fft_vec_2[0], 0))
        n = np.cross(a1_s, a2_s)/np.linalg.norm(np.cross(a1_s, a2_s))
        a1 = np.cross(a2_s, n)/np.linalg.norm(np.cross(a1_s,a2_s))
        a2 = np.cross(n, a1_s)/np.linalg.norm(np.cross(a1_s,a2_s))

        return (np.array((a1[1], a1[0])), np.array((a2[1], a2[0])))

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
            void sampleImage(int_16* average_uc, int_16* image, int_8* mask, int_16 shape_y, int_16 shape_x,
                             int_8 sample_rate, int_16 shape_uc_y, int_16 shape_uc_x, double base_vec_1_a,
                             double base_vec_1_b, double base_vec_2_a, double base_vec_2_b)
        """
        assert self._c_sampler is not None, 'Sampler not loaded, please call "load_c_sampler" first.'
        assert self.image is not None, 'No image to process'
        assert self.average_unitcell_shape is not None, 'Need to know the shape of the resulting unitcell'
        if self.mask is None:
            self.mask = np.ones(self.image.shape)
        assert self.image.shape == self.mask.shape, 'Image and mask have to have the same shape'


        c_int32_p = ctypes.POINTER(ctypes.c_int32)
        image_p = self.image.ctypes.data_as(c_int32_p)

        c_uint8_p = ctypes.POINTER(ctypes.c_uint8)
        mask_p = self.mask.ctypes.data_as(c_uint8_p)

        c_int32_p2 = ctypes.POINTER(ctypes.c_int32)
        self.average_unitcell = np.zeros(tuple(self.average_unitcell_shape), dtype=np.int32)
        average_uc_p = self.average_unitcell.ctypes.data_as(c_int32_p2)

        shape_y = ctypes.c_int32()
        shape_y.value = self.image.shape[0]

        shape_x = ctypes.c_int32()
        shape_x.value = self.image.shape[1]

        uc_shape_y = ctypes.c_int32()
        uc_shape_y.value = self.average_unitcell_shape[0]

        uc_shape_x = ctypes.c_int32()
        uc_shape_x.value = self.average_unitcell_shape[1]

        sample_rate = ctypes.c_int32()
        sample_rate.value = self.sample_rate

        offset_y = ctypes.c_double()
        offset_y.value = self.offset[0]

        offset_x = ctypes.c_double()
        offset_x.value = self.offset[1]

        base_vec_1_a = ctypes.c_double()
        base_vec_1_a.value = self.base_vec_1[1]

        base_vec_1_b = ctypes.c_double()
        base_vec_1_b.value = self.base_vec_1[0]

        base_vec_2_a = ctypes.c_double()
        base_vec_2_a.value = self.base_vec_2[1]

        base_vec_2_b = ctypes.c_double()
        base_vec_2_b.value = self.base_vec_2[0]

        return self._c_sampler.sampleUnitCell(c_int32_p(), mask_p, shape_x, shape_y, average_uc_p, uc_shape_x, uc_shape_y,
                                              base_vec_1_a, base_vec_1_b, base_vec_2_a, base_vec_2_b, offset_x, offset_y,
                                              sample_rate)