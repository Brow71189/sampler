import numpy as np
import ctypes


class Sampler(object):
    c_sampler = ''

    def __init__(self, **kwargs):
        self.base_vec_1 = kwargs.get('base_vec_1', None)
        self.base_vec_2 = kwargs.get('base_vec_1', None)
        self.angle = kwargs.get('angle', None)
        self.offset = kwargs.get('offset', None)
        self.image = kwargs.get('image', None)
        self.mask = kwargs.get('mask', None)
        self.average_unitcell = None

    def calculate_base_from_fft(self, fft_vec_1, fft_vec_2):
        pass

    def load_c_sampler(self):
        pass

    def sample_image(self):
        pass