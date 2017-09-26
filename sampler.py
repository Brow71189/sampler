import numpy as np
import ctypes
import platform

class Sampler(object):
    c_sampler_path = '/home/mittelberger2/git/sampler/sampleUnitCell.so'

    def __init__(self, **kwargs):
        self.base_vec_1 = kwargs.get('base_vec_1', None)
        self.base_vec_2 = kwargs.get('base_vec_1', None)
        self.sample_rate = kwargs.get('sample_rate', None)
        self.offset = kwargs.get('offset', np.array((0,0)))
        self._image = None
        self.image = kwargs.get('image', None)
        self._image_shape = None
        self.mask = kwargs.get('mask', None)
        self._average_unitcell = None
        self.average_unitcell_shape = kwargs.get('average_unitcell_shape', None)
        self._pretty_unitcell = None
        self.periodic_repeats = kwargs.get('periodic_repeats', None)
        self._c_sampler = None
        self.number_unitcells = None
        self.median_unitcell = None
        self._unitcell_data = None
        self.unitcell_stack = None

    @property
    def image(self):
        if self._image is not None:
            return np.reshape(self._image, self._image_shape)

    @image.setter
    def image(self, image):
        if image is not None:
            if image.dtype != np.int32:
                raise ValueError('Image has to be of type np.int32')
            self._image_shape = image.shape
            self._image = np.ascontiguousarray(image.ravel())

    @property
    def average_unitcell(self):
        return np.reshape(self._average_unitcell, self.average_unitcell_shape)

    @average_unitcell.setter
    def average_unitcell(self, average_unitcell):
        self._average_unitcell = np.ascontiguousarray(average_unitcell.ravel())

    @property
    def pretty_unitcell(self):
        #output_shape = np.array(self.average_unitcell_shape, dtype=np.float) * self.periodic_repeats
        #output_shape = (int(output_shape[0]), int(output_shape[1]))
        output_shape = np.array((256,256),int)        
        return np.reshape(self._pretty_unitcell, output_shape)

    @pretty_unitcell.setter
    def pretty_unitcell(self, pretty_unitcell):
        self._pretty_unitcell = np.ascontiguousarray(pretty_unitcell.ravel())

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
        mask = self.mask
        if mask is None:
            mask = np.ascontiguousarray(np.ones(self._image.shape, dtype=np.uint8))
        assert self._image.shape == mask.shape, 'Image and mask have to have the same shape'

        c_int32_p = ctypes.POINTER(ctypes.c_int32)
        image_p = self._image.ctypes.data_as(c_int32_p)

        c_uint8_p = ctypes.POINTER(ctypes.c_uint8)
        mask_p = mask.ctypes.data_as(c_uint8_p)

        c_float_p2 = ctypes.POINTER(ctypes.c_float)
        self.average_unitcell = np.zeros(tuple(self.average_unitcell_shape), dtype=np.float32)
        average_uc_p = self._average_unitcell.ctypes.data_as(c_float_p2)

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

        import time
        starttime = time.time()
        res =  self._c_sampler.sampleUnitCell(image_p, mask_p, shape_x, shape_y, average_uc_p, uc_shape_x, uc_shape_y,
                                              base_vec_1_a, base_vec_1_b, base_vec_2_a, base_vec_2_b, offset_x, offset_y,
                                              sample_rate)
        print('Runtime: {:f} s'.format(time.time() - starttime))
        self.number_unitcells = res
        return res

    def make_pretty_output(self, order, sym):
        output_shape = np.array((256,256),int)
        #enforce boolean
        sym = bool(sym)
        c_float_p = ctypes.POINTER(ctypes.c_float)
        self.pretty_unitcell = np.zeros(output_shape, dtype=np.float32)
        pretty_uc_p = self._pretty_unitcell.ctypes.data_as(c_float_p)

        uc_shape_x = ctypes.c_int32()
        uc_shape_x.value = output_shape[1]
        uc_shape_y = ctypes.c_int32()
        uc_shape_y.value = output_shape[0]

        sample_rate = ctypes.c_int32()
        sample_rate.value = self.sample_rate

        zoom = ctypes.c_double()
        zoom.value = 2.0
        
        c_order = ctypes.c_int32()
        c_order.value = order
        
        mirrors = ctypes.c_bool()
        mirrors.value = sym

        res = self._c_sampler.viewUnitCell(pretty_uc_p, uc_shape_x, uc_shape_y, sample_rate, zoom, c_order, mirrors)
        if res == -1:
            raise RuntimeError('You have to run a sampling before displaying the unitcell')

#        angle = np.cos(np.dot(self.base_vec_1,self.base_vec_2)/(np.linalg.norm(self.base_vec_1)*np.linalg.norm(self.base_vec_2)))
#        aff_mat = np.array(
#                           [[     1, 0 ],
#                            [-angle, 1]]
#                           )
#        output_shape = np.array(self.average_unitcell_shape, dtype=np.int) * self.periodic_repeats
#        img = scipy.ndimage.affine_transform(self.average_unitcell, aff_mat,
#                                             output_shape=(int(output_shape[0]), int(output_shape[1])), mode='wrap')
#        return img

    def view_moment(self, order=1):

        if not (hasattr(self, '_moment_{:.0f}'.format(order)) and hasattr(self, 'moment_{:.0f}'.format(order))):
            setattr(self, '_moment_{:.0f}'.format(order), None)
            def fget(self):
                return getattr(self, ('_moment_{:.0f}'.format(order))).reshape(self.average_unitcell_shape)

            def fset(self, moment):
                setattr(self, '_moment_{:.0f}'.format(order), np.ascontiguousarray(moment.ravel()))

            setattr(Sampler, 'moment_{:.0f}'.format(order), property(fget=fget, fset=fset))

        c_float_p = ctypes.POINTER(ctypes.c_float)
        setattr(self, 'moment_{:.0f}'.format(order), np.zeros(self.average_unitcell_shape, dtype=np.float32))
        moment_p = getattr(self, '_moment_{:.0f}'.format(order)).ctypes.data_as(c_float_p)
        c_order = ctypes.c_int32()
        c_order.value = order
        res = self._c_sampler.viewMoment(c_order, moment_p)
        if res == -1:
            raise RuntimeError('You have to run a sampling before displaying a moment')

    def get_unitcell_stack(self):
        assert self.number_unitcells > 0, 'You need at least one successful sampling run before finding mean'
        c_int32_p = ctypes.POINTER(ctypes.c_int32)
        self._unitcell_data = np.empty(self.number_unitcells*np.prod(self.average_unitcell_shape), dtype=np.int32)
        unitcell_data_p = self._unitcell_data.ctypes.data_as(c_int32_p)
        res = self._c_sampler.getUnitCells(unitcell_data_p)
        assert res == self.number_unitcells*np.prod(self.average_unitcell_shape), res
        self.unitcell_stack = self._unitcell_data.reshape(tuple(self.average_unitcell_shape) + (self.number_unitcells,))
        self.unitcell_stack = np.swapaxes(self.unitcell_stack, 0, 2)
        self.unitcell_stack = np.swapaxes(self.unitcell_stack, 1, 2)

    def get_median(self):
        self.get_unitcell_stack()
        self.median_unitcell = np.median(self.unitcell_stack, axis=0)
        
    def get_internal_offset(self):
        self._c_sampler.get_offset_X.restype = ctypes.c_double
        coffX = self._c_sampler.get_offset_X()
        self._c_sampler.get_offset_Y.restype = ctypes.c_double
        coffY = self._c_sampler.get_offset_Y()
        return (coffX, coffY)

def calculate_counts(image, threshold=1e-9):
    """
    Returns the divisor to translate float values in "image" to actual counts.
    """
    #set all values <0 to 0
    image[image<0] = 0.0
    #flatten and sort image by pixel values
    sort_im = np.sort(np.ravel(image))
    #find "steps" in intensity

    differences = sort_im[1:] - sort_im[0:-1]
    steps = differences[differences>threshold]
    #int_steps = []

    min_step = np.amin(steps)

    int_steps = steps[steps<1.5*min_step]
    counts_divisor = np.mean(int_steps)
    image[image<0]=0.0
    image = np.asarray(np.rint(image/counts_divisor), dtype=np.int32)

    return image
