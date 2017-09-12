#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:06:50 2017

@author: mittelberger2
"""
import sampler
import tifffile
import numpy as np
import matplotlib.pyplot as plt

sampler.Sampler.c_sampler_path = './sampleUnitCell.so'
sampl = sampler.Sampler()
sampl.load_c_sampler()
sampl.image = sampler.calculate_counts(tifffile.imread('Crop_of_SuperScan.tif'))
sampl.base_vec_1 = np.array((-41.1391, -17.5251))
sampl.base_vec_2 = np.array((-26.8973, 5.7381))


sampl.average_unitcell_shape=(32,32)
sampl.sample_rate=4

print('Number of unitcells: {:d}'.format(sampl.sample_image()))

sampl.periodic_repeats=8
sampl.make_pretty_output()

plt.gray()
plt.close()

tifffile.imsave('Average unitcell.tif', sampl.average_unitcell.astype(np.uint16))
uc = plt.matshow(sampl.average_unitcell)
num = uc.figure.number
uc.figure.canvas.set_window_title('Average unitcell ({:d})'.format(num))

tifffile.imsave('Pretty unitcell.tif', sampl.pretty_unitcell.astype(np.uint16))
puc = plt.matshow(sampl.pretty_unitcell)
num = puc.figure.number
puc.figure.canvas.set_window_title('Pretty unitcell ({:d})'.format(num))

for i in range(1, 5):
    sampl.view_moment(i)
    tifffile.imsave('Moment {:d}'.format(i), getattr(sampl, 'moment_{:d}'.format(i)))
    mom = plt.matshow(getattr(sampl, 'moment_{:d}'.format(i)))
    num = mom.figure.number
    mom.figure.canvas.set_window_title('Moment {:d} ({:d}).tif'.format(i, num))

sampl.get_median()
tifffile.imsave('Median unitcell.tif', sampl.median_unitcell.astype(np.uint16))
med = plt.matshow(sampl.median_unitcell)
num = med.figure.number
med.figure.canvas.set_window_title('Median unitcell ({:d})'.format(num))

plt.show()