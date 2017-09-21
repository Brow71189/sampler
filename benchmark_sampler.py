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
import sys

sampler.Sampler.c_sampler_path = './sampleUnitCell.so'
sampl = sampler.Sampler()
sampl.load_c_sampler()


img = sampler.calculate_counts(tifffile.imread('0000_SuperScan (MAADF).tif'))


sampl.base_vec_1 = np.array((-40.97, -17.71))
sampl.base_vec_2 = np.array((-26.88, 5.68))


sampl.average_unitcell_shape=(32,32)
sampl.sample_rate=128

'''
match_map = np.zeros((61,61),dtype=np.float32)
for y in range(0, 61):
    top = 64 * y
    for x in range(0, 61):
        left = 64 * x
        bottom = top + 256
        right = left + 256
        sampl.image = img[top:bottom,left:right]
        num_cells = sampl.sample_image()
        contrast = np.std(sampl.average_unitcell)/np.mean(sampl.average_unitcell)
        match_map[y,x] = contrast
        sys.stdout.write('.')
        if x == 60:
            print(' y: {:d}'.format(y))
        else:
            sys.stdout.flush()
        #print('top: {:d},left: {:d} Number of unitcells: {:d} contrast: {:f} '.format(top,left,num_cells,contrast))

tifffile.imsave('fine_match_map.tif', match_map.astype(np.float32))

sys.exit(0)
'''
'''    
left = 24*128
right = 26*128+512

top = 17*128 
bottom =20*128 + 512
'''

left = 3072
right = 4096
top = 2048
bottom = 3072

'''
left = 25*128
right = 25*128+512
top = 19*128 
bottom =19*128 + 512
'''
sampl.image = img[top:bottom,left:right]
num_cells = sampl.sample_image()
contrast = np.std(sampl.average_unitcell)/np.mean(sampl.average_unitcell)
print('top: {:d},left: {:d},bottom: {:d}, right: {:d} Number of unitcells: {:d} contrast: {:f} '.format(top,left,bottom,right,num_cells,contrast))
sampl.periodic_repeats=8

plt.gray()
plt.close()

tifffile.imsave('intSuperScan.tif', sampl.image.astype(np.uint16))

tifffile.imsave('Average unitcell.tif', sampl.average_unitcell.astype(np.float32))


for i in range(0, 3):
    sampl.make_pretty_output(i)
    tifffile.imsave('Pretty_Moment_{:d}'.format(i), sampl.pretty_unitcell.astype(np.float32))
'''    
    mom = plt.matshow(sampl.pretty_unitcell)
    num = mom.figure.number
    mom.figure.canvas.set_window_title('Pretty Moment {:d} ({:d}).tif'.format(i, num))
'''
plt.show()
