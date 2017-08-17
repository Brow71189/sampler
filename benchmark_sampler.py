#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 17:06:50 2017

@author: mittelberger2
"""
import sampler
import tifffile
import numpy as np

sampler.Sampler.c_sampler_path = '/home/mittelberger2/git/sampler/sampleUnitCell/sampleUnitCell.so'
sampl = sampler.Sampler()
sampl.load_c_sampler()
v1=np.array((956.2,1135.7))-1024
v2=np.array((1018.25,943.7))-1024
sampl.average_unitcell_shape=(32,32)
sampl.image = tifffile.imread('0000_128_nm_zoom_out_crystal.tif')
sampl.sample_rate=4
sampl.offset = np.array((0,0))
base=sampl.calculate_base_from_fft(v1,v2)
sampl.base_vec_2=base[1]*2048
sampl.base_vec_1=base[0]*2048
print('Number of unitcells: {:d}'.format(sampl.sample_image()))