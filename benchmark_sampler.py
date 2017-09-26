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
import math

sampler.Sampler.c_sampler_path = './sampleUnitCell.so'
sampl = sampler.Sampler()
sampl.load_c_sampler()


img = sampler.calculate_counts(tifffile.imread('0000_SuperScan (MAADF).tif'))

#img = tifffile.imread('FullSuperScanIntM.tif').astype(np.int32)

#sampl.base_vec_1 = np.array((-40.97, -17.71))
a2 = sampl.base_vec_2 = np.array((-26.88, 5.68))
a1 = sampl.base_vec_1 = np.array((14.09,23.39))

#sampl.offset = 49.5/54*a1 + 34/54*a2 + np.array((-0.50,-0.25))


sampl.offset = np.array((536.949722+2048,521.195463+3072))

r1 = (2*(a1[0]**2+a1[1]**2)**0.5)
r2 = (2*(a2[0]**2+a2[1]**2)**0.5)

w1 = math.atan2(a1[0],a1[1])
w2 = math.atan2(a2[0],a2[1])

len1 = math.ceil(r1)
len2 = math.ceil(r2)

sampl.average_unitcell_shape=(len1,len2)
sampl.sample_rate=256

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
sampl.offset-=np.array ((top,left))
offx = sampl.offset[1]
offy = sampl.offset[0]
print('a1(r,w): ({:f},{:f})\ta2(r,w): ({:f},{:f})'.format(r1,w1,r2,w2))
print('a1(x,y): ({:f},{:f})\ta2(x,y):({:f},{:f})\tuc: {:d}x{:d}'.format(a1[1],a1[0],a2[1],a2[0],len1,len2))
print('offsetx: {:f} from {:d}\toffsety: {:f} from {:d}'.format(offx,left,offy,top))
print('top: {:d},left: {:d},bottom: {:d}, right: {:d}'.format(top,left,bottom,right)) 
num_cells = sampl.sample_image()
contrast = np.std(sampl.average_unitcell)/np.mean(sampl.average_unitcell)
print('Number of unitcells: {:d} contrast: {:f} '.format(num_cells,contrast))
noff = sampl.get_internal_offset()
print("suggesting offsets: {:f},{:f}".format(noff[1],noff[0]))

tifffile.imsave('intSuperScan.tif', sampl.image.astype(np.uint16))

tifffile.imsave('Average unitcell.tif', sampl.average_unitcell.astype(np.float32))

for i in range(1, 2):
    sampl.make_pretty_output(i,True)    
    pretty=sampl.pretty_unitcell.copy()     
    tifffile.imsave('Pretty_Moment_{:d}'.format(i), sampl.pretty_unitcell.astype(np.float32))
    sampl.make_pretty_output(i,False)
    ugly=sampl.pretty_unitcell.copy()  
    tifffile.imsave('Ugly_Moment_{:d}'.format(i), sampl.pretty_unitcell.astype(np.float32))
    diff = pretty - ugly
    tifffile.imsave('Diff_Moment_{:d}'.format(i), diff.astype(np.float32))    
'''    
    mom = plt.matshow(sampl.pretty_unitcell)
    num = mom.figure.number
    mom.figure.canvas.set_window_title('Pretty Moment {:d} ({:d}).tif'.format(i, num))
'''

