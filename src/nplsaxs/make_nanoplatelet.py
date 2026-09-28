#!/usr/bin/env python

import sys
import os
import shutil
import math
import numpy as np
from nanoplatelet import Nanoplatelet
import matplotlib.pyplot as plt

#helix_length = 1040
xtal_length = 1000
xtal_width = 100
xtal_thickness = None
xtal_nlayers = 7

radius = 102.41
pitch = 125.25

#delta = xtal_width*math.sqrt(pitch**2+4*np.pi**2*radius**2)/(2*np.pi*radius)
#xtal_length = (helix_length-delta)*math.sqrt(
#        pitch**2+4*np.pi**2*radius**2)/pitch
#print(delta, xtal_length)

npl = Nanoplatelet()

npl.set_xtal_unit_cell('CdSe', 'zincblende', 6.08, 6.08, 6.08)

#TODO: Check angle theta
npl.set_xtal_extents(xtal_length, xtal_width, xtal_thickness, bc='nnn', 
                     nlayers=xtal_nlayers)

#npl.set_shape('Plane')
#npl.set_shape('Helicoid', pitch=pitch)
#npl.set_shape('Helicoid', pitch=pitch)
npl.set_shape('HelicalRibbon', radius=radius, pitch=pitch)
#npl.set_shape('General', radius=radius, pitch=pitch, profile='Line')

#npl.create(orient_along=[0,0,1])
npl.create(orient_along=None)

print(f"R = {npl.radius:g}")
print(f"P = {npl.pitch:g}")
print(f"kg = {npl.kg:g}")
print(f"km = {npl.km:g}")
print(f"theta = {npl.theta:g}")

#fn = f"ml3-test_inf.xyz"
#fn = f"ml2-test_hc_{radius/10:.0f}_{pitch/10:.0f}.xyz"
#npl.to_xyz(fn)


#fn = f"ml3-test_inf.txt"
#fn = f"ml2-test_hc_{radius/10:.0f}_{pitch/10:.0f}.txt"
Q, pattern = npl.calc_scattering_pattern(
                calculator='AESDebye', fn='test.csv', pattern_type='x', beg=0.001,
                end=10.0, step=0.002, nthreads=-1, ncells=15,
                debyer_cmd='debyer')
#print(Q, pattern)
#npl_vol = npl.xtal_length * npl.xtal_width * npl.xtal_thickness
#pattern *= (0.001/(npl_vol*1e-3))
#plt.loglog(Q*10, pattern)
#plt.show()
