#!/usr/bin/env python

import sys
import os
import shutil
import math
import numpy as np
from nanoplatelet import Nanoplatelet
import matplotlib.pyplot as plt

#helix_length = 1040
xtal_length = 220
xtal_width = 70
xtal_thickness = None
xtal_nlayers = 11

radius = 110 #xtal_length/(2*math.pi)
pitch = 220

#delta = xtal_width*math.sqrt(pitch**2+4*np.pi**2*radius**2)/(2*np.pi*radius)
#xtal_length = (helix_length-delta)*math.sqrt(
#        pitch**2+4*np.pi**2*radius**2)/pitch
#print(delta, xtal_length)

npl = Nanoplatelet()

npl.set_xtal_unit_cell('CdSe', 'zincblende', 6.08, 6.08, 6.08)

#TODO: Check angle theta
npl.set_xtal_extents(xtal_length, xtal_width, xtal_thickness, bc='nnn', 
                     nlayers=xtal_nlayers)

npl.set_shape('Plane')
#npl.set_shape('Helicoid', pitch=pitch)
#npl.set_shape('Cylinder', radius=radius)
#npl.set_shape('HelicalRibbon', radius=radius, pitch=pitch)
#npl.set_shape('Helicoid', pitch=pitch)
#npl.set_shape('HelicalRibbon', radius=radius, pitch=pitch)
#npl.set_shape('General', radius=radius, pitch=pitch, profile='Line')

#npl.create(orient_along=[0,0,1])
npl.create(orient_along=None)
print(f"R = {npl.radius:g}")
print(f"P = {npl.pitch:g}")
print(f"kg = {npl.kg:g}")
print(f"km = {npl.km:g}")
print(f"theta = {npl.theta:g}")

fn = f"ml5-plane.xyz"
npl.to_xyz(fn)
num_atoms = len(npl._atoms)
print(f"Number of atoms = {num_atoms}")

fn = f"ml5-plane-aesdb.csv"
Q, pattern = npl.calc_scattering_pattern(
                calculator='AESDebye', fn=None, pattern_type='x', beg=1.0,
                end=8.0, step=0.002, nthreads=-1, ncells=15,
                debyer_cmd='debyer')
#print(Q, pattern)
npl_vol = npl.xtal_length * npl.xtal_width * npl.xtal_thickness
pattern *= (0.001/(num_atoms*npl_vol*1e-3))
plt.plot(Q*10, pattern, ls='None', marker='x', label='AESDebye')

fn = f"ml5-plane-dbr.txt"
Q, pattern = npl.calc_scattering_pattern(
                calculator='Debyer', fn=None, pattern_type='x', beg=1.0,
                end=8.0, step=0.002, nthreads=-1, ncells=15,
                debyer_cmd='debyer')
#print(Q, pattern)
npl_vol = npl.xtal_length * npl.xtal_width * npl.xtal_thickness
pattern *= (0.001/(npl_vol*1e-3))
plt.plot(Q*10, pattern, label='Debyer')

plt.legend(loc='best')
plt.xlabel('Q (nm$^{-1}$)')
plt.show()
