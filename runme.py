from __future__ import division
from scipy.io import loadmat
from numpy import *
from matplotlib.pyplot import *
from scipy.ndimage.filters import median_filter

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.font_manager import FontProperties
import matplotlib.pyplot as plt

import cjlib
import bloch

# Water and MT data taken from Seth's 2006 MRM paper for 3T WM

##  protons
amide = bloch.Pool('amide', t1 = 1.7, t2 = 0.033, lifetime = 1/28, chemical_shift = 3.5, concentration = 71.9)
mm = bloch.Pool('macromolecular', t1 = 1.7, t2 = 0.00001, lifetime = 1/4, chemical_shift = -2.43, concentration = 2*55000*0.13)
freewater = bloch.Pool('freewater', t1 = 1.700, t2 = 0.080, lifetime = 0, chemical_shift = 0.0, concentration = 2*55000)

#freqdata = linspace(-20,20,201)
freqdata = loadtxt('/Users/craig/work/Research/Data/20120209-mr3-cest-many-ik/data/cestfreq_121.txt')
freqdata = freqdata[3:-2] / 7 / 42.57

#================================================
##
## Setup the Pulse
##
##================================================

sg_100_100_0_base = array([[0, 413, 859, 1337, 1847, 2389, 2962, 3568, 4204, 4870,
	5565, 6288, 7038, 7814, 8613, 9435, 10277, 11138, 12015, 12906,
	13809, 14721, 15640, 16563, 17487, 18409, 19327, 20237, 21136, 22023,
	22892, 23742, 24570, 25372, 26146, 26889, 27598, 28271, 28905, 29497,
	30046, 30549, 31005, 31411, 31767, 32070, 32320, 32515, 32655, 32739,
	32767, 32739, 32655, 32515, 32320, 32070, 31767, 31411, 31005, 30549,
	30046, 29497, 28905, 28271, 27598, 26889, 26146, 25372, 24570, 23742,
	22892, 22023, 21136, 20237, 19327, 18409, 17487, 16563, 15640, 14721,
	13809, 12906, 12015, 11138, 10277, 9435, 8613, 7814, 7038, 6288,
	5565, 4870, 4204, 3568, 2962, 2389, 1847, 1337, 859, 413,
	0 ]]).transpose()

duration = 25./1000.

s_sat = zeros( (len(freqdata), N ) )

b1 = 0.6 * 1 * 42.57 # ( Hz )

sg_100_100_0 = concatenate( ( ones(sg_100_100_0_base.shape), sg_100_100_0_base / sg_100_100_0_base.max(), ones(sg_100_100_0_base.shape) ), axis=1 )
sg_100_100_0[:,1] = sg_100_100_0[:,1] * b1;
sg_100_100_0[:,2] = sg_100_100_0[:,2] * duration/len(sg_100_100_0_base) # in s
sg_100_100_0 = concatenate( (sg_100_100_0, array([[1, 0, 45/1000]]) ) )

pulse = sg_100_100_0

hard_pulse = zeros((1,3))
hard_pulse[0,0] = 1
hard_pulse[0,1] = 0.4 * 42.57
hard_pulse[0,2] = 2

##================================================
##
##  Loop over the delays 
##
##================================================

s_sat = bloch.solve((freewater, amide, mm), freqdata, pulse, magnetic_field=11.7, verbose=1, pulse_repeat=200, post_dynamic_delay=0)

inds = find( abs(freqdata) <= 11 )
freq = freqdata[ inds ]

fitinds = r_[ find( abs(freq) <= 1 ), find( abs(freq) >= 5 ) ] 

newfreq, mm_fixed, lorentzian_fixed, A, x0, w, b, k = cjlib.cestFit( freq, s_sat[inds,-1], fitinds, newfreq=freq )

figure(1)
clf()
plot( freqdata, s_sat[:,-1] )
plot( freq, lorentzian_fixed, 'g--' )
xlabel('Frequency Offset From Water (ppm)')
ylabel('Relative Intensity (%)')
grid('on')
xlim((10,-10))
ylim((0,1.1))
#text(5.5, 0.6, freewater.textify(), fontsize=12)
#text(5.5, 0.5, mm.textify(), fontsize=12)
#text(5.5, 0.4, amide.textify(), fontsize=12)
figure(1).canvas.draw()

figure(2)
clf()
plot( freq, lorentzian_fixed - mm_fixed )
plot( [-10,10], [0,0], 'k', linewidth=2)
xlabel('Frequency Offset From Water (ppm)')
ylabel('Relative Intensity (%)')
grid('on')
xlim((10,-10))
ylim((-0.2,0.2))
figure(2).canvas.draw()

