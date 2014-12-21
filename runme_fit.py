from __future__ import division
import cjlib
import cjnifti
import scipy.integrate
from scipy import interpolate
from scipy.optimize import fmin_l_bfgs_b, fmin_tnc

a = cjnifti.cjnifti()

T2m = 10/1000/1000

gm_table={}
for blah in linspace(-50,50,1001):
	gm_table['%3.1f'%blah] = scipy.integrate.quad( f, 0, pi/2, args=(T2m, blah) )[0]
print 'done with table'

##================================================================
##
##  Create the function to fit to here
##
##================================================================

def wm( R1m, T2w, T2m, delta_mw, R,  RM0m_o_R1w, one_over_R1wT2w, freq, w1, f0_shift, A, delta_w):

	# passed in R1m = 1			# relaxation rate of bulk water pool (s-1)
	#R1w = 1		# relaxation rate of solid-like macromolecule pool (s-1 )
	# passed in T2w = 0.100 		# relaxation time of bulk water (s)
	#T2m = 10/1000/1000 	# relaxation time of macromolecular pool (s)
	#delta_w = 0.0		# frequency offset for bulk water pool (ppm) 
	#delta_m = 0 		# frequency offset of solid like macromolecular pool (ppm)
	# passed in delta_mw = 2.19*3*42.57	# frequency difference between solid like macromolecular pool and bulk water pool (ppm -> Hz)
	# passed in R = 41 			# Exchnage rate between bulk water pool and solid like macromolecular pool (s-1)
	#M0m = 0.14*100		# Equilibrium magnetization of solid-like macromolecule pool
	# passed in w1 = 1  *42.57*2*pi 	# RF irradtion power level (uT -> radians/sec)

	# passed in RM0m_o_R1w = 2.75
	# passed in one_over_R1wT2w = 61

	f = lambda theta,T2m,delta_m: sin(theta)*sqrt(2/pi)*T2m/abs((3*cos(theta)**2-1))*exp( -2*((2*pi*delta_m*T2m)/(abs(3*cos(theta)**2-1)))**2)

	y = zeros( freq.shape )
	Rrfm = zeros( freq.shape )

	for ii,dw in enumerate( freq ):

		delta_w = dw
		delta_m = delta_w + delta_mw

		#gm = scipy.integrate.quad( f, 0, pi/2, args=(T2m, delta_m) )[0]
		gm = gm_table[ '%3.1f' % round(delta_m/7/42.57,1) ]
		Rrfm[ii] = w1**2*pi*gm

		num = R1m * (RM0m_o_R1w) + Rrfm[ii] + R1m + R
		denom = (RM0m_o_R1w) * (Rrfm[ii] + R1m ) + ( 1 + (w1/(2*pi*delta_w))**2*one_over_R1wT2w ) * (Rrfm[ii] + R1m + R)

		y[ii] = num / denom

	ff = extrap1d( interpolate.interp1d( freq-f0_shift, y, kind='linear') )
	y = ff( freq )

	return A*y

def fun( x, freq, mm, w1, fitinds, delta_w, mm_fixed, lorentzian_fixed, scale):

	x = x * scale

	print x

	R1m, T2w, delta_mw, R,  RM0m_o_R1w, one_over_R1wT2w, f0_shift, A = x

	y = wm( R1m, T2w, T2m, delta_mw, R,  RM0m_o_R1w, one_over_R1wT2w, freq, w1, f0_shift, A, delta_w)

	figure(2)
	clf()

	subplot(2,1,1)
	plot( freq/7/42.57, mm, 'b-o', label='acquired' )
	plot( freq/7/42.57, y, 'g--', linewidth=2, label='fitted' )
	plot( freq/7/42.57, lorentzian_fixed, 'c--', linewidth=2, label='Lorentzian' )
	plot( freq[fitinds]/7/42.57, 0.01*ones(fitinds.shape), 'r.', label='acquired' )
	grid('on')
	xlim((40,-40))
	legend(loc='best')

	subplot(2,1,2)
	plot( freq/7/42.57, y-mm, 'g-' )
	plot( freq/7/42.57, lorentzian_fixed - mm_fixed, 'c--', linewidth=2, label='Lorentzian' )
	plot([-40,40], [0,0],'k')
	grid('on')
	xlim((40,-40))
	figure(2).canvas.draw()

	ie = nonzero( freq[fitinds] >= 3.5 )[0] 
	iw = nonzero( freq[fitinds] <= 3.5 )[0] 

	ae = y[fitinds[ie]] - mm[fitinds[ie]]
	aw = y[fitinds[iw]] - mm[fitinds[iw]]

	return sqrt( sum( ae**2 ) + (len(ie)/len(iw))**2* sum(aw**2) )

def extrap1d(interpolator):
	xs = interpolator.x
	ys = interpolator.y

	def pointwise(x):
		if x < xs[0]:
			return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
		elif x > xs[-1]:
			return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
		else:
			return interpolator(x)

	def ufunclike(xs):
		return array(map(pointwise, array(xs)))

	return ufunclike


##================================================================
##
##  Read in the data here
##
##================================================================


## Load the data
if 'd' not in globals():
	d = a.read('fitdata/20120105-mr3-cest-long_4_1_reg.nii.gz')

sl = d.shape[1]/2 - 3

## Load the frequencies
freq = loadtxt('fitdata/cestfreq_121.txt') / 7 / 42.57
inds_sat = find( abs( freq ) < 200 )
inds_nosat = find( abs( freq ) > 200 )

inds_pos = find( freq >= 0)
inds_neg = find( freq <= 0)

freq = freq[inds_sat]

fitinds = r_[ find( abs(freq) < 2 ), find( abs(freq)>9) ]

figure(1)
clf()
cjlib.mimagecb( d[0,sl] )
#mask = cjlib.roipoly()

mm,ss = cjlib.applymask( d[:,sl], mask )
mm = mm[inds_sat] / mean( mm[inds_nosat] ) * 0.98

##================================================================
##
##  Do fitting here
##
##================================================================

R1m = 1			# relaxation rate of bulk water pool (s-1)
T2w = 0.100 		# relaxation time of bulk water (s)
delta_mw = 2.4*7*42.57	# frequency difference between solid like macromolecular pool and bulk water pool (ppm -> Hz)
R = 41 			# Exchnage rate between bulk water pool and solid like macromolecular pool (s-1)
RM0m_o_R1w = 0.95
one_over_R1wT2w = 5.3
delta_w = 0 * 7 * 42.57
T2m = 10 / 1000 / 1000

w1 = 0.4  *42.57*2*pi 	# RF irradtion power level (uT -> radians/sec)

# Find the minimum of the data and use as the offset
indexm = find( min(mm) == mm )
f0_shift = (freq[indexm]+0.01) * 7*42.57

# extra global scaling amplitude
A = 1.0 

# Initial guess
x0 = r_[R1m, T2w, delta_mw, R, RM0m_o_R1w, one_over_R1wT2w, f0_shift, A]
scale = array( abs(x0) )

# for comparison
newfreq, mm_fixed, lorentzian_fixed, tt_A, tt_x0, tt_w, tt_b, tt_k = cjlib.cestFit( freq, mm, fitinds, freq )

# Fit the data
x = fmin_l_bfgs_b(fun, fprime=[], x0=x0/scale, args=(freq*7*42.57,mm, w1, fitinds, delta_w, mm_fixed, lorentzian_fixed, scale), bounds=( r_[0.8, 1.2]/scale[0], r_[0.07, 0.2]/scale[1], r_[-3*7*42.57, 3*7*42.57]/scale[2], r_[20,70]/scale[3], r_[0,5]/scale[4], r_[0, 100]/scale[5], r_[freq[indexm-2]*7*42.57, freq[indexm+2]*7*42.57]/scale[6], r_[0.95, 1.05]/scale[7]), approx_grad=True)[0]

# Rescale
x = x * scale 
R1m, T2w, delta_mw, R,  RM0m_o_R1w, one_over_R1wT2w, f0_shift, A = x

#  Calculate to plot
y = wm( R1m, T2w, T2m, delta_mw, R,  RM0m_o_R1w, one_over_R1wT2w, freq*7*42.57, w1, f0_shift, A, delta_w )

figure(2)
clf()

subplot(2,1,1)
plot( freq, mm, 'bo', label='acquired' )
plot( freq, y, 'g--x', label='fitted' )
plot( freq, lorentzian_fixed, 'c--x', label='fitted' )
grid('on')
xlim((40,-40))
legend()

subplot(2,1,2)
plot( freq, y-mm, 'g-x' )
plot( freq, lorentzian_fixed-mm_fixed, 'c-x' )
grid('on')
xlim((40,-40))
figure(2).canvas.draw()

