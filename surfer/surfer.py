import sys
import scipy.interpolate as interpolate
import sobol_seq
import numpy as np
#from matplotlib.pyplot import *
from numpy.random import uniform, seed
from numpy import linspace, concatenate, meshgrid, isnan, array, loadtxt
from math import acos, pi
from copy import deepcopy
try:
	import qrng
except:
	pass 
	#print("Quantum RNG not available")
	
class ChaosRNG(object):

	def __init__(self, R=4, seed=None):
		self._R = R
		if seed is not None:
			self._X = seed
		else:
			self._X = uniform()

	def generate(self):
		self._X = self.logisticmap(self._X, self._R)
		return acos(1.-2*self._X)/pi

	def generate_array(self, N):
		assert(N>0)
		arr = []
		for i in range(N):
			arr.append( self.generate() )
		return array(arr)

	def logisticmap(self, x, r):
	    return x * r * (1 - x)


class QuantumRNG(object):

	def __init__(self, backend='', token=''):
		if backend=='' and token=='':
			qrng.set_provider_as_IBMQ('')
			qrng.set_backend()
		else:
			print ( " * Connecting to %s..." % backend)
			qrng.set_provider_as_IBMQ('')  ### TODO: put token
			qrng.set_backend(backend)
			print ( " * Connection completed")

	def generate(self):
		res = qrng.get_random_double(0,1)
		#print ( " * Measurement performed! (%f)" % res)
		return res

	def generate_array(self, N):
		assert(N>0)
		print ( " * Performing %d measurements of states in superposition, please wait..." % N)
		arr =[]
		for i in range(N):
			arr.append( self.generate() )
		return array(arr) 




class surF(object):

	def __init__(self):
		self._fitness_function = None
		self._fun_samples = None
		self._boundaries = None
		self._result_interpolation = None
		self._grid = None
		self._last_coefficients = None

		self.fig = figure()
		self.ax = self.fig.add_subplot(1, 1, 1, projection="3d")


	def specify_fitness(self, fun):
		self._fitness_function = fun

	def specify_search_space(self, bounds):
		self._boundaries = np.array(bounds)

	def build_model(self, numpoints=100, resolution=100, coefficients=10,
		sampling_method="PRNG"):
		print (" * Taking %d samples of the search space" % numpoints) 
		if sampling_method=="PRNG": 
			print ("   using preudo-random random numbers")
		elif sampling_method=="QRNG":
			print ("   using quasi-random random numbers")
		elif sampling_method=="CRNG": 
			print ("   using chaotic random numbers")
		elif sampling_method=="PPRG":
			print ("   using point packing")
		elif sampling_method=="QUANTUM":
			print ("   using quantum-generated random numbers")

		self._sample_from_fun(points=numpoints, sampling_method=sampling_method)

		# create grid
		aggr = [ linspace(x[0], x[1], resolution) for x in self._boundaries ]
		self._grid = np.array(meshgrid( *aggr , indexing='ij'))
		
		self._calc_approx_from_points()
		print (" * Creating interpolation out of samples, using resolution %d points" % resolution)
		self._interpolate()
		print (" * Building the smoothed surrogate model using %d coefficients" % coefficients)
		self._create_model(coeffs=coefficients)

		self._last_coefficients = coefficients

	def _interpolate(self):
		self._result_interpolation = interpol(self._fun_samples, self._f_approx_values, tuple([g for g in self._grid]))

	def _sample_from_fun(self, points=100, sampling_method="PRNG"):

		D = len(self._boundaries)

		if sampling_method=="PRNG":
			# random uniform
			hyperrand = np.array([uniform(self._boundaries[d][0], self._boundaries[d][1], points) for d in range(D)]).T
			self._fun_samples = hyperrand.reshape((points, D))
		
		elif sampling_method=="QRNG":
			# random Sobol
			sobol = sobol_seq.i4_sobol_generate(D, points)
			sobol = self._boundaries.T[0] + sobol * (self._boundaries.T[1] - self._boundaries.T[0])
			self._fun_samples = np.array(sobol).reshape((points, D))
		
		elif sampling_method=="CRNG":
			C = ChaosRNG()
			all_points = []
			for d in range(D):
				p = C.generate_array(points)
				all_points.append(p)
			all_points = array(all_points)
			fignew, axnew = subplots(1,1)
			res = self._boundaries.T[0] + all_points.T * (self._boundaries.T[1] - self._boundaries.T[0])
			self._fun_samples = deepcopy(res)

		elif sampling_method=="PPRG":
			raise Exception("Generation method not implemented: "+sampling_method)
		
		elif sampling_method=="QUANTUM":
			Q = QuantumRNG( '', '' ) # TODO: implement auth interface on Q experience
			all_points = []
			for d in range(D):
				p = Q.generate_array(points)
				all_points.append(p)
			all_points = array(all_points)
			fignew, axnew = subplots(1,1)
			res = self._boundaries.T[0] + all_points.T * (self._boundaries.T[1] - self._boundaries.T[0])
			self._fun_samples = deepcopy(res)


		else:
			raise Exception("Generation method not supported: "+sampling_method)

	def _calc_approx_from_points(self, use_multiproc=False):
		from multiprocessing import Pool
		p = Pool(processes=4)
		if use_multiproc:
			all_evaluations = p.map(self._fitness_function, self._fun_samples) 
			print (all_evaluations)
			#exit()
		else:
			all_evaluations = [self._fitness_function(v) for v in self._fun_samples]
		self._f_approx_values = np.array(all_evaluations)
		
	def _create_model(self, coeffs=10):
		self._f_approx = approximate_fitn(
			self._result_interpolation,
			self._boundaries,
			n_coeffs=coeffs,
		)

	def approximate(self, v):
		return self._f_approx(v)		

	def update_model(self, new_samples, fitness_values):
		print (" * Updating model with %d new samples" % len(new_samples))
		self._fun_samples = concatenate((self._fun_samples, new_samples))
		self._f_approx_values = concatenate((self._f_approx_values, fitness_values))
		self._interpolate()
		self._create_model(coeffs=self._last_coefficients)

	def plot_approximate_landscape(self):
		obj_fft_approx_values = array([self.approximate(tu) for tu in zip(*self._grid)])
		self.ax.plot_surface(self._grid[0], self._grid[1], obj_fft_approx_values, label="Surrogate model", alpha=0.6)
		#show()

	def plot_real_landscape(self, res=100):
		res_true_fit = []
		for x in linspace(self._boundaries[0][0], self._boundaries[0][1], res):
			for y in linspace(self._boundaries[1][0], self._boundaries[1][1], res):
				res_true_fit.append(self._fitness_function([x,y]))
		res_true_fit = array(res_true_fit).reshape((res,res))
		self.ax.plot_surface(self._grid[0], self._grid[1], res_true_fit, label="Real function", alpha=0.6)
		#show()


def interpol(pnts, f_app_values, lattice):
	"Use linear interpolation where possible and nearest interpolation otherwise"
	z1 = interpolate.griddata(pnts, f_app_values, lattice, method="linear")
	z2 = interpolate.griddata(pnts, f_app_values, lattice, method="nearest")
	z = np.where(isnan(z1), z2, z1)
	return z

def approximate_fitn(f_values, space_size, n_coeffs=None):
	n = len(f_values)
	d = len(space_size)
	if n_coeffs is None:
		rem_coeffs = n
	else:
		rem_coeffs = n - n_coeffs
	coeffs = np.fft.fftn(f_values)
	coeffs = np.fft.fftshift(coeffs)

	for idx_tuple in np.ndindex(coeffs.shape):
		for i in idx_tuple: 
			if i < rem_coeffs//2 or i >= n - rem_coeffs//2:
				coeffs[idx_tuple] = 0

	coeffs = np.fft.ifftshift(coeffs)
	f_approx_vals = np.fft.ifftn(coeffs).real

	points = []
	for i in range(0, d):
		points.append(np.linspace(space_size[i][0], space_size[i][1], n))

	interp = interpolate.RegularGridInterpolator(tuple(points), f_approx_vals, method='linear')

	def f_approx(coords):
		return interp(tuple(coords))

	return f_approx



if __name__ == '__main__':
	pass