############################################################################
# File         : robin_projection.py
# Author       : Neil Massey
# Created      : 08/03/12
# Purpose      : Robinson projection of data
#                
############################################################################

import numpy

from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib import spines
from matplotlib import axis
from matplotlib.ticker import FixedLocator, Formatter
from matplotlib.transforms import BboxTransformTo, Transform, Affine2D, Bbox, \
	 IdentityTransform, TransformWrapper
from matplotlib.patches import Circle

from proj_common import *
import numpy

############################################################################

def plen_pdfe_plat_arrays():
	# get the plen and pdfe from the table
	plat = numpy.arange(-95, 100, 5)
	plen = numpy.array([0.5322,  0.5322,  0.5722, 0.6213, 0.6732, 0.7186, 0.7597, 
						0.7986,  0.8350, 0.8679, 0.8962, 0.9216, 0.9427, 
						0.9600,  0.9730, 0.9822, 0.9900, 0.9954, 0.9986,
						1.0000,  0.9986, 0.9954, 0.9900, 0.9822, 0.9730,
					    0.9600,  0.9427, 0.9216, 0.8962, 0.8679, 0.8350,
					    0.7986,  0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 
						0.5322,  0.5322])
	# This table has been translated from the original (on wikipedia) table
	# by (pdfe + 1) * 0.5
	pdfe = numpy.array([-0.05,   0.00000, 0.01195, 0.03030, 0.05320, 0.07825, 0.10485,
					    0.13270, 0.16155, 0.19120, 0.22145, 0.25210, 0.28300,
						0.31400, 0.34500, 0.37600, 0.40700, 0.43800, 0.46900,
						0.50000, 0.53100, 0.56200, 0.59300, 0.62400, 0.65500,
						0.68600, 0.71700, 0.74790, 0.77855, 0.80880, 0.83845,
						0.86730, 0.89515, 0.92175, 0.94680, 0.96970, 0.98805,
						1.00000, 1.05])
	return plen, pdfe, plat

############################################################################

def get_plen_pdfe(lats):
	plen, pdfe, plat = plen_pdfe_plat_arrays()
	plen_1 = numpy.interp(lats, plat, plen)
	pdfe_1 = numpy.interp(lats, plat, pdfe)
	return plen_1, pdfe_1

############################################################################

def get_lat_len_from_pdfe(in_pdfe):
	plen, pdfe, plat = plen_pdfe_plat_arrays()
	# do the interpolation backwards
	out_lat = numpy.interp(in_pdfe, pdfe, plat)
	out_plen = numpy.interp(in_pdfe, pdfe, plen)
	return out_lat, out_plen

############################################################################

class RobinAxes(Axes):
	"""
		Robinson projection
	"""
	name = "robinson"

	class RobinTransform(Transform):
		input_dims = 2
		output_dims = 2
		is_separable = False

		def transform(self, ll):
			x = ll[:, 0:1]
			y = ll[:, 1:2]
			plen, pdfe = get_plen_pdfe(y)
			x_proj = x*(1.0/360)*plen + 0.5
			y_proj = pdfe
			xy = numpy.concatenate((x_proj, y_proj),1)
			return xy

		def transform_path(self, path):
			vertices = path.vertices
			ipath = path.interpolated(path._interpolation_steps)
			return Path(self.transform(ipath.vertices), ipath.codes)

		def inverted(self):
			return RobinAxes.RobinInverseTransform()

	########################################################################

	class RobinInverseTransform(Transform):
		input_dims = 2
		output_dims = 2
		is_separable = False

		def transform(self, ll):
			x = ll[:, 0:1]
			y = ll[:, 1:2]
			x_proj = []
			y_proj = []
			for y0 in range(0, y.shape[0]):
				yp, plen = get_lat_len_from_pdfe(y[y0][0])
				xp = (x[y0][0] - 0.5)/plen * 360
				y_proj.append([yp])
				x_proj.append([xp])
			xy = numpy.concatenate((x_proj, y_proj),1)
			return xy

		def transform_path(self, path):
			vertices = path.vertices
			ipath = path.interpolated(path._interpolation_steps)
			return Path(self.transform(ipath.vertices), ipath.codes)

		def inverted(self):
			return RobinAxes.RobinTransform()
	
	########################################################################

	def __init__(self, *args, **kwargs):
		Axes.__init__(self, *args, **kwargs)
		self.set_aspect(1.0, adjustable='box', anchor = 'C')
		self.cla()

	########################################################################

	def _init_axis(self):
		self.xaxis = axis.XAxis(self)
		self.yaxis = axis.YAxis(self)
		self._update_transScale()
		
	########################################################################

	def cla(self):
		"""
			Override to set up sensible defaults
		"""
		Axes.cla(self)
		self.set_longitude_grid(45)
		self.set_latitude_grid(20)
		self.set_xlim([-180, 180])
		self.set_ylim([-90, 90])
		# Do not display ticks -- we only want gridlines and text
		self.xaxis.set_ticks_position('none')
		self.yaxis.set_ticks_position('none')

	########################################################################

	def set_xlim(self, *args, **kwargs):
		"""
			Override set_xlim and set_ylim to be between -180,180 and -90,90
		"""
		Axes.set_xlim(self, -180, 180)
		Axes.set_ylim(self, -90, 90)

	def set_ylim(self, *args, **kwargs):
		self.set_xlim()

	########################################################################

	def set_longitude_grid(self, degrees):
		"""
            Set the number of degrees between each longitude grid.
		"""
		# Set up a FixedLocator at each of the points, evenly spaced
		# by degrees.
		number = (360.0 / degrees) + 1
		self.xaxis.set_major_locator(FixedLocator(
			numpy.linspace(-180, 180, number, True)))
		self.xaxis.set_major_formatter(EastWestDegreeFormatter(degrees))

	########################################################################

	def set_latitude_grid(self, degrees):
		"""
			Set the number of degrees between each longitude grid.
		"""
		# Set up a FixedLocator at each of the points, evenly spaced
		# by degrees.
		number = (160.0 / degrees) + 1
		self.yaxis.set_major_locator(FixedLocator(
			numpy.linspace(-80, 80, number, True)))
		self.yaxis.set_major_formatter(NorthSouthDegreeFormatter(degrees))

	########################################################################

	def _set_lim_and_transforms(self):
		self.transProjection = self.RobinTransform()
		self.transAxes = BboxTransformTo(self.bbox)
		self.transData = (
			self.transProjection +
			self.transAxes )

		self._xaxis_transform = (
			# projection works in lat / lon space so translate to this
			Affine2D().translate(0.0,-0.5).scale(1.0,180.0) +			
			self.transData )

		self._xaxis_text1_transform = (
			Affine2D().translate(0.0, -90) + 
			self.transData +
			Affine2D().translate(0.0, -50.0) )	# what are these coordinates in?

		self._yaxis_transform = (
			# projection works in lat / lon space so translate to this
			Affine2D().translate(-0.5,0.0).scale(360.0,1.0) +
			self.transData)

		self._yaxis_text1_transform = (
			Affine2D().translate(-180,0.0) +
			self.transData +
			Affine2D().translate(-50.0,0.0))

	########################################################################

	def get_xaxis_transform(self,which='grid'):
		assert which in ['tick1','tick2','grid']
		return self._xaxis_transform

	def get_xaxis_text1_transform(self, pixelPad):
		return self._xaxis_text1_transform, 'bottom', 'center'

	def get_yaxis_transform(self,which='grid'):
		assert which in ['tick1','tick2','grid']
		return self._yaxis_transform

	def get_yaxis_text1_transform(self, pixelPad):
		return self._yaxis_text1_transform, 'center', 'right'

	########################################################################

	def _gen_axes_spines(self):
		lats = numpy.arange(-90, 90, 1)
		plen, pdfe = get_plen_pdfe(lats)
		path = []
		# draw top
		xp = 0.5 - 0.5 * plen[0]
		yp = pdfe[0]
		path.append((xp, yp))
		# draw right hand side
		for l in range(0, lats.shape[0]):
			xp = 0.5*plen[l] + 0.5
			yp = pdfe[l]
			path.append((xp, yp))
		# draw left hand side
		for l in range(lats.shape[0]-1, -1, -1):
			xp = 0.5 - 0.5 * plen[l]
			yp = pdfe[l]
			path.append((xp, yp))
		axes_path = Path(path)
		return {'custom_robin':spines.Spine(self, 'linear', axes_path)}
