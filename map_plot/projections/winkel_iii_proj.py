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

def winkel_iii_proj(lon, lat, max_x = 1e20, max_y = 1e20):
	# winkel projection taken from Wikipedia:
	# http://en.wikipedia.org/wiki/Winkel_tripel_projection
	lam = numpy.radians(lon)
	hlam = lam * 0.5
	sin_hlam = numpy.sin(hlam)
	phi = numpy.radians(lat)
	cos_phi = numpy.cos(phi)

	pi_div_2 = 0.636619772
	phi_1 = numpy.arccos(pi_div_2)
	cos_phi_1 = pi_div_2

	alpha = numpy.arccos(cos_phi_1*numpy.cos(hlam))
	# need to use unnormalised sinc - numpy.sinc returns the normalised version
	# sinc(x) = sin(x) / x
	sinc_alpha = numpy.sin(alpha) / alpha

	x = (lam * cos_phi_1 + (2*cos_phi*sin_hlam/sinc_alpha))
	y = (phi + (numpy.sin(phi)/sinc_alpha))

	# do a transform to 0..1 if coordinates are supplied
	if max_x != 1e20:
		x = x/max_x * 0.5 + 0.5
	if max_y != 1e20:
		y = y/max_y * 0.5 + 0.5

	lat_idxs = numpy.where(lat < -90)
	if len(lat_idxs) != 0.0:
		y[lat_idxs] = 0.0

	lat_idxs = numpy.where(lat > 90)
	if len(lat_idxs) != 0.0:
		y[lat_idxs] = 1.0

	return x, y

############################################################################

class WinkelIIIAxes(Axes):
	"""
		Winkel-III projection
	"""
	name = "winkel-iii"

	class WinkelTransform(Transform):
		input_dims = 2
		output_dims = 2
		is_separable = False
	
		def __init__(self, *args, **kwargs):
			Transform.__init__(self, *args, **kwargs)
			self.x_max, dmax_y = winkel_iii_proj(180, 0)
			dmax_x, self.y_max = winkel_iii_proj(180, 90)

		def transform(self, ll):
			x = ll[:, 0:1]
			y = ll[:, 1:2]
			x_proj, y_proj = winkel_iii_proj(x, y, self.x_max, self.y_max)
			xy = numpy.concatenate((x_proj, y_proj),1)
			return xy

		def transform_path(self, path):
			vertices = path.vertices
			ipath = path.interpolated(path._interpolation_steps)
			return Path(self.transform(ipath.vertices), ipath.codes)
	
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
#		self.xaxis.set_major_formatter(EastWestDegreeFormatter(degrees))

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
		self.transProjection = self.WinkelTransform()
		self.transAxes = BboxTransformTo(self.bbox)
		self.transData = (
			self.transProjection +
			self.transAxes )

		self._xaxis_transform = (
			# projection works in lat / lon space so translate to this
			Affine2D().translate(0.0,-0.5).scale(1.0,180.0) +			
			self.transData )

		self._xaxis_text1_transform = (
			Affine2D().translate(0.0, -90.1) + 
			self.transData +
			Affine2D().translate(0.0, -50) )

		self._yaxis_transform = (
			# projection works in lat / lon space so translate to this
			Affine2D().translate(-0.5,0.0).scale(360.0,1.0) +
			self.transData )

		self._yaxis_text1_transform = (
			Affine2D().translate(-180, 0.0) +
			self.transData +
			Affine2D().translate(-50, 0.0) )

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
		lons = numpy.arange(-180, 180, 1)
		x_max, dmax_y = winkel_iii_proj(180, 0)
		dmax_x, y_max = winkel_iii_proj(180, 90)

		path = []
		# draw bottom
		for l in range (0, lons.shape[0]):
			xp, yp = winkel_iii_proj(lons[l], lats[0], x_max, y_max)
			path.append((xp, yp))
		# draw right hand side
		for l in range(0, lats.shape[0]):
			xp, yp = winkel_iii_proj(lons[-1], lats[l], x_max, y_max)
			path.append((xp, yp))
		# draw the bottom
		for l in range(lons.shape[0]-1, -1, -1):
			xp, yp = winkel_iii_proj(lons[l], lats[-1], x_max, y_max)
			path.append((xp, yp))
		# draw left hand side
		for l in range(lats.shape[0]-1, -1, -1):
			xp, yp = winkel_iii_proj(lons[0], lats[l], x_max, y_max)
			path.append((xp, yp))
		axes_path = Path(path)
		return {'custom_winkel_iii':spines.Spine(self, 'linear', axes_path)}
