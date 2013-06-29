############################################################################
# File         : cylin_projection.py
# Author       : Neil Massey
# Created      : 08/03/12
# Purpose      : Simple cylindrical	projection of data
#                
############################################################################

import numpy

from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib.ticker import FixedLocator, Formatter
from matplotlib.transforms import BboxTransformTo, Transform, Affine2D, Bbox, \
	 IdentityTransform, TransformWrapper
from proj_common import *

############################################################################

class CylinAxes(Axes):
	"""
		Cylindrical projection
	"""
	name = "cylindrical"

	class CylinTransform(Transform):
		input_dims = 2
		output_dims = 2
		is_separable = False

		def transform(self, ll):
			return ll

		def transform_path(self, path):
			vertices = path.vertices
			ipath = path.interpolated(path._interpolation_steps)
			return Path(self.transform(ipath.vertices), ipath.codes)

		def inverted(self):
			return CylinAxes.CylinTransform()
	
	########################################################################

	def __init__(self, *args, **kwargs):
		Axes.__init__(self, *args, **kwargs)
		self.set_aspect(1.0, adjustable='box', anchor = 'C')
		self.cla()
		
	########################################################################

	def cla(self):
		"""
			Override to set up sensible defaults
		"""
		Axes.cla(self)
		self.set_longitude_grid(45)
		self.set_latitude_grid(20)
		# Do not display ticks -- we only want gridlines and text
		self.xaxis.set_ticks_position('none')
		self.yaxis.set_ticks_position('none')

		self.x_lim = [-180, 180]
		self.y_lim = [-90,  90]
		self.set_xlim(self.x_lim)
		self.set_ylim(self.y_lim)

	########################################################################

	def set_longitude_grid(self, degrees):
		"""
            Set the number of degrees between each longitude grid.
		"""
		# Set up a FixedLocator at each of the points, evenly spaced
		# by degrees.
		number = (360 / degrees) + 1
		self.lon_grid_deg = degrees
		self.xaxis.set_major_locator(FixedLocator(
			numpy.linspace(-180, 180, number, False)))
		self.xaxis.set_major_formatter(EastWestDegreeFormatter(degrees))

	########################################################################

	def set_latitude_grid(self, degrees):
		"""
			Set the number of degrees between each longitude grid.
		"""
		# Set up a FixedLocator at each of the points, evenly spaced
		# by degrees.
		number = (160 / degrees) + 1
		self.lat_grid_deg = degrees
		self.yaxis.set_major_locator(FixedLocator(
			numpy.linspace(-80, 80, number, True)))
		self.yaxis.set_major_formatter(NorthSouthDegreeFormatter(degrees))
