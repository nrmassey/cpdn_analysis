############################################################################
# File      : north_polar_proj.py
# Author    : Neil Massey
# Created   : 08/03/12
# Purpose   : Polar projection for the north pole only
# Reference : http://stackoverflow.com/questions/2417794/how-to-make-the-angles-in-a-matplotlib-$
############################################################################

import numpy

from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D, Bbox, IdentityTransform
from proj_common import *

############################################################################

class NorthPolarAxes(PolarAxes):
	"""Stereographic polar projection for North Pole only.  
	   Longitude, latitude coordinates"""
	name = 'northpolar'
	class NorthPolarAxisTransform(PolarAxes.PolarTransform):
		def transform(self, tr):
			xy	 = numpy.zeros(tr.shape, numpy.float_)
			t    = tr[:, 0:1] + numpy.pi * 3.0/2.0
			r    = tr[:, 1:2]
			x	 = xy[:, 0:1]
			y	 = xy[:, 1:2]
			x[:] = r * numpy.cos(t)
			y[:] = r * numpy.sin(t)
			return xy

	class NorthPolarTransform(PolarAxes.PolarTransform):
		def transform(self, tr):
			xy	 = numpy.zeros(tr.shape, numpy.float_)
			t    = tr[:, 0:1] + numpy.pi * 3.0/2.0
			r    = tr[:, 1:2]
			x	 = xy[:, 0:1]
			y	 = xy[:, 1:2]
			x[:] = (1-r) * numpy.cos(t)
			y[:] = (1-r) * numpy.sin(t)
			return xy

		transform_non_affine = transform

	########################################################################

	def __init__(self, *args, **kwargs):
		self.lon_degs = 45
		self.lat_degs = 15
		PolarAxes.__init__(self, *args, **kwargs)

    ########################################################################

	def cla(self):
		"""
			Override to set up sensible defaults
		"""
		Axes.cla(self)
		self.set_longitude_grid(45)
		self.set_latitude_grid(15)
		self.xaxis.set_ticks_position('none')
		self.yaxis.set_ticks_position('none')

    ########################################################################

	def _set_lim_and_transforms(self):
		PolarAxes._set_lim_and_transforms(self)
		self.transProjection = self.NorthPolarTransform()
		self.axisProjection = self.NorthPolarAxisTransform()
		self.transData = (
			Affine2D().scale(numpy.pi/180.0, 1.0/90.0) + 
			self.transScale + 
			self.transProjection + 
			(self.transProjectionAffine + self.transAxes))
		self._xaxis_pretransform = (
			Affine2D().scale(1.0, 1.0) )
		self._xaxis_transform = (
			self._xaxis_pretransform + 
			self.axisProjection +
			self.PolarAffine(IdentityTransform(), Bbox.unit()) +
			self.transAxes)
		self._xaxis_text1_transform = (
			self._theta_label1_position +
			self._xaxis_transform)
		self._yaxis_transform = (
			Affine2D().scale(numpy.pi*2.0, 1.0) + 
			self.transScale + 
			self.axisProjection + 
			(self.transProjectionAffine + self.transAxes))
		self._yaxis_text1_transform = (
			self._r_label1_position +
			Affine2D().translate(0.0, -0.051) +
			self._yaxis_transform)

	########################################################################

	def set_longitude_grid(self, degrees=-1):
		if degrees == -1:
			degrees = self.lon_degs
		self.lon_degs = degrees
		number = (360.0 / degrees)
		ts = numpy.linspace(0, 360-degrees, number, True)
		th_labels = []
		for t in range(0, ts.shape[0]):
			th_labels.append(StringFormatEastWest(ts[t]))
		self.set_thetagrids(ts, th_labels)
		self.set_latitude_grid(-1)

	########################################################################

	def set_latitude_grid(self, degrees=-1):
		if degrees == -1:
			degrees = self.lat_degs
		self.lat_degs = degrees
		number = (90.0 / degrees) + 1
		rs = numpy.linspace(1e-10, 1.0, number, True)
		r_labels = [""]
		for r in range(1, rs.shape[0]-1):
			r_labels.append(StringFormatNorthSouth(90-rs[r]*90))
		self.set_rgrids(rs, r_labels, angle=1.0-self.lon_degs/360.0)
		# set the last transformation
		lat_mv = (rs[-1] - rs[-2])
		lat_sc = (rs[-1] - rs[1])
		self._xaxis_pretransform.clear().scale(1.0, lat_sc).translate(0.0, lat_mv)

    ########################################################################

	def set_xlim(self, *args, **kwargs):
		"""
		Override set_xlim and set_ylim to be between -pi, pi and 0, 1.0
		"""
		PolarAxes.set_xlim(self, -2*numpy.pi, 2*numpy.pi)
		PolarAxes.set_ylim(self, 0, 1)

    ########################################################################

	def set_ylim(self, *args, **kwargs):
		self.set_xlim()
