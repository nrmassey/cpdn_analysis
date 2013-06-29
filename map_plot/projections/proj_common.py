############################################################################
# File         : proj_common.py
# Author       : Neil Massey
# Created      : 08/03/12
# Purpose      : Common functions for projections
#
############################################################################

from matplotlib.ticker import FixedLocator, Formatter
from matplotlib.axes import Axes
import numpy

########################################################################

def StringFormatEastWest(x):
	degrees = round(x)
	# \u00b0 : degree symbol
	if degrees > 180.0:
		degrees -= 360
	if degrees == 0.0 or degrees == 180.0:
		hemi = ""
	elif degrees == -180.0:
		hemi = ""
		degrees *= -1
	elif degrees < 0.0: 
		degrees *= -1
		hemi = "W"
	else:
		hemi = "E"
	return u"%d\u00b0" % degrees + hemi

########################################################################

def StringFormatNorthSouth(x):
	degrees = round(x)
	# \u00b0 : degree symbol
	if degrees == 0.0:
		hemi = ""
	elif degrees < 0.0: 
		degrees *= -1
		hemi = "S"
	else:
		hemi = "N"
	return u"%d\u00b0" % degrees + hemi

########################################################################

class EastWestDegreeFormatter(Formatter):
	"""Formatter to display coords as (-)180W to 180E"""
	def __init__(self, round_to=1.0):
		self._round_to = round_to

	def __call__(self, x, pos=None):
		return StringFormatEastWest(x)

########################################################################

class NorthSouthDegreeFormatter(Formatter):
	"""Formatter to display coords as (-)180W to 180E"""
	def __init__(self, round_to=1.0):
		self._round_to = round_to

	def __call__(self, x, pos=None):
		return StringFormatNorthSouth(x)
