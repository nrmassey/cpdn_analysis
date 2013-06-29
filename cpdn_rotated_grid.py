###############################################################################
# File         : cpdn_rotated_grid.py
# Author       : Neil Massey
# Created      : 29/05/13
# Purpose      : class to represent a rotated grid
# Changes      : 
###############################################################################

from conv_rot_grid import *
import numpy

###############################################################################

class cpdn_rotated_grid:
	"""Class to represent a rotated grid, its pole and coordinates"""
	
	###########################################################################

	def __init__(self, rotated_pole_lat, rotated_pole_lon, lat_vals, lon_vals,
				 slice_with_global_coords = True):
		# initialise with all the available information
		self.__rotated_pole_lat = rotated_pole_lat
		self.__rotated_pole_lon = rotated_pole_lon
		self.__rotated_lats = lat_vals
		self.__rotated_lons = lon_vals
		# create the global latitude and longitude coordinates
		self.__global_lats = numpy.zeros([len(lat_vals), len(lon_vals)], 'f')
		self.__global_lons = numpy.zeros([len(lat_vals), len(lon_vals)], 'f')
		self.__slice_with_global_coords = slice_with_global_coords
		# calculate the global latitude and longitude coordinates
		for y in range(0, len(lat_vals)):
			for x in range(0, len(lon_vals)):
				v = rot2reg(lon_vals[x], lat_vals[y],
							rotated_pole_lon, rotated_pole_lat)
				self.__global_lons[y,x] = v[0]
				self.__global_lats[y,x] = v[1]

	###########################################################################

	def get_global_latitude_values(self):
		return self.__global_lats
	
	###########################################################################

	def get_global_longitude_values(self):
		return self.__global_lons
	
	###########################################################################

	def get_rotated_latitude_values(self):
		return self.__rotated_lats
	
	###########################################################################

	def get_rotated_longitude_values(self):
		return self.__rotated_lons
	
	###########################################################################

	def get_rotated_pole_latitude(self):
		return self.__rotated_pole_lat

	###########################################################################
	
	def get_rotated_pole_longitude(self):
		return self.__rotated_pole_lon
	
	###########################################################################

	def translate_to_grid(self, lat_val, lon_val):
		"""Translate a global (true) latitude and longitude onto this rotated
		   grid."""
		return reg2rot(lon_val, lat_val,
					   self.__rotated_pole_lon, self.__rotated_pole_lat)
	
	###########################################################################

	def set_slice_method(self, sm):
		"""Set the slicing method used in the subsetting of variables.
		   Really only of use for rotated grids.
		   sm = "rotated" : use rotated grid coordinates around the rotated pole
		   sm = "global" : use global grid coordinates 
		   					- i.e. the true latitude and longitude"""
		if sm == "rotated":
			self.__slice_with_global_coords = False
		elif sm == "global":
			self.__slice_with_global_coords = True
	
	###########################################################################

	def slice_with_global_coords(self):
		return self.__slice_with_global_coords

