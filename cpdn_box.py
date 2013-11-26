###############################################################################
# File         : cpdn_box.py
# Author       : Neil Massey
# Created      : 18/05/12
# Purpose      : class to represent a netCDF variable, along with its associated
#                dimensions and attributes.  Allows slicing based on contents
#                of the dimension variables
# Changes      : 
###############################################################################

# import the netCDF reading library
#from scipy.io.netcdf import *		# use this line to use the correct SciPy netCDF input / output
from netcdf_file import *			# my scipy is so out of date I've got a local copy of this
from collections import Iterable
import numpy
from daytime import *
from datetime import datetime
from copy import *
from types import *
from cpdn_rotated_grid import cpdn_rotated_grid
import sys

###############################################################################

def get_attributes(ncvar):
	# return the attributes of a variable as a dictionary
	attr = ncvar._attributes
	attr_dict = {}
	for a in attr:
		if a not in ["assignValue", "getValue", "typecode"]:
			attr_dict[a] = getattr(ncvar, a)
	return attr_dict

###############################################################################

class cpdn_boxdim:
	"""Class to represent a netCDF dimension, along with its axis, bounds, name,
		attributes and values"""

	###########################################################################

	def __init__(self, name="", vals=[], attrs={}, axis="", bounds=[],
				 start_date=None, time_units=None, n_days_per_year=None):
		# allow an empty constructor as well as a full one
		if name != "":
			self.__dim_name = name
		if vals != "":
			self.__dim_vals = numpy.array(vals)	# use numpy.array to copy values so as
												# to free up the reference to the ncfile
		self.__dim_attrs = attrs
		if axis != "":
			self.__dim_axis = axis
		elif name != "" or attrs != {}:		# have to have something to guess from
			self.__load_dimension_axis()
		# check for attributes vs values for time types
		if axis=="T" and \
		   (isinstance(n_days_per_year, NoneType) or \
 		    isinstance(start_date, NoneType) or \
            isinstance(time_units, NoneType)):
			self.__interpret_time_units()
		else:
			self.__n_days_per_year = n_days_per_year
			self.__start_date = start_date
			self.__time_units = time_units
		if bounds != []:
			self.__dim_bounds = bounds
		elif vals != []:					# have to have something to guess from
			self.__guess_dimension_bounds()

	###########################################################################

	def __load_dimension_axis(self):
		# get the axis of the dimension in the order they appear
		# axis are "T", "X", "Y", "Z"
		self.__dim_axis = ""
		#
		axis_std_names = { "X" : ["X", "longitude"],
						   "Y" : ["Y", "latitude"],
						   "Z" : ["Z", "height", "pressure", "level", "surface", "depth"],
						   "T" : ["T", "t", "time"] }
		#
		axis_dim_names = { "X" : ["X", "x", "longitude", "lon"],
						   "Y" : ["Y", "y", "latitude", "lat"],
						   "Z" : ["Z", "z", "height", "pressure", "level", "surface", "depth"],
						   "T" : ["T", "t", "time"] }
		#
		# try the "axis" attribute first			
		axis = None
		if "axis" in self.__dim_attrs.keys():
			axis = self.__dim_attrs["axis"]
		# else try to guess on the standard name
		elif axis == None and "standard_name" in self.__dim_attrs.keys():
			std_name = self.__dim_attrs["standard_name"]
			for ak in axis_std_names.keys():
				for an in axis_std_names[ak]:
					if an == std_name:
						axis = ak
							
		# finally guess on the dimension name
		if axis == None:
			for ak in axis_dim_names.keys():
				for an in axis_dim_names[ak]:
					if self.__dim_name[0:len(an)] == an:
						axis = ak
		self.__dim_axis = axis

		# add the axis attribute if it isn't present
		if not "axis" in self.__dim_attrs.keys():
			self.__dim_attrs["axis"] = axis

	###########################################################################

	def __guess_dimension_bounds(self):
		# guess the bounds
		l = self.__dim_vals.shape[0]		# length of dimension
		bounds_vals = numpy.zeros((l, 2), 'f')
		for e in range(1, l-1):
			# calculate offsets between this grid point and previous 
			# grid_point and next grid_point
			o1 = 0.5*(self.__dim_vals[e] - self.__dim_vals[e-1])
			o2 = 0.5*(self.__dim_vals[e+1] - self.__dim_vals[e])
			# create from offsets
			bounds_vals[e,0] = self.__dim_vals[e] - o1
			bounds_vals[e,1] = self.__dim_vals[e] + o2
		if l >= 2:
			# special case for first ...
			of =  0.5*(self.__dim_vals[1] - self.__dim_vals[0])
			bounds_vals[0,0] = self.__dim_vals[0] - of
			bounds_vals[0,1] = self.__dim_vals[0] + of
			# ... and last
			ol =  0.5*(self.__dim_vals[l-1] - self.__dim_vals[l-2])
			bounds_vals[l-1,0] = self.__dim_vals[l-1] - ol
			bounds_vals[l-1,1] = self.__dim_vals[l-1] + ol
		else:
			# if we're the time axis then we'll have a go at guessing
			# the bounds
			if self.__dim_axis == "T":
				for e in range(0, l):
					# convert to day time
					dv = self.__dim_vals[e] / self.__time_units
					dt = float_to_daytime(dv, self.__n_days_per_year)
					# monthly mean - note need to modify this for calendar != 360 day
					if dt.day in [14,15,16] and dt.hour == 0:
						bounds_vals[e,0] = (self.__dim_vals[e] - dt.day) + 1
						bounds_vals[e,1] = self.__dim_vals[e] + (30 - dt.day) + 1
					elif dt.hour == 12:		# daily mean
						bounds_vals[e,0] = self.__dim_vals[e] - 0.5
						bounds_vals[e,1] = self.__dim_vals[e] + 0.5
					else:					# give up !
						bounds_vals[e,0] = self.__dim_vals[e]
						bounds_vals[e,1] = self.__dim_vals[e]
			else:
				# no chance of guessing - bounds are the same at start and end
				for e in range(0, l):
					bounds_vals[e,0] = self.__dim_vals[e]
					bounds_vals[e,1] = self.__dim_vals[e]
		self.__dim_bounds = bounds_vals

	###########################################################################

	def __load_dimension_bounds(self, ncfile):
		# get the bounds for the dimensions.
		# if the dimension does not have the "bounds" or "bnds" attribute
		# then guess at the bounds
		#
		# get the name of the variable containing the dimension bounds
		if "bounds" in self.__dim_attrs.keys():
			bounds_var_name = self.__dim_attrs["bounds"]
		elif "bnds" in self.__dim_attrs.keys():
			bounds_var_name = self.__dim_attrs["bnds"]
		else:
			bounds_var_name = None
		# load the bounds from the bounds_dim_var
		if not isinstance(bounds_var_name, NoneType):
			# check the bounds name actually exists as a variable in the file
			if bounds_var_name in ncfile.variables.keys():
				bounds_var = ncfile.variables[bounds_var_name]
				self.__dim_bounds = bounds_var[:]
			else:
				bounds_var_name = None
				self.__guess_dimension_bounds()
		else:
			self.__guess_dimension_bounds()

	###########################################################################

	def __interpret_time_units(self):
		# dn = dimension number
		# get the unit string
		if "units" in self.__dim_attrs.keys():
			unit_string = self.__dim_attrs["units"]
		# get the number of calendar days in the year
		self.__n_days_per_year = 365.25
		if "calendar" in self.__dim_attrs.keys():
			if self.__dim_attrs["calendar"] == "360_day":
				self.__n_days_per_year = 360
			if self.__dim_attrs["calendar"] == "noleap":
				self.__n_days_per_year = 365

		# determine the number of days between each time stamp
		unit_types = {"seconds since" : 86400,
					  "minutes since" : 1440,
					  "hours since" : 24.0,
					  "days since"  : 1.0,
					  "years since" : 1.0/self.__n_days_per_year}
		for k in unit_types.keys():
			if k in unit_string:
				start_date_string = unit_string[len(k):].strip(" ")
				# two possible formats for daytime
				try:
					do = datestring_to_daytime(start_date_string)
				except:
					raise Exception("Unknown date format : " + start_date_string + \
	                                ", in attribute : units, in variable :" + \
					                self.__dim_name)
				# have to reinterpret the daytime as the year might be 360 days, not 365
				self.__start_date = daytime_to_float(do, self.__n_days_per_year)
				self.__time_units = unit_types[k]

	###########################################################################

	def __interpret_index_item(self, si):
		# interpret an item used to index the array
		idx_type = "None"
		if isinstance(si, slice):
			idx_type = "Slice"
		elif isinstance(si, basestring):
			idx_type = "String"
		elif isinstance(si, daytime):
			idx_type = "Daytime"
		elif isinstance(si, int):
			idx_type = "Integer"
		elif isinstance(si, float):
			idx_type = "Float"
		return idx_type

	###########################################################################

	def __get_idx_in_dim(self, val):
		if self.__interpret_index_item(self.__dim_bounds[0,1] ) != "String" and \
		   self.__dim_bounds[0,1] - self.__dim_bounds[0,0] > 0:
			idx = numpy.intersect1d(\
					numpy.where(val >= self.__dim_bounds[:,0])[0],\
					numpy.where(val <  self.__dim_bounds[:,1])[0])
		else:
			# need to do reverse search on negatively increasing dimensions
			idx = numpy.intersect1d(\
					numpy.where(val <= self.__dim_bounds[:,0])[0],\
					numpy.where(val >  self.__dim_bounds[:,1])[0])
		if idx.size == 0:
			return None
		else:
			return idx

	###########################################################################

	def __get_idx(self, si):
		# get a single index corresponding to the si
		# see __get_idxs below to determine the type
		idx = None
		#
		idx_type = self.__interpret_index_item(si)
		# check that a nested slice hasn't been passed in !
		if idx_type == "Slice":
			raise Exception("Failed indexing dimension " + self.__dim_name +\
                            " nested slice attempted with: " + str(si) + ".")
		#
		# check whether it's a date-formatted string
		elif idx_type == "String":
			# can only access date / time axis with a string
			if self.__dim_axis != "T":
				raise Exception("Failed indexing dimension " + self.__dim_name +\
								" with " + idx_type + ", as it is not a time dimension.")
			# other wise convert string to daytime - two possible formats
			try:
				do = datestring_to_daytime(si)
			except:
				raise Exception("Unknown date format : " + si)
			# now get the date as the offset to the start date
			ff = daytime_to_float(do, self.__n_days_per_year)
			idx_float_date = self.__time_units * \
							 (daytime_to_float(do, self.__n_days_per_year) - self.__start_date)
			# get the index
			idx = self.__get_idx_in_dim(idx_float_date)
		#
		# daytime next
		elif idx_type == "Daytime":
			if self.__dim_axis != "T":
				raise Exception("Failed indexing dimension " + self.__dim_name +\
								" with " + idx_type + ", as it is not a time dimension.")
			idx_float_date = self.__time_units * \
							 (daytime_to_float(si, self.__n_days_per_year) - self.__start_date)
			# get the index
			idx = self.__get_idx_in_dim(idx_float_date)
		#
		# if it's an integer then try to get an index, followed by a value
		elif idx_type == "Integer":
			# try to get the index first
			if si >= 0 and si <= self.__dim_vals.shape[0]:
				idx = si
			# otherwise get it in the axis
			else:
				idx = self.__get_idx_in_dim(si)
				if idx == None:
					# weird thing here to do with default slice of [:] giving slice(0,2**31-1,0)!
					if si == 2**31-1:
						idx = self.__dim_vals.shape[0]
					else:
						raise Exception("Failed indexing dimension " + self.__dim_name +\
    	                                " with " + str(si) + " as it is out of range of the dimension")
		#
		# if it's a floating point then just look in the values for the bounds
		elif idx_type == "Float":
			idx = self.__get_idx_in_dim(si)
			if idx == None:
				raise Exception("Failed indexing dimension " + self.__dim_name +\
                                " with " + str(si) + " as it is out of range of the dimension")
		return idx

	###########################################################################

	def get_idxs(self, si):
		# get an index corresponding to the si, which might be an index or a dimension value
		# have to determine: 
		# 1. if it's a slice (will have to determine what each slice component is)
		# 2. if it's a daytime object (will be a dimension value)
		# 3. if it's a string formatted date (will be a dimension value)
		# 4. if it's an integer (may be an index or a dimension value)
		# 5. if it's a floating point number (will be a dimension value)

		idx_type = self.__interpret_index_item(si)

		start = None
		end = None
		step = None
		# if it's a slice then call __get_idx on each component
		if idx_type == "Slice":
			start = self.__get_idx(si.start)
			end   = self.__get_idx(si.stop)
			if end != None:
				end += 1
			step  = self.__get_idx(si.step)
		else:
			start = self.__get_idx(si)
			if start != None:
				end   = start+1
		return start, end, step

	###########################################################################

	def load_name_axis_attributes(self, ncfile, ncvar, ncdim):
		# get the dimensions for the variable in the order in which they occur
		#
		# set the name
		self.__dim_name = ncdim
		# set the dimension values
		dim_var = ncfile.variables[ncdim]
		vals = dim_var[:]
		self.__dim_vals = numpy.array(vals) # use numpy.array to copy values so as
											# to free up the reference to the ncfile
		# set the dimension attributes
		self.__dim_attrs = get_attributes(dim_var)
		# load the dimension attribute
		self.__load_dimension_axis()

		# get the time startdate - if a time dimension
		# get a time dimension
		if self.__dim_axis == "T":
			self.__start_date = 0.0
			self.__time_units = 1
			# get the units - if it exists
			if "units" in self.__dim_attrs.keys():
				self.__interpret_time_units()

		self.__load_dimension_bounds(ncfile)

	###########################################################################

	def get_name(self):
		return self.__dim_name

	###########################################################################

	def get_axis(self):
		return self.__dim_axis

	###########################################################################

	def get_attributes(self):
		return self.__dim_attrs

	###########################################################################

	def get_values(self, t_mode="string"):
		"""Get the values for the dimension.
		   returns  : dimension values as a 1-D array"""
		# if this is a time dimension then convert to date times and return
		if self.__dim_axis == "T" and t_mode != "true_value": # only used for saving out!
			time_values = []
			for t in range(0, self.__dim_vals.shape[0]):			
				tv = self.__start_date + self.__dim_vals[t] * 1.0/self.__time_units
				dtv = float_to_daytime(tv, self.__n_days_per_year)
				if t_mode == "daytime":
					time_values.append(dtv)
				elif t_mode == "string":
					time_values.append(dtv.isoformat(" "))
				else:
					time_values.append(tv)
			if t_mode == "value":
				return numpy.array(time_values, 'f')
			else:
				return time_values
		else:
			return self.__dim_vals

	###########################################################################

	def set_values(self, values):
		"""Set the dimensions values to have the input values.
			Note for time dimension this is the "True" values."""
		self.__dim_vals = values

	###########################################################################

	def get_bounds(self, t_mode="string"):
		"""Get the bounds of the dimension.
		   returns  : bounds in array [n,2] - low bound is [m,0], high bound is [m,1]"""
		if self.__dim_axis == "T" and t_mode != "true_value":
			time_bounds = []
			for t in range(0, self.__dim_bounds.shape[0]):
				# lower bound
				tv1 = self.__start_date + self.__dim_bounds[t,0] * self.__time_units
				dtv1 = float_to_daytime(tv1, self.__n_days_per_year)
				# upper bound
				tv2 = self.__start_date + self.__dim_bounds[t,1] * self.__time_units
				dtv2 = float_to_daytime(tv2, self.__n_days_per_year)
				if t_mode == "daytime":
					time_bounds.append([dtv1, dtv2])
				elif t_mode == "string":
					time_bounds.append([dtv1.isoformat(" "), dtv2.isoformat(" ")])
				else:
					time_bounds.append([tv1, tv2])
			return numpy.array(time_bounds)
		else:
			return self.__dim_bounds

	###########################################################################

	def get_len(self):
		return len(self.__dim_vals)

	###########################################################################

	def get_length(self):
		return len(self.__dim_vals)

	###########################################################################

	def get_time_details(self):
		if self.__dim_axis != "T":
			raise "Trying to get time details for non-time dimension"
		return self.__start_date, self.__time_units, self.__n_days_per_year

###############################################################################

class cpdn_box_iterator:
	"""Class to allow iteration over a box whilst keeping some dimensions
       constant.  e.g. iterate over all the time and Z dimensions, returning
	   the X & Y as fields"""
	#
	###########################################################################
	#
	def __init__(self, box, axes=[]):
		"""Initialise iterator with an already created box
		   Inputs : box  - cpdn_box.
					axes - list of axes to iterate over"""
		self.__box = box
		self.__it_axes = axes
		# maximum iteration point
		self.__m_it_pt = 1
		for a in self.__box.get_dimension_axes():
			if a in self.__it_axes:
				self.__m_it_pt *= box.get_dimension(a).get_len()
		self.begin()

	###########################################################################

	def begin(self, get_target=False):
		"""Reset the iterator to the beginning."""
		# current index
		self.__c_idx = []
		self.__c_it_pt = 0
		bd = self.__box.get_dimension_axes()
		for a in bd:
			# use slice(:) for dimensions not in iterator
			if not (a in self.__it_axes):
				self.__c_idx.append(slice(None, None, None))
			else:
				self.__c_idx.append(0)
		if get_target:
			# convert __c_idx into an index into the target object
			tgt_idx = []
			for a in range(0, len(bd)):
				if bd[a] in self.__it_axes:
					tgt_idx.append(self.__c_idx[a])
			return tuple(self.__c_idx), tuple(tgt_idx)
		else:
			return tuple(self.__c_idx)

	###########################################################################

	def end(self):
		"""Check whether the iterator is at the end"""
		return self.__c_it_pt == self.__m_it_pt

	###########################################################################

	def next(self, get_target=False):
		"""Advance to the next item in the iteration sequence.
		   Returns : index into source array, index into target array."""
		#
		# determine where to add the index to
		bd = self.__box.get_dimension_axes()
		for a in range(len(bd)-1, -1, -1):
			if bd[a] in self.__it_axes:
				if self.__c_idx[a] < self.__box.get_dimension(bd[a]).get_len()-1:
					self.__c_idx[a] += 1
					break
				else:
					self.__c_idx[a] = 0
		# increase the current iteration point and return
		self.__c_it_pt += 1
		if get_target:
			# convert __c_idx into an index into the target object
			tgt_idx = []
			for a in range(0, len(bd)):
				if bd[a] in self.__it_axes:
					tgt_idx.append(self.__c_idx[a])
			return tuple(self.__c_idx), tuple(tgt_idx)
		else:
			return tuple(self.__c_idx)

###############################################################################

class cpdn_box:
	"""Class to represent a netCDF variable, along with its dimensions and attributes"""

	###########################################################################

	def __init__(self, ncfile=None, ncvar=None, name=None, dims=None, 
				 var_attrs={}, glob_attrs={}, off=0.0, sf=1.0, data=None,
				 rotated_grid=None):
		self.__ncfile = ncfile
		self.__ncvar  = ncvar
		self.__name   = name
		self.__dims   = dims
		self.__var_attrs = var_attrs
		self.__glob_attrs = glob_attrs
		self.__off    = off
		self.__sf     = sf
		self.__data   = data
		self.__rotated_grid = rotated_grid
		# reshape data to match dimensions
		if self.__data != None and self.__dims != None:
			dim_lens = []
			for d in self.__dims:
				dim_lens.append(d.get_len())
			self.__data = self.__data.reshape(dim_lens)

		# reset offset and scale factor if in attributes
		if not isinstance(self.__var_attrs, NoneType):
			self.__var_attrs["scale_factor"] = sf
			self.__var_attrs["add_offset"] = off

			# convert mv into float
			if "missing_value" in self.__var_attrs.keys():
				self.__var_attrs["missing_value"] = numpy.float32(self.__var_attrs["missing_value"])
			if "_FillValue" in self.__var_attrs.keys():
				self.__var_attrs["_FillValue"] = numpy.float32(self.__var_attrs["_FillValue"])

	###########################################################################

	def load(self, filename, varname):
		"""Init prepares the box using the field from the filename and varname"""
		# open the file
		try:
			self.__ncfile = netcdf_file(filename, "r")
		except:
			raise IOError(filename + " not found")

		# get the variable - raise an exception if not found in the file
		if not varname in self.__ncfile.variables.keys():
			raise IOError(varname + " not found in " + filename)
		self.__ncvar  = self.__ncfile.variables[varname]
		self.__name = varname
		# get data from __ncvar
		self.__data = None
		# get the attributes - do this before loading the dimensions so that we have the
		# possibility of having rotated grid coordinates
		self.__var_attrs = get_attributes(self.__ncvar)
		# get the global attributes
		self.__glob_attrs = get_attributes(self.__ncfile)
		# get the dimensions
		self.__load_dimensions()
		# set scaling factor and offset
		if "scale_factor" in self.__var_attrs.keys():
			self.__sf = self.__var_attrs["scale_factor"]
		else:	
			self.__sf = 1.0
		if "add_offset" in self.__var_attrs.keys():
			self.__off = self.__var_attrs["add_offset"]
		else:
			self.__off = 0.0		
		# convert mv into float
		if "missing_value" in self.__var_attrs.keys():
			self.__var_attrs["missing_value"] = numpy.float32(self.__var_attrs["missing_value"])
		if "_FillValue" in self.__var_attrs.keys():
			self.__var_attrs["_FillValue"] = numpy.float32(self.__var_attrs["_FillValue"])

	###########################################################################

	def __del__(self):
		if self.__ncfile != None:
			self.__ncfile.close()

	###########################################################################

	def __calc_rotated_grid_idxs__(self, idx_args):
		# Translate the coordinates given on a global grid to those on a rotated
		# grid so that they can just be fed into the calc_idx method of the dimensions
		rg = self.get_rotated_grid()
		# first determine if x and y dimensions exist and get the original index if they do
		x_idx = -1
		y_idx = -1
		for i in range(0, len(idx_args)):
			# get the x and y original indices - these might be a slice
			if self.__dims[i].get_axis() == "X":
				x_oidx = idx_args[i]
				x_idx = i
			elif self.__dims[i].get_axis() == "Y":
				y_oidx = idx_args[i]
				y_idx = i
		# not found
		if x_idx == -1 or y_idx == -1:
			return idx_args
		# now translate the x_oidx and y_oidx into the rotated grid
		if isinstance(x_oidx, slice) and isinstance(y_oidx, slice):
			# generate 2 slices each with 2 translated coordinates for the start and 
			# end and just the spacing for the skip values
			# First check whether we should try to interpret these dimensions as slices with
			# latitude / longitude coordinates (floats) or indexes into the array (integer)
			translate = isinstance(y_oidx.start, float) & isinstance(x_oidx.start, float) & \
						isinstance(y_oidx.stop, float) & isinstance(x_oidx.stop, float)
			if translate:
				s0 = rg.translate_to_grid(y_oidx.start, x_oidx.start)
				s1 = rg.translate_to_grid(y_oidx.stop,  x_oidx.stop)
				x_slice = slice(s0[0], s1[0], x_oidx.step)
				y_slice = slice(s0[1], s1[1], y_oidx.step)
			else:
				if isinstance(x_oidx.step, NoneType):
					x_slice = slice(int(x_oidx.start), int(x_oidx.stop), None)
				else:
					x_slice = slice(int(x_oidx.start), int(x_oidx.stop), int(x_oidx.step))
				if isinstance(y_oidx.step, NoneType):
					y_slice = slice(int(y_oidx.start), int(y_oidx.stop), None)
				else:
					y_slice = slice(int(y_oidx.start), int(y_oidx.stop), int(y_oidx.step))
			#
			idx_args[x_idx], idx_args[y_idx] = x_slice, y_slice
		elif isinstance(x_oidx, slice) and not isinstance(y_oidx, slice):
			# check whether we should translate
			translate = isinstance(y_oidx, float) & isinstance(x_oidx.start, float) & \
						isinstance(x_oidx.stop, float)
			if translate:
				s0 = rg.translate_to_grid(y_oidx, x_oidx.start)
				s1 = rg.translate_to_grid(y_oidx, x_oidx.stop)
				x_slice = slice(s0[0], s1[0], x_oidx.step)
				y_slice = slice(s0[1], s1[1], None)
			else:
				if isinstance(x_oidx.step, NoneType):
					x_slice = slice(int(x_oidx.start), int(x_oidx.stop), None)
				else:
					x_slice = slice(int(x_oidx.start), int(x_oidx.stop), int(x_oidx.step))
				y_slice = slice(y_oidx, y_oidx, None)
			idx_args[x_idx], idx_args[y_idx] = x_slice, y_slice
		elif isinstance(y_oidx, slice) and not isinstance(x_oidx, slice):
			# check whether we should translate
			translate = isinstance(x_oidx, float) & isinstance(y_oidx.start, float) & \
						isinstance(y_oidx.stop, float)
			if translate:
				s0 = rg.translate_to_grid(y_oidx.start, x_oidx)
				s1 = rg.translate_to_grid(y_oidx.stop, x_oidx)
				x_slice = slice(s0[0], s1[0], None)
				y_slice = slice(s0[1], s1[1], y_oidx.step)
			else:
				if isinstance(y_oidx.step, NoneType):
					y_slice = slice(int(y_oidx.step), int(y_oidx.stop), None)
				else:
					y_slice = slice(int(y_oidx.step), int(y_oidx.stop), int(y_oidx.step))
				x_slice = (x_oidx, x_oidx, None)
			idx_args[x_idx], idx_args[y_idx] = x_slice, y_slice
		else:
			translate = isinstance(x_oidx, float) & isinstance(y_oidx, float)
			if translate:
				idx_args[x_idx], idx_args[y_idx] = rg.translate_to_grid(y_oidx, x_oidx)
			else:
				idx_args[x_idx], idx_args[y_idx] = x_oidx, y_oidx
		return idx_args

	###########################################################################

	def __calc_idxs__(self, *args):
		"""Calculate and return indices from the arguments."""
		idx_list = []
		if not isinstance(args[0], Iterable):
			idx_args = list(args)
		else:
			idx_args = list(args[0])
		# if this is a rotated grid and slicing is to be performed on the global grid
		# then we have to convert the indices in idx_args into those to be used on the
		# rotated grid
		if self.has_rotated_grid() and self.get_rotated_grid().slice_with_global_coords():
			idx_args = self.__calc_rotated_grid_idxs__(idx_args)
		#
		for i in range(0, len(idx_args)):
			# check whether it's a single index or not
			start, stop, step = self.__dims[i].get_idxs(idx_args[i])
			idx_list.append(slice(start,stop,step))
		# append the rest of the dimensions as slices so that we can use
		# them to access the dimension values
		for i in range(len(idx_list), len(self.__dims)):
			idx_list.append(slice(0, self.__dims[i].get_len(), None))
		# convert to tuple
		idx_list = tuple(idx_list)
		return idx_list

	###########################################################################

	def __getitem__(self, *args):
		"""Return a box containing the subsection of data, the relevant dimension
		   coordinates and retaining the box's attributes"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to read in a netCDF file")
		# get the indices
		idx_list = self.__calc_idxs__(*args)
		# get the data
		# if this box is a result of subsetting then __data will exist
		if self.__data != None:
			new_data = self.__data.__getitem__(idx_list)
		# otherwise it has been loaded in so use self.__ncvar
		else:
			new_data = self.__ncvar.__getitem__(idx_list)
		# create the new dimensions for the box
		new_dims = []
		new_dim_lens = []
		for d in range(0, len(self.__dims)):
			# get the dimension
			c_dim = self.__dims[d]
			# sub-set the dimension values according to the above
			new_dim_vals = c_dim.get_values(t_mode="true_value").__getitem__(idx_list[d])
			new_dim_bounds = c_dim.get_bounds(t_mode="true_value").__getitem__(idx_list[d])
			# start date and time units for T dim
			if c_dim.get_axis() == "T":
				start_date, time_units, n_days_per_year = c_dim.get_time_details()
			else:
				start_date, time_units, n_days_per_year = None, None, None
			# create a new dimension with the same name, axis and attributes, but with
			# the subsetted dimension values and bounds
			new_boxdim = cpdn_boxdim(c_dim.get_name(), new_dim_vals, 
									 c_dim.get_attributes(), c_dim.get_axis(), 
									 new_dim_bounds, start_date, time_units,
									 n_days_per_year)
			new_dims.append(new_boxdim)

		# does the original box have a rotated grid?
		if self.has_rotated_grid():
			rg = self.get_rotated_grid()
			# find the X dimension and Y dimension
			for nd in new_dims:
				if nd.get_axis() == "Y":
					y_dim = nd
				if nd.get_axis() == "X":
					x_dim = nd
			new_rotated_grid = cpdn_rotated_grid(rg.get_rotated_pole_latitude(),
												 rg.get_rotated_pole_longitude(),
												 y_dim.get_values(),
												 x_dim.get_values())
		else:
			new_rotated_grid = None

		# create a new box with the same name and attributes as the current box
		# but subsetted dimensions - no netCDF file or variable is associated
		# with the new box
		new_box = cpdn_box(None, None, self.__name, new_dims, self.__var_attrs, 
						   self.__glob_attrs, self.__off, self.__sf, new_data,
						   new_rotated_grid)
		return new_box

	###########################################################################

	def __setitem__(self, *args):
		"""Return a box containing the subsection of data, the relavent dimension
			coordinates and retaining the box's attributes"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		# get the indices
		idx_list = self.__calc_idxs__(args)
		# if the data hasn't been loaded yet then we'll have to grab it
		if self.__data == None:
			self.__data = copy(self.__ncvar[:])
		self.__data.__setitem__(idx_list, args[1])

	###########################################################################

	def __merge_meta_data(self, left, right, exception_list=[]):
		# essentially merge two dictionaries together, not including the 
		# exception list
		new_meta_data = {}
		# new_meta_data just a copy of left without the exception_list
		for k in left.keys():
			if not k in exception_list:
				new_meta_data[k] = left[k]

		# now merge the right hand side metadata
		for k in right.keys():
			if not k in exception_list:
				# check whether the key is already in the metadata
				if k in new_meta_data.keys():
					# in meta_data but does it have the same value?
					m_data_equal = (right[k] == new_meta_data[k])
					if isinstance(m_data_equal, numpy.ndarray):
						m_data_equal = m_data_equal.any()
					if not m_data_equal:
						# create a new key
						suffix = 0
						while k+"_"+str(suffix) in new_meta_data.keys():
							suffix += 1
						# add to the data
						new_meta_data[k+"_"+str(suffix)] = right[k]
					# if it does have the same value then it's already in the 
					# dictionary so doesn't have to be added again
				else:
					# not in new meta data so just add
						new_meta_data[k] = right[k]
		return new_meta_data

	###########################################################################

	def __parse_meta_data(self, other, std_name_format, cell_method):
		# get attributes together
		self_attrs = self.get_attributes()
		other_attrs = other.get_attributes()

		# merge the meta-data apart from the stuff we're going to manipulate
		new_attrs = self.__merge_meta_data(self_attrs, other_attrs, 
									["cell_methods", "units", "standard_name", "history"])

		# construct the standard name
		if "standard_name" in self_attrs.keys() and \
		   "standard_name" in other_attrs.keys():
			self_std_name = self_attrs["standard_name"]
			other_std_name = other_attrs["standard_name"]
		else:
			self_std_name = self.get_name()
			other_std_name = other.get_name()
		new_attrs["standard_name"] = std_name_format % (self_std_name, other_std_name)

		# construct the units - nothing clever like unit conversion
		if "units" in self_attrs.keys() and "units" in other_attrs.keys():
			new_attrs["units"] = std_name_format % (self_attrs["units"], other_attrs["units"])

		# do the cell methods
		method_string = ""
		for d in self.get_dimension_names():
			method_string += d + ": "
		method_string += cell_method
		if "cell_methods" in self_attrs.keys():
			method_string += "; " + self_attrs["cell_methods"]
		if "cell_methods" in other_attrs.keys():
			method_string += "; " + other_attrs["cell_methods"]
		new_attrs["cell_methods"] = method_string

		# do the history
		history_string = ""
		if "history" in self_attrs.keys():
			history_string += self_attrs["history"]
		if "history" in other_attrs.keys():
			if "history" in self_attrs.keys():
				if other_attrs["history"] != self_attrs["history"]:
					history_string += "; " + other_attrs["history"]
			else:
				history_string += "; " + other_attrs["history"]
		# add new history
		if history_string != "":
			history_string += "; "
		history_string += datetime.now().isoformat() + " altered by CPDN: " + \
						  std_name_format % (self.get_name(), other.get_name())
		new_attrs["history"] = history_string
		return new_attrs

	###########################################################################

	def __add__(self, other):
		"""Add two boxes together and return a new box containing the sum of
			the two boxes."""
		# get the data for this box
		this_data = self.get_values()
		# get the missing value indices for this box data
		mv_t = self.get_missing_value()
		mv_t_idx = numpy.where(this_data == mv_t)
		# add integer or float to data
		if isinstance(other, int) or isinstance(other, float):
			# add the scalar to the data
			new_box_data = this_data + other
			# add the missing value back in
			new_box_data[mv_t_idx] = mv_t
			# create new box, everything is the same
			new_box = cpdn_box(dims=self.get_dimensions(), 
							   var_attrs=self.get_attributes(),
							   name=self.get_name(),
							   glob_attrs=self.get_global_attributes(),
							   off=0.0, sf=1.0, data=new_box_data)
		# add two boxes together
		elif isinstance(other, cpdn_box):
			# first check that the data is compatible
			if self.get_dimension_lengths() != other.get_dimension_lengths():
				raise Exception("Cannot add boxes : dimension lengths differ")

			# get the other box data and missing value indices
			other_data = other.get_values()
			mv_o = other.get_missing_value()
			mv_o_idx = numpy.where(other_data == mv_o)

			# add the data together
			new_box_data = this_data + other_data
			# restore missing values - to this boxes missing values
			new_box_data[mv_t_idx] = mv_t
			new_box_data[mv_o_idx] = mv_t

			# construct the new name and meta data
			new_box_name = self.get_name() + "_add_" + other.get_name()
			new_box_attrs = self.__parse_meta_data(other, "sum_of_%s_and_%s", "sum")

			# new global meta data just a merge between the two
			new_glob_attrs = self.__merge_meta_data(self.get_global_attributes(), 
													other.get_global_attributes())

			# create a new box with the new data, name and metadata
			new_box = cpdn_box(dims=self.get_dimensions(), var_attrs=new_box_attrs,
							   name=new_box_name, glob_attrs=new_glob_attrs,
							   off=0.0, sf=1.0, data=new_box_data)
		else:
			raise Exception("Cannot add boxes : type error for right hand side")
		return new_box

	###########################################################################

	def __sub__(self, other):
		"""Subtract two boxes from each other and return a new box containing 
			the difference of the two boxes."""
		# get the data for this box
		this_data = self.get_values()
		# get the missing value indices for this box data
		mv_t = self.get_missing_value()
		mv_t_idx = numpy.where(this_data == mv_t)

		# sub integer or float from data
		if isinstance(other, int) or isinstance(other, float):
			# add the scalar to the data
			new_box_data = this_data - other
			# put the missing value back in
			new_box_data[mv_t_idx] = mv_t
			# create new box, everything is the same
			new_box = cpdn_box(dims=self.get_dimensions(), 
							   var_attrs=self.get_attributes(),
							   name=self.get_name(),
							   glob_attrs=self.get_global_attributes(),
							   off=0.0, sf=1.0, data=new_box_data)
		# sub two boxes together
		elif isinstance(other, cpdn_box):
			# first check that the data is compatible
			if self.get_dimension_lengths() != other.get_dimension_lengths():
				raise Exception("Cannot add boxes : dimension lengths differ")

			# get the other box data and missing value indices
			other_data = other.get_values()
			mv_o = other.get_missing_value()
			mv_o_idx = numpy.where(other_data == mv_o)

			# sub the data
			new_box_data = this_data - other_data

			# restore missing values - to this box's missing values
			new_box_data[mv_t_idx] = mv_t
			new_box_data[mv_o_idx] = mv_t

			# construct the new name and meta data
			new_box_name = self.get_name() + "_sub_" + other.get_name()
			new_box_attrs = self.__parse_meta_data(other, "difference_of_%s_and_%s", "difference")

			# new global meta data just a merge between the two
			new_glob_attrs = self.__merge_meta_data(self.get_global_attributes(), 
												other.get_global_attributes())

			# create a new box with the new data, name and metadata
			new_box = cpdn_box(dims=self.get_dimensions(), var_attrs=new_box_attrs,
							   name=new_box_name, glob_attrs=new_glob_attrs,
							   off=0.0, sf=1.0, data=new_box_data)
		else:
			raise Exception("Cannot subtract boxes : type error for right hand side")
		return new_box

	###########################################################################

	def __mul__(self, other):
		"""Multiply two boxes together and return a new box containing the 
			product of the two boxes."""
		# multiply the data

		# get the data for this box
		this_data = self.get_values()
		# get the missing value indices for this box data
		mv_t = self.get_missing_value()
		mv_t_idx = numpy.where(this_data == mv_t)

		# multiply data by integer or float to data
		if isinstance(other, int) or isinstance(other, float):
			# add the scalar to the data
			new_box_data = this_data * other
			# put the missing value back in
			new_box_data[mv_t_idx] = mv_t
			# create new box, everything is the same
			new_box = cpdn_box(dims=self.get_dimensions(), 
							   var_attrs=self.get_attributes(),
							   name=self.get_name(),
							   glob_attrs=self.get_global_attributes(),
							   off=0.0, sf=1.0, data=new_box_data)
		# multply two boxes together
		elif isinstance(other, cpdn_box):
			# first check that the data is compatible
			if self.get_dimension_lengths() != other.get_dimension_lengths():
				raise Exception("Cannot multiply boxes : dimension lengths differ")

			# get the other box data and missing value indices
			other_data = other.get_values()
			mv_o = other.get_missing_value()
			mv_o_idx = numpy.where(other_data == mv_o)

			new_box_data = this_data * other_data
			# restore missing values - to this box's missing values
			new_box_data[mv_t_idx] = mv_t
			new_box_data[mv_o_idx] = mv_t

			# construct the new name and meta data
			new_box_name = self.get_name() + "_mul_" + other.get_name()
			new_box_attrs = self.__parse_meta_data(other, "product_of_%s_and_%s", "product")

			# new global meta data just a merge between the two
			new_glob_attrs = self.__merge_meta_data(self.get_global_attributes(), 
													other.get_global_attributes())

			# create a new box with the new data, name and metadata
			new_box = cpdn_box(dims=self.get_dimensions(), var_attrs=new_box_attrs,
							   name=new_box_name, glob_attrs=new_glob_attrs,
							   off=0.0, sf=1.0, data=new_box_data)
		else:
			raise Exception("Cannot multiply boxes : type error for right hand side")

		return new_box

	###########################################################################

	def __div__(self, other):
		"""Divide a box by another box and return a new box containing the 	
			result of the division."""
		# get the data for this box
		this_data = self.get_values()
		# get the missing value indices for this box data
		mv_t = self.get_missing_value()
		mv_t_idx = numpy.where(this_data == mv_t)

		# divide data by integer or float to data
		if isinstance(other, int) or isinstance(other, float):
			# add the scalar to the data
			new_box_data = this_data / other
			# put the missing value back in
			new_box_data[mv_t_idx] = mv_t
			# create new box, everything is the same
			new_box = cpdn_box(dims=self.get_dimensions(), 
							   var_attrs=self.get_attributes(),
							   name=self.get_name(),
							   glob_attrs=self.get_global_attributes(),
							   off=0.0, sf=1.0, data=new_box_data)
		# multply two boxes together
		elif isinstance(other, cpdn_box):
			# first check that the data is compatible
			if self.get_dimension_lengths() != other.get_dimension_lengths():
				raise Exception("Cannot add boxes : dimension lengths differ")

			# get the other box data and missing value indices
			other_data = other.get_values()
			mv_o = other.get_missing_value()
			mv_o_idx = numpy.where(other_data == mv_o)

			# divide the data
			new_box_data = this_data / other_data
			# restore missing values - to this box's missing values
			new_box_data[mv_t_idx] = mv_t
			new_box_data[mv_o_idx] = mv_t

			# construct the new name and meta data
			new_box_name = self.get_name() + "_div_" + other.get_name()
			new_box_attrs = self.__parse_meta_data(other, "ratio_of_%s_to_%s", "ratio")

			# new global meta data just a merge between the two
			new_glob_attrs = self.__merge_meta_data(self.get_global_attributes(), 
													other.get_global_attributes())

			# create a new box with the new data, name and metadata
			new_box = cpdn_box(dims=self.get_dimensions(), var_attrs=new_box_attrs,
							   name=new_box_name, glob_attrs=new_glob_attrs,
							   off=0.0, sf=1.0, data=new_box_data)
		else:
			raise Exception("Cannot multiply boxes : type error for right hand side")

		return new_box

	###########################################################################

	def __load_dimensions(self):
		# get the dimensions details - bounds etc.
		self.__dims = []
		for dim in self.__ncvar.dimensions:
			cpdn_dim = cpdn_boxdim()
			cpdn_dim.load_name_axis_attributes(self.__ncfile, self.__ncvar, dim)
			self.__dims.append(cpdn_dim)
			
		# rotated grid - does this variable have a rotated grip mapping
		if "grid_mapping" in self.__var_attrs.keys():
			grid_map_name = self.__var_attrs["grid_mapping"]
			# get the variable containing the grid mapping detail
			grid_map_var = self.__ncfile.variables[grid_map_name]
			# get the type of grid mapping - this could be extended to other grid mappings
			grid_map_type = grid_map_var._attributes["grid_mapping_name"]
			if "rotated_latitude_longitude" in grid_map_type:
				# now get the location of the rotated pole
				rotated_pole_lon = grid_map_var._attributes["grid_north_pole_longitude"]
				rotated_pole_lat = grid_map_var._attributes["grid_north_pole_latitude"]
				self.__rotated_grid = cpdn_rotated_grid(rotated_pole_lat, rotated_pole_lon,
														self.get_dimension("Y").get_values(),
														self.get_dimension("X").get_values())
			
	###########################################################################

	def get_ndims(self):
		"""Return the number of dimensions"""
		return len(self.__dims)

	###########################################################################

	def get_dimensions(self):
		"""Return all dimensions in a list"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		dims = []
		for d in range(0, len(self.__dims)):
			dims.append(self.__dims[d])
		return dims

	###########################################################################

	def get_dimension_names(self):
		"""Return all the dimension names"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		names = []
		for d in range(0, len(self.__dims)):
			names.append(self.__dims[d].get_name())
		return names

	###########################################################################

	def get_dimension_axes(self):
		"""Return all the dimension axes"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		axes = []
		for d in range(0, len(self.__dims)):
			axes.append(self.__dims[d].get_axis())
		return axes

	###########################################################################

	def get_dimension_lengths(self):
		"""Return the lengths of all dimensions"""
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		lengths = []
		for d in range(0, len(self.__dims)):
			lengths.append(self.__dims[d].get_len())
		return lengths

	###########################################################################

	def get_dimension(self, si):
		"""Get the dimension either by the name, or the axis, or the index.
		   Returns a dimension object"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")
		if isinstance(si, int):
			if si < 0 or si >= len(self.__dims):
				raise Exception("Index out of range : " + str(si))
			return self.__dims[si]
		else:
			axes = []
			names = []
			for d in range(0, len(self.__dims)):
				axes.append(self.__dims[d].get_axis())
				names.append(self.__dims[d].get_name())
			try:
				idx = axes.index(si)
			except:
				try:
					idx = names.index(si)
				except:
					raise Exception("Axis name or dimension name not found : " + str(si))
			return self.__dims[idx]

	###########################################################################
	
	def has_rotated_grid(self):
		return not isinstance(self.__rotated_grid, NoneType)
		
	###########################################################################
	
	def get_rotated_grid(self):
		if not self.has_rotated_grid():
			raise Exception("get_rotated_grid called on a box without a rotated grid")
		return self.__rotated_grid
	
	###########################################################################

	def get_values(self, mode="scale"):
		"""Return all the data for the box.  Use [] to subset the data first,
		   which will return another box for you to call get_data on.
		   mode=="scale" : apply the scale factor and offset
		   mode=="true_value" : do not apply the scale factor - just return as is"""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")

		# if this is a box that has previously been subset then data will != None
		# and so we should use that
		mv = self.get_missing_value()
		if self.__data != None:
			v = self.__data[:] 
			# get missing value
			mv_i = numpy.where(v == mv)
			# do the scaling
			if mode == "scale":
				v = v * self.__sf + self.__off
			else:
				v = v * 1.0			# copy array
			# restore missing value
			v[mv_i] = mv
		else:
		# otherwise it has been loaded in and we should use __ncvar
			v = self.__ncvar[:]
			if not isinstance(mv, NoneType):
				mv_i = numpy.where(v == mv)
			if mode == "scale":
				v = v * self.__sf + self.__off
			else:
				v = v * 1.0			# copy array
			if not isinstance(mv, NoneType):
				v[mv_i] = mv
		return v

	###########################################################################

	def set_values(self, data):
		"""Simple in place replacement of values."""
		if list(data.shape) != self.get_dimension_lengths():
			raise Exception("Data has different dimension lengths to the box")
		self.__data = data
		self.__sf = 1.0
		self.__off = 0.0

	###########################################################################

	def get_missing_value(self):
		"""Get the missing value from the attributes.  Return None if not
			present."""
		mv = numpy.inf
		# missing value has precedent over _FillValue - also scale and offset
		if "missing_value" in self.__var_attrs.keys():	
			mv = self.__var_attrs["missing_value"]
		elif "_FillValue" in self.__var_attrs.keys():
			mv = self.__var_attrs["_FillValue"]
		return mv

	###########################################################################

	def get_sf(self):
		return self.__sf

	###########################################################################

	def get_off(self):
		return self.__off

	###########################################################################

	def get_attributes(self):
		"""Return the variable attributes for the box."""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")

		return self.__var_attrs

	###########################################################################

	def set_attributes(self, attrs):
		self.__var_attrs = attrs

	###########################################################################

	def get_global_attributes(self):
		"""Return the variable attributes for the box."""
		# check the box is not empty
		if self.__data == None and self.__ncvar == None:
			raise Exception("Box is empty - use .load to load in a netCDF file")

		return self.__glob_attrs

	###########################################################################

	def get_name(self):
		"""Return the name of the variable / field associated with the box"""
		return self.__name

	###########################################################################

	def set_name(self, name):
		"""Set the name of the variable"""
		self.__name = name

	###########################################################################

	def _write_to_netcdf_file(self, out_file):
		"""Function to write the box out to an already open netcdf file."""
		dimension_names = []
		create_dim = False
		for d in self.__dims:
			# create the dimensions - first check whether the dimension already exists
			d_vals = d.get_values("true_value")
			if not d.get_name() in out_file.dimensions.keys():
				dim_name = d.get_name()
				create_dim = True
			else:
				# the dimension already exists, check whether the values are the same
				# if they are then assume the dimension is the same
				existing_dim_var = out_file.variables[d.get_name()]
				if numpy.alltrue(d_vals == existing_dim_var[:]):
					dim_name = d.get_name()			# they are the same
				else:
					# create a new dimension and have a new dimension name
					dn = 0
					dim_name = d.get_name()
					while dim_name in out_file.dimensions.keys():
						dim_name = d.get_name() + "_" + str(dn)
						dn += 1
					create_dim = True

			if create_dim:
				out_file.createDimension(dim_name, d.get_len())

				# create the variable for the dimension and set its values
				dim_var = out_file.createVariable(dim_name, d_vals.dtype, \
												  (dim_name,))
				dim_var[:] = d_vals
				# add the dimension attributes
				dim_attrs = d.get_attributes()
				for k in dim_attrs.keys():
					dim_var.__setattr__(k, dim_attrs[k])

				# if the bounds exist in the attributes then create them as variables
				if "bounds" in dim_attrs.keys():
					# if the bounds dimension has not been created yet then create it
					if not "bnds" in out_file.dimensions.keys():
						out_file.createDimension("bnds", 2)
					# get the values
					bnd_vals = d.get_bounds("true_value")
					# get the name of the bounds
					bnd_name = dim_attrs["bounds"]
					# check that the bounds name doesn't already exist as a variable
					bn = 0
					new_bnd_name = bnd_name
					while new_bnd_name in out_file.variables.keys():
						new_bnd_name = bnd_name + "_" + str(bn)
						bn += 1

					# create the variable - dimensions are the same as the dimension 
					# variable (lat, lon, time, etc.) and "bnds" - which is always 2
					bnd_var = out_file.createVariable(new_bnd_name, bnd_vals.dtype, \
            	                                      (dim_name, "bnds",))
					bnd_var[:] = bnd_vals
			# add to the dimension names for this variable
			dimension_names.append(dim_name)

		# written the dimensions - now write the variable
		# get the values
		if self.__ncvar != None:
			vals = self.__ncvar[:]
		else:
			vals = self.__data[:]
		# create the variable - make sure the name is unique
		vn = 0
		new_var_name = self.__name
		while new_var_name in out_file.variables.keys():
			new_var_name = self.__name + "_" + str(vn)
			vn += 1

		# add the global attributes
		glob_attrs = self.get_global_attributes()
		if glob_attrs != None:
			for k in glob_attrs.keys():
				out_file.__setattr__(k, glob_attrs[k])

		out_var = out_file.createVariable(new_var_name, vals.dtype, \
										  tuple(dimension_names))
		# add the variable attributes
		var_attrs = self.get_attributes()
		if var_attrs != None:
			for k in var_attrs.keys():
				out_var.__setattr__(k, var_attrs[k])
				
		# if this is a rotated grid then write out the grid mapping
		if self.has_rotated_grid():
			# new variable - name is variable attributes
			grid_map_var = out_file.createVariable(var_attrs["grid_mapping"], 'c', ())
			# say what type of projection this is
			grid_map_var.__setattr__("grid_mapping_name", "rotated_latitude_longitude")
			# add the attributes to the variable saying where the pole is located
			grid_map_var.__setattr__("grid_north_pole_latitude",
									 self.get_rotated_grid().get_rotated_pole_latitude())
			grid_map_var.__setattr__("grid_north_pole_longitude",
									 self.get_rotated_grid().get_rotated_pole_longitude())
									
		
		# write the values - do this last for efficiency
		out_var[:] = vals

	###########################################################################

	def save(self, filename):
		"""Save the box out to a netCDF file.
			Input : filename - name of file to write to."""
		# check for .nc extension
		if filename[-3:] != ".nc":
			filename += ".nc"
		# open file for writing
		out_file = netcdf_file(filename, "w")
		self._write_to_netcdf_file(out_file)
		out_file.close()

###############################################################################

def save_multi_box(filename, list_of_boxes):
	"""Save all the boxes in list_of_boxes into the same file."""
	# check for .nc extension
	if filename[-3:] != ".nc":
		filename += ".nc"
	# open file for writing
	out_file = netcdf_file(filename, "w")
	# write each box in turn
	for b in list_of_boxes:
		b._write_to_netcdf_file(out_file)
		out_file.flush()
	out_file.close()

###############################################################################

def amend_attributes(attrs, method_str=None, history_str=None):
	# add the cell_methods and history to the attributes
	tgt_attrs = copy(attrs)
	# amend the attributes to add the mean to the cell methods and to update the history
	if method_str != None:
		if "cell_methods" in tgt_attrs.keys():
			tgt_attrs["cell_methods"] += "; " + method_str
		else:
			tgt_attrs["cell_methods"] = method_str
	#
	if history_str != None:
		if "history" in tgt_attrs.keys():
			tgt_attrs["history"] += "; " + history_str
		else:
			tgt_attrs["history"] = history_str

	if "source" in tgt_attrs.keys():
		tgt_attrs["source"] = "climateprediction.net"

	return tgt_attrs
