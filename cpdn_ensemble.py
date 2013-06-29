###############################################################################
# File         : cpdn_ensemble.py
# Author       : Neil Massey
# Created      : 01/06/12
# Purpose      : routines to generate an ensemble of files from directory
#				 listings and return subsets of the ensemble based on the
#				 date or metadata
# Changes      : 
###############################################################################

from cpdn_box import *			# use the box for I/O and to get the time dim and metadata
from glob import glob
from types import *
from copy import *
import gzip
import os
import cPickle as pickle

###############################################################################

class cpdn_ensemble_member:
	"""Single member of the ensemble"""

	###########################################################################

	def __init__(self, filename, varname, t_dim=None, var_attrs=None, glob_attrs=None):
		self.__filename = filename
		self.__varname = varname
		# allowing t_dim, var_attrs and glob_attrs to be None allows overloading 
		# and creating from pickle
		if isinstance(t_dim, NoneType) or isinstance(var_attrs, NoneType) or \
		   isinstance(glob_attrs, NoneType):
			temp_box = cpdn_box()			
			temp_box.load(self.__filename, self.__varname)
			# quality check
			temp_vals = temp_box.get_values()
			self.__valid = True
			self.__valid &= numpy.isfinite(temp_vals).any()
			self.__valid &= (numpy.max(temp_vals) < 2e20)
			self.__valid &= (numpy.min(temp_vals) > -2e20)
		
		if isinstance(t_dim, NoneType):
			self.__t_dim = temp_box.get_dimension("T")			
		else:
			self.__t_dim = t_dim

		if isinstance(var_attrs, NoneType):
			self.__var_attrs = temp_box.get_attributes()
		else:
			self.__var_attrs = var_attrs

		if isinstance(glob_attrs, NoneType):
			self.__glob_attrs = temp_box.get_global_attributes()
		else:
			self.__glob_attrs = glob_attrs

	###########################################################################

	def save(self, fh):
		"""Save in pickle format to already opened filehandle."""
		pickle.dump(self.__filename, fh)
		pickle.dump(self.__varname, fh)
		pickle.dump(self.__t_dim, fh)
		pickle.dump(self.__var_attrs, fh)
		pickle.dump(self.__glob_attrs, fh)

	###########################################################################

	def is_valid(self):
		return self.__valid

	###########################################################################

	def get_t_dim(self):
		"""Get the time dim of the ensemble member"""
		return self.__t_dim

	###########################################################################

	def get_box(self):
		"""Return the box associated with the ensemble member"""
		# have to load the box first
		temp_box = cpdn_box()
		temp_box.load(self.__filename, self.__varname)
		return temp_box

	###########################################################################

	def get_filename(self):
		"""Return the filename associated with the ensemble member"""
		return self.__filename

	###########################################################################

	def get_varname(self):
		"""Return the variable name associated with the ensemble member"""
		return self.__varname
	
	###########################################################################

	def get_attributes(self):
		"""Return the variable attributes for the ensemble member."""
		return self.__attrs

	###########################################################################

	def get_global_attributes(self):
		"""Return the global attributes for the ensemble member."""
		return self.__glob_attrs

	###########################################################################

	def __getitem__(self, idx):
		"""Return the box's getitem as the ensemble member getitem"""
		# have to load the box to get the item
		temp_box = cpdn_box()
		temp_box.load(self.__filename, self.__varname)
		return temp_box.__getitem__(idx)

	###########################################################################

	def get_values(self):
		"""Return the box's values."""
		temp_box = cpdn_box()
		temp_box.load(self.__filename, self.__varname)
		return temp_box.get_values()

###############################################################################

class cpdn_ensemble:
	"""Class to store an ensemble of climate data files"""

	###########################################################################

	def __init__(self, varname=None, ensemble_list=[]):
		self.__varname = varname
		self.__ensemble_list = copy(ensemble_list)

	###########################################################################

	def add_file(self, filepath, limit=-1, recursive_level=0):
		"""Add files to the ensemble.  Supports * wildcards (e.g. *.pa*).
			Files will contain varname, as specified in the constructor.
			filepath = path to files, including wildcard
			limit = maximum number of members in ensemble"""
		# first check that a name has been given
		if self.__varname == None:
			raise Exception("Ensemble has no varname.  Either create when invoking constructor or load from previously created ensemble.")
		# get a list of files that match the wildcard and add the directories to it
		last_slash = filepath.rfind("/")
		path_only = filepath[:last_slash]
		wildcard = filepath[last_slash:]
		f_list = glob(filepath)
		# add directories to the f_list
		for d in os.listdir(path_only):
			if os.path.isdir(path_only + "/" + d):
				f_list.append(path_only + "/" + d)
		# loop through the files and create a box for each one
		for f in f_list:
			if limit != -1 and len(self.__ensemble_list) >= limit:
				break
			# if recursive and f is a directory then recurse into the directory
			# remove wildcard first
			if os.path.isdir(f):
				if recursive_level > 0:
					self.add_file(f + wildcard, limit, recursive_level-1)
			else:
				# if self.__varname does not exist in the file then an exception
				# will be raised and the file will not be added
				try:
					ensemble_member = cpdn_ensemble_member(f, self.__varname)
					if ensemble_member.is_valid():
						self.__ensemble_list.append(ensemble_member)
					else:
						print f + " : Invalid"
				except:
					pass

	###########################################################################

	def __date_in_slice_to_days(self, date_in):
		"""Convert the date_in slice to a list of days since 0000 in 365.25
		   day format (for consistency)"""
		# if it's not a slice then create one
		if isinstance(date_in, slice):
			date_in_idxs = [date_in.start, date_in.stop, date_in.step]
		else:
			date_in_idxs = [date_in, None, None]

		# loop through the date_slice and create a slice with floating point
		# representation of date since 0000 in 365.25 day format (for consistency)
		date_vals = []
		for	i in date_in_idxs:
			if isinstance(i, NoneType):
				date_vals.append(None)
			elif isinstance(i, daytime):
				date_vals.append(daytime_to_float(i))
			elif isinstance(i, basestring):
				day_date = datestring_to_daytime(i)
				date_vals.append(daytime_to_float(day_date))
			elif isinstance(i, float) or isinstance(i, numpy.float32)\
			  or isinstance(i, numpy.float64):
				date_vals.append(date_in)
		# amend the values if end is None - make end = start + 23/24
		if isinstance(date_vals[1], NoneType) and not\
			isinstance(date_vals[0], NoneType):
			date_vals[1] = date_vals[0] + 23.0/24
		return date_vals

	###########################################################################

	def subset_by_date(self, date_in):
		"""Get members that have valid data on the supplied date.  
			Date can either be a daytime object, a string or floating point
			representation or a slice thereof."""
		# convert the slice to days since 0000
		date_val = self.__date_in_slice_to_days(date_in)
		sd, tu, n_days_py = self.__ensemble_list[0].get_t_dim().get_time_details()

		# members to return
		return_members = []
		# loop through each ensemble member
		for e in self.__ensemble_list:
			try:
				# get date bounds in daytime format
				t_vals = e.get_t_dim().get_bounds("daytime")
				# convert to floating point of days since 0000
				start_t = daytime_to_float(t_vals[0][0], n_days_py)		# +		date_val[0] = |
				end_t = daytime_to_float(t_vals[-1][1], n_days_py)   	# +		date_val[1] = |
				# six possible cases for date inclusion in the range
				# first case : |+ +|        # successful
				# 2nd case   : +| +|        # successful
				# 3rd case   : |+ |+        # successful
				# 4th case   : +| |+        # successful
				# 5th case   : || ++        # not successful
				# 6th case   : ++ ||        # not successful
				if (date_val[0] <= start_t and date_val[1] >= end_t) or\
				   (date_val[0] > start_t and date_val[0] <= end_t and date_val[1] > end_t) or\
				   (date_val[0] <= start_t and date_val[1] <= end_t and date_val[1] > start_t) or\
				   (date_val[0] > start_t and start_t <= end_t and date_val[1] <= end_t and date_val[1] >= start_t):
					return_members.append(e)
			except:
				# error raised if date out of range - ignore
				pass

		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def subset_by_year(self, year_n):
		"""Get the members of the ensemble that occur within a certain year"""
		# loop through ensemble
		return_members = []
		for e in self.__ensemble_list:
			# get date bounds in daytime format
			t_vals = e.get_t_dim().get_values("daytime")
			# get start and end year
			start_y = t_vals[0].year
			end_y = t_vals[-1].year
			# yeat_n lies between start_date and end_date of ensemble member
			if year_n >= start_y and year_n <= end_y:
				return_members.append(e)
		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def subset_by_month(self, month_in):
		"""Get the members of the ensemble that occur within a certain month"""
		month_map = {"Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
					 "Jul": 7, "Aug": 8, "Sep": 9, "Oct":10, "Nov":11, "Dec":12}
		# check for type
		if isinstance(month_in, basestring):
			try:
				month_n = month_map[month_in]
			except:
				raise Exception("Unknown month: " + str(month_in) + ".  Use Jan, Feb, Mar, Apr,..., Dec")
		elif isinstance(month_in, int):
			if month_in < 1 or month_in > 12:
				raise Exception("Month number is not between 1 and 12")
			else:
				month_n = month_in
		else:
			raise Exception("Unknown month type")

		# loop through ensemble
		return_members = []
		for e in self.__ensemble_list:
			# get date bounds in daytime format
			t_vals = e.get_t_dim().get_values("daytime")
			# get start and end month
			start_m = t_vals[0].month
			end_m = t_vals[-1].month
			# month_n lies between start_date and end_date of ensemble member
			if month_n >= start_m and month_n <= end_m:
				return_members.append(e)
		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def subset_by_season(self, season_in):
		"""Get the members of the ensemble that occur within a certain season"""
		season_map = {"DJF":[12,1,2], "MAM":[3,4,5], "JJA":[6,7,8], "SON":[9,10,11]}
		# check for type and get list of month numbers containing season
		if isinstance(season_in, basestring):
			try:
				season_list = season_map[season_in]
			except:
				raise Exception("Unknown season " + str(season_in) + ".  Use DJF, MAM, JJA or SON.")
		elif isinstance(season_in, int):
			if season_in < 1 and season_in > 4:
				raise Exception("Unknown season number " + str(season_in) + ".  Use 1=DJF, 2=MAM, 3=JJA or 4=SON.")
			else:
				season_list = season_map[["DJF", "MAM", "JJA", "SON"][season_in-1]]

		return_members = []
		for e in self.__ensemble_list:
			# get date bounds in daytime format
			t_vals = e.get_t_dim().get_values("daytime")
			start_m = t_vals[0].month
			end_m = t_vals[-1].month
			# if start_m or end_m in season_list then append member
			if start_m in season_list or end_m in season_list:
				return_members.append(e)
		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def subset_by_attribute(self, attr_key, attr_value):
		"""Get the members that have the attribute attr_key and attr_value for 
			the key."""
		
		return_members = []
		for e in self.__ensemble_list:
			# get the attributes
			var_attrs = e.get_attributes()
			# if the attr_key exists and the value for the attribute is equal
			# add to the return list
			if attr_key in var_attrs.keys() and var_attrs[attr_key] == attr_value:
				return_members.append(e)
		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def subset_by_global_attribute(self, attr_key, attr_value):
		"""Get the members that have the global attribute attr_key and attr_value
			for the key."""

		return_members = []
		for e in self.__ensemble_list:
			# get the attributes
			glob_attrs = e.get_global_attributes()
			# if the attr_key exists and the value for the attribute is equal
			# add to the return list
			if attr_key in glob_attrs.keys() and glob_attrs[attr_key] == attr_value:
				return_members.append(e)
		return cpdn_ensemble(self.__varname, return_members)

	###########################################################################

	def get_members(self):
		return self.__ensemble_list

	###########################################################################

	def get_size(self):
		"""Return the size of the ensemble."""
		return len(self.__ensemble_list)

	###########################################################################

	def get_length(self):
		return self.get_size()

	###########################################################################

	def get_start_date(self, t_mode="string"):
		"""Get the first date in the ensemble."""
		c_date = self.__ensemble_list[0].get_t_dim().get_values("value")[0]
		sd, tu, n_days_py = self.__ensemble_list[0].get_t_dim().get_time_details()
		for e in self.__ensemble_list[1:]:
			t_date = e.get_t_dim().get_values("value")[0]
			if t_date < c_date:
				c_date = t_date
		if t_mode == "daytime":
			c_date = float_to_daytime(c_date, n_days_py)
		if t_mode == "string":
			c_date = float_to_daytime(c_date, n_days_py).isoformat(" ")
		return c_date

	###########################################################################

	def get_end_date(self, t_mode="string"):
		"""Get the first date in the ensemble."""
		c_date = self.__ensemble_list[0].get_t_dim().get_values("value")[-1]
		sd, tu, n_days_py = self.__ensemble_list[0].get_t_dim().get_time_details()
		for e in self.__ensemble_list[1:]:
			t_date = e.get_t_dim().get_values("value")[-1]
			if t_date > c_date:
				c_date = t_date
		if t_mode == "daytime":
			c_date = float_to_daytime(c_date, n_days_py)
		if t_mode == "string":
			c_date = float_to_daytime(c_date, n_days_py).isoformat(" ")
		return c_date

	###########################################################################

	def get_time_period(self):
		"""Get the number of days between ensemble members, or within ensemble
			members."""
		# get from the bounds
		bnds = self.__ensemble_list[0].get_box().get_dimension("T").get_bounds(t_mode="value")
		tp = bnds[0][1] - bnds[0][0]
		return tp

	###########################################################################

	def save(self, name):
		"""Save as a pickle."""
		fh = open(name, "w")
		# loop through and save each ensemble member
		pickle.dump(self.__varname, fh)
		pickle.dump(len(self.__ensemble_list), fh)	# write out ensemble size
		for e in self.__ensemble_list:
			e.save(fh)
		fh.close()

	###########################################################################

	def load(self, name):
		"""Load from a pickle.  Supports gzipped pickles."""
		if ".gz" in name:
			fh = gzip.open(name, "r")
		else:
			fh = open(name, "r")
		# loop through and save each ensemble member
		self.__varname = pickle.load(fh)
		e_size = pickle.load(fh)		# read ensemble size
		for i in range(0, e_size):
			filename = pickle.load(fh)
			varname  = pickle.load(fh)
			t_dim    = pickle.load(fh)
			var_attrs = pickle.load(fh)
			glob_attrs = pickle.load(fh)
			e_mem = cpdn_ensemble_member(filename, varname, t_dim, 
										 var_attrs, glob_attrs)
			self.__ensemble_list.append(e_mem)
		fh.close()
