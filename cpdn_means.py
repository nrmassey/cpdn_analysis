###############################################################################
# File         : cpdn_means.py
# Author       : Neil Massey
# Created      : 25/05/12
# Purpose      : routines to take a cpdn_box and produce area-weighted means
#                or standard deviations from them
# Changes      : 
###############################################################################

from cpdn_box import *
from cpdn_helpers import *
import numpy
import numpy.ma as ma
import math
from datetime import datetime
from daytime import days_elapsed
from copy import *
from types import *
import random

debug=True

###############################################################################

def calculate_grid_area(lat1, lon1, lat2, lon2, r=1.0):
	"""Calculate the surace area of a grid box on a sphere defined by two points, 
		given as two latitude longitude pairs.
	   Parameters : lat1, lon1 : latitude and longitude of left hand corner
       		      : lat2, lon2 : latitude and longitude of right hand corner
	   Returns    : area in metres squared (on a unit sphere if r=1.0)"""

	lat1_rad = math.radians(lat1)
	lon1_rad = math.radians(lon1)
	lat2_rad = math.radians(lat2)
	lon2_rad = math.radians(lon2)

	area = r**2*math.fabs(lon2_rad-lon1_rad) *\
            math.fabs(math.sin(lat2_rad) - math.sin(lat1_rad))
	return area

###############################################################################

def calc_aa_weights(X_bounds, Y_bounds, mask=None):
	"""Calculate the Area Average weights, given an X dimension's bounds and
	   a Y dimension's bounds.
	   Returns : numpy array of weights."""

	# create the weights area
	weights = numpy.zeros([Y_bounds.shape[0], X_bounds.shape[0]], 'f')

	# now loop through the Y_bounds only and X bounds
	for x in range(0, X_bounds.shape[0]):
		for y in range(0, Y_bounds.shape[0]):
			weights[y,x] = calculate_grid_area(Y_bounds[y][0], X_bounds[x][0],
								    		   Y_bounds[y][1], X_bounds[x][1])

	# check mask and multiply if it exists
	if mask != None:
		weights = weights * mask

	# normalise the weights
	weights = weights / numpy.max(weights)
	return weights

###############################################################################

# helper functions for creating the output data

def create_tgt_dims(box, meaned_axes):
	# create the target dimensions
	tgt_dims = []
	for sd in box.get_dimensions():
		if not sd.get_axis() in meaned_axes:
			tgt_dims.append(sd)
		else:
			# create a dimension of length one but incorporating the extents of
			# the original dimension as the bounds
			sdv = sd.get_values("true_value")
			
			tgt_dim_vals = numpy.array([sdv[0] + 0.5*(sdv[-1]-sdv[0])])
			tgt_dim_bnds = numpy.array([[sdv[0], sdv[-1]]])
			# create the dimensions
			if sd.get_axis() == "T":
				start_date, time_units, n_days = sd.get_time_details()
				tgt_dim_vals = tgt_dim_vals.astype(int)
			else:
				start_date, time_units, n_days = None, None, None
			# add the bounds to the attributes
			tgt_attrs = sd.get_attributes()
			tgt_attrs["bounds"] = sd.get_name() + "_bounds"
			tgt_dim = cpdn_boxdim(name=sd.get_name(), vals=tgt_dim_vals, 
								  attrs=sd.get_attributes(), axis=sd.get_axis(),
								  bounds=tgt_dim_bnds, start_date=start_date,
								  time_units=time_units, n_days_per_year=n_days)
			tgt_dims.append(tgt_dim)
	return tgt_dims

###############################################################################

def cpdn_global_mean(box, mask=None):
	"""Create an area-weighted global mean from a cpdn_box.  Weights are
	   automatically calculated from the box dimensions.  The meaning will always
	   work on the X & Y dimensions, operating over the other dimensions.  For
	   example, if the box contains a time dimensions, then a time-series will 
	   be created.  If the box contains a Z dimension then a multi-level 
	   time-series will be created.
	   Returns : a box containing the area-average means"""

	# get the the X and Y bounds
	X_bounds = box.get_dimension("X").get_bounds()
	Y_bounds = box.get_dimension("Y").get_bounds()
	# calculate the weights from these bounds
	wts = calc_aa_weights(X_bounds, Y_bounds, mask)
	wts = wts.squeeze()

	# get the position of the X and Y axes
	dim_axes = box.get_dimension_axes()
	X_axis = dim_axes.index("X")
	Y_axis = dim_axes.index("Y")

	# get the dimensions to iterate over
	it_dims = []
	dim_lens = []
	for d in box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)
			dim_lens.append(d.get_len())

	# create the target data
	tgt_data = numpy.zeros(dim_lens, 'f')

	# create an iterator
	box_it = cpdn_box_iterator(box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)

	# iterate over the box dimensions that are not X and Y
	while not box_it.end():
		src_data = box[src_idx].get_values().squeeze()
		tgt_data[tgt_idx] = numpy.average(src_data, weights=wts)
		src_idx, tgt_idx = box_it.next(True)

	tgt_dims = create_tgt_dims(box, ["X","Y"])

	# now create a box to put the data in - all dimensions except X & Y, same attributes
	# as original
	# get the attributes first and the scale factor and offset (if any) from them
	method_str = "longitude: latitude: area-weighted mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: area-weighted mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=box.get_name(),
					   glob_attrs=tgt_glob_attrs, off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_zonal_mean(box, mask=None):
	"""Create a zonal mean from a cpdn_box. The meaning will always
       work on the X dimension, operating over the other dimensions.  For
       example, if the box contains a time dimensions, then a time-series will 
       be created.  If the box contains a Z dimension then a multi-level 
       time-series will be created.
       Returns : a box containing the zonal means"""

	# No need to calculate weights as the Y axis is not meaned
	# get the axis number used for the "X" dimension
	X_axis = box.get_dimension_axes().index("X")
	src_data = box.__getitem__().get_values()
	tgt_data = numpy.average(src_data, axis=X_axis)

	# create target dimensions
	tgt_dims = create_tgt_dims(box, ["X"])

	# now create a box to put the data in - all dimensions except X
	# get the attributes first and the scale factor and offset (if any) from them
	method_str = "longitude: zonal mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: zonal mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)

	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=box.get_name(),
					   off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_meridional_mean(box, mask=None):
	"""Create a meridional mean from a cpdn_box. The meaning will always
       work on the Y dimension, operating over the other dimensions.  For
       example, if the box contains a time dimensions, then a time-series will 
       be created.  If the box contains a Z dimension then a multi-level 
       time-series will be created.
       Returns : a box containing the meridional means"""

	# get the the X and Y bounds - only need the first X bounds as the meaning
	# is over a longitude line
	X_bounds = box.get_dimension("X").get_bounds()[0:1]
	Y_bounds = box.get_dimension("Y").get_bounds()
	# calculate the weights from these bounds
	wts = calc_aa_weights(X_bounds, Y_bounds, mask)

	# squeeze the weights as we will be squeezing the data later - no need to
	# transpose in this case
	wts = wts.squeeze()
	# get the source data
	src_data = box.__getitem__().get_values()

	# get the Y_axis index
	Y_axis = box.get_dimension_axes().index("Y")
	# create the target data
	tgt_data = numpy.average(src_data, axis=Y_axis, weights=wts)
	# create the target dimensions
	tgt_dims = create_tgt_dims(box, ["Y"])

	# now create a box to put the data in - all dimensions except X
	# amend the attributes to add the mean to the cell methods and to update the history
	method_str = "latitude: meridional mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: meridional mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=box.get_name(),
					   glob_attrs=tgt_glob_attrs, off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_temporal_mean(box):
	"""Create a temporal mean from a cpdn_box."""
	# get the axis number used for the "T" dimension
	T_axis = box.get_dimension_axes().index("T")
	src_data = box.__getitem__().get_values().squeeze()
	tgt_data = numpy.average(src_data, axis=T_axis)

	# create target dimensions
	tgt_dims = create_tgt_dims(box, ["T"])

	# now create a box to put the data in - all dimensions except T
	# get the attributes first and the scale factor and offset (if any) from them
	# amend the attributes to add the mean to the cell methods and to update the history
	method_str = "time: mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: temporal mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, glob_attrs=tgt_glob_attrs,
					   name=box.get_name(), off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_period_mean(box, period):
	# get the axis number used for the "T" dimension
	T_axis = box.get_dimension_axes().index("T")
	src_data = box.__getitem__().get_values().squeeze()
	src_times = box.get_dimension("T").get_values("true_values")
	src_T_dim = box.get_dimension("T")

	# create target dimension sizes
	tgt_data_size = list(src_data.shape)
	tgt_data_size[T_axis] = tgt_data_size[T_axis]/period
	# create target data
	tgt_data = numpy.zeros(tgt_data_size, 'f')
	tgt_times = numpy.zeros([tgt_data_size[T_axis]], 'f')

	# create a source index
	src_idx = []
	tgt_idx = []
	for d in range(0, src_data.ndim):
		src_idx.append(slice(None, None, None))
		tgt_idx.append(slice(None, None, None))

	# amend the T dimension
	sd, tu, nd = src_T_dim.get_time_details()
	# now loop over the data
	for t in range(0, tgt_data_size[T_axis]):
		t0 = t * period
		t1 = t0 + period
		src_idx[T_axis] = slice(t0, t1)
		tgt_idx[T_axis] = t
		tgt_data[tgt_idx] = numpy.mean(src_data.__getitem__(src_idx), axis=T_axis)
		tgt_times[t] = src_times[t0] * tu + (src_times[t1-1] - src_times[t0]) * tu / 2.0 - sd

	# create the target dimensions
	tgt_dims = create_tgt_dims(box, ["T"])
	tgt_T_dim = cpdn_boxdim(name=src_T_dim.get_name(), vals=tgt_times, 
						    attrs=src_T_dim.get_attributes(), axis="T", 
 			                start_date=sd, time_units=tu, n_days_per_year=nd)
	tgt_dims[T_axis] = tgt_T_dim

	# amend the attributes to add the mean to the cell methods and to update the history
	method_str = "time: period mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: temporal period mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=box.get_name(),
					   glob_attrs=tgt_glob_attrs, off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_temporal_sum(box):
	"""Create a temporal mean from a cpdn_box."""
	# get the axis number used for the "T" dimension
	T_axis = box.get_dimension_axes().index("T")
	src_data = box.__getitem__().get_values()
	tgt_data = numpy.sum(src_data, axis=T_axis)

	# create target dimensions
	tgt_dims = create_tgt_dims(box, ["T"])

	# now create a box to put the data in - all dimensions except T
	# get the attributes first and the scale factor and offset (if any) from them
	# amend the attributes to add the mean to the cell methods and to update the history
	method_str = "time: sum"
	history_str = datetime.now().isoformat() + " altered by CPDN: temporal sum."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, glob_attrs=tgt_glob_attrs,
					   name=box.get_name(), off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_monthly_mean(box):
	"""Create a monthly mean from a cpdn_box, 
	   e.g. mean of Jan 1969, mean of Feb 1969, ..., mean of Dec 1969
	        mean of Jan 1970, mean of Feb 1970, ..., etc."""
	T_axis = box.get_dimension_axes().index("T")
	T_dim = box.get_dimension("T")
	T_vals = T_dim.get_values("daytime")
	M_sd = T_vals[0].month
	Y_sd = T_vals[0].year
	sd, tu, nd = T_dim.get_time_details()
	# get missing value
	mv = box.get_missing_value()
	# loop over the data
	c_pos = 0
	# create the array
	tgt_shape = []
	new_t_len = int(T_dim.get_len() / nd) * 12
	for d in box.get_dimensions():
		if d.get_axis() != "T":
			tgt_shape.append(d.get_len())
		else:
			tgt_shape.append(new_t_len)
	tgt_data = numpy.ones(tgt_shape, 'f') * -mv
	
	# do the loop
	t_pos = 0
	while c_pos < len(T_vals) and t_pos < new_t_len:
		s_pos = c_pos
		# find the start and end of the month
		while c_pos < len(T_vals) and T_vals[c_pos].month == T_vals[s_pos].month:
			c_pos += 1
		# do the mean
		data = ma.masked_equal(box[s_pos:c_pos-1].get_values(), mv)
		tgt_data[t_pos] = data.mean(axis=T_axis).filled(mv)
		t_pos += 1
	
	# create the target time values and the bounds
	vals = []
	bnds = []
	c_yr = Y_sd
	c_mn = M_sd
	for m in range(0, tgt_shape[T_axis]):
		if nd == 360:
			DV = m*30*un+15 + T_sd
			vals.append(DV)
			bnds.append([DV-15, DV+15])
		else:
			if (c_mn == 13):
				c_yr += 1
				c_mn = 1
			DM = int(0.5 + float(days_elapsed[c_mn] - days_elapsed[c_mn-1])/2)
			cd = daytime(c_yr, c_mn, DM)
			DV = daytime_to_float(cd) - sd
			if c_yr % 4 == 0 and c_mn >= 3:
				DV += 1
			vals.append(DV)
			bnds.append([DV-DM, DV+DM])
			c_mn += 1

	method_str = "time: monthly mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: monthly mean."
	tgt_box = cpdn_clone_box(box, new_data=tgt_data, T_vals=vals, T_bnds=bnds, 
				   method_string=method_str, history_string=history_str)
	return tgt_box
	
###############################################################################

def cpdn_climatological_month_mean(box):
	"""Create a climatological month mean from a cpdn_box. i.e. all Jans meaned,
	   all Febs meaned etc."""
	# get the time dimension and its values
	T_axis = box.get_dimension_axes().index("T")
	T_vals = box.get_dimension("T").get_values("daytime")
	time_vals_list = [[] for m in range(0, 12)]
	# get missing value
	mv = box.get_missing_value()
	# get missing value indices - assume same across all time periods
	mv_idx = numpy.where(box[0].get_values().squeeze() == mv)

	# loop through the T_vals and build 12 lists (one for each month) of time values
	for t_val in T_vals:
		time_vals_list[t_val.month-1].append(t_val)

	# now produce the month means, need an array to store them in - first 
	# determine shape
	tgt_shape = []
	for d in box.get_dimensions():
		if d.get_axis() != "T":
			tgt_shape.append(d.get_len())
		else:
			tgt_shape.append(12)
	tgt_data = numpy.zeros(tgt_shape, 'f')

	# now loop through and add to the target data, weighted by 1.0/number of jans etc.
	src_idx = [slice(None, None, None) for d in range(0,box.get_ndims())]
	for m in range(0, 12):
		mw = 1.0 / len(time_vals_list[m])
		for sm in range(0, len(time_vals_list[m])):
			src_idx[T_axis] = time_vals_list[m][sm]
			v = numpy.mean(box.__getitem__(src_idx).get_values(), axis=T_axis)
			tgt_data[m] += mw * v
			# reinstate missing values
			tgt_data[m][mv_idx] = mv

	# create the target dimensions
	tgt_dims = []
	start_year = T_vals[0].year
	end_year   = T_vals[-1].year
	for d in box.get_dimensions():
		if d.get_axis() == "T":
			# get the start date, units and number of days in the year
			sd, un, nd = d.get_time_details()
			# create the values and the bounds
			vals = []
			for m in range(0, 12):
				if nd == 360:
					vals.append(m*30*un)
				else:
					vals.append(days_elapsed[m]*un)
			bnds = numpy.array([[vals[x], vals[x]+(end_year-start_year)*nd+30] for x in range(0, 12)])
			# create the time dimension and append
			td = cpdn_boxdim(name=d.get_name(), vals=vals, attrs=d.get_attributes(), 
							 axis="T", bounds=bnds, start_date=sd, time_units=un, 
							 n_days_per_year=nd)
			tgt_dims.append(td)
		else:
			tgt_dims.append(d)

	# create the return box
	method_str = "time: climatological month mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: climatological month mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, glob_attrs=tgt_glob_attrs,
					   name=box.get_name(), off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_climatological_seasonal_mean(box):
	"""Create a climatological seasonal mean from a cpdn_box. i.e. all DJF meaned,
	   all MAM, all JJA, all SON meaned."""
	# get the time dimension and its values
	T_axis = box.get_dimension_axes().index("T")
	T_vals = box.get_dimension("T").get_values("daytime")
	time_vals_list = [[] for m in range(0, 4)]

	# loop through the T_vals and build 4 lists (one for each month) of time values
	for t_val in T_vals:
		# DJF
		if t_val.month in [12, 1, 2]:
			time_vals_list[0].append(t_val)
		# MAM
		elif t_val.month in [3, 4, 5]:
			time_vals_list[1].append(t_val)
		# JJA
		elif t_val.month in [6, 7, 8]:
			time_vals_list[2].append(t_val)
		# DJF
		elif t_val.month in [9, 10, 11]:
			time_vals_list[3].append(t_val)		

	# now produce the seasonal means, need an array to store them in - first 
	# determine shape
	tgt_shape = []
	z_dim = "Z" in box.get_dimension_axes()
	for d in box.get_dimensions():
		if d.get_axis() != "T":
			tgt_shape.append(d.get_len())
		else:
			tgt_shape.append(4)
	tgt_data = numpy.zeros(tgt_shape, 'f')
	tgt_sum = numpy.zeros(tgt_shape, 'f')
	mv = box.get_missing_value()

	# now loop through and add to the target data, weighted by 1.0/number of jans etc.
	src_idx = [slice(None, None, None) for d in range(0,box.get_ndims())]
	for m in range(0, 4):
		for sm in range(0, len(time_vals_list[m])):
			# get the data
			src_idx[T_axis] = time_vals_list[m][sm]
			data = box.__getitem__(src_idx).get_values().squeeze()
			# get where the data does not equal the mv
			non_mv_idx = numpy.where(data != mv)
			# add the data and sum
			if z_dim:
				tgt_data[m][0][non_mv_idx] += data[non_mv_idx]
				tgt_sum[m][0][non_mv_idx] += 1
			else:
				tgt_data[m][non_mv_idx] += data[non_mv_idx]
				tgt_sum[m][non_mv_idx] += 1
		# divide through by the sum to get the mean
		tgt_data[m] = tgt_data[m] / tgt_sum[m]
		# set where the sum is 0 to be the mv
		if not numpy.isinf(mv):
			tgt_data[m][tgt_sum[m] == 0] = mv

	# create the target dimensions
	tgt_dims = []
	start_year = T_vals[0].year
	end_year   = T_vals[-1].year
	for d in box.get_dimensions():
		if d.get_axis() == "T":
			# get the start date, units and number of days in the year
			sd, un, nd = d.get_time_details()
			# create the values and the bounds
			vals = [(x*90*un)+15 for x in range(0, 4)]
			bnds = numpy.array([[vals[x], vals[x]+(end_year-start_year)*nd+90*un] for x in range(0, 4)])
			# create the time dimension and append
			td = cpdn_boxdim(name=d.get_name(), vals=vals, attrs=d.get_attributes(),
							 axis="T", bounds=bnds, start_date=sd, time_units=un, 
							 n_days_per_year=nd)
			tgt_dims.append(td)
		else:
			tgt_dims.append(d)

	# create the return box
	method_str = "time: climatological seasonal mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: climatological seasonal mean."
	tgt_attrs = amend_attributes(box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, glob_attrs=tgt_glob_attrs,
					   name=box.get_name(), off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def cpdn_ensemble_mean(ensemble):
	"""Create an ensemble mean from a cpdn_ensemble.  Each ensemble member is
		equally weighted."""

	# get the start date, end date and time period from the ensemble
	sd = ensemble.get_start_date(t_mode="value")
	ed = ensemble.get_end_date(t_mode="value")
	pd = ensemble.get_time_period()

	out_data = None
	time_data = []

	# iterate over the days
	for d in numpy.arange(sd, ed+pd, pd):
		sub_ens = ensemble.subset_by_date(d)
		ens_mems = sub_ens.get_members()
		# need the start date of the ensemble member
		sde, tue, n_days_py = ens_mems[0].get_t_dim().get_time_details()
		if debug==True:	# debug output
			print "Date " + float_to_daytime(d, n_days_py).isoformat(" ") +\
				  " ensemble members " + str(len(ens_mems))
		# get the weight - 1.0 / number of ensemble members
		ew = 1.0 / len(ens_mems)
		# get the time dimension axis number
		T_axis = ens_mems[0].get_box().get_dimension_axes().index("T")
		# now do for the rest of the ensemble members - but adding to the weighted sum
		c_data = None
		for e in range(0, len(ens_mems)):
			# let the current index be the ensemble member
			c_idx = e
			# set a flag that a value hasn't been assigned
			val_assigned = False
			# keep going until a value is assigned
			while not val_assigned:
				# try getting a value
				try:
					val = ew * ens_mems[c_idx][d-sde].get_values()
					val_assigned = True
				except:
					# if there is an exception then the value could not be assigned
					# and so the ensemble member is invalid
					# try again with a random ensemble member
					c_idx = int(random.uniform(0, len(ens_mems)))

			# assign to c_data
			if isinstance(c_data, NoneType):
				c_data = val
			else:
				c_data += val

		# do we have to create the data (in the timeseries) or just append to it?
		if isinstance(out_data, NoneType):
			out_data = c_data
		else:
			out_data = numpy.concatenate([out_data, c_data], axis=T_axis)
		time_data.append(d)

	# now create a box containing the data and the new T axis
	# use first ensemble members box as the template
	src_box = ensemble.subset_by_date(sd).get_members()[0].get_box()
	# get the dimensions - and find the time axis dimension
	tgt_dims = src_box.get_dimensions()
	T_axis = src_box.get_dimension_axes().index("T")

	##### steps to create a new time dimension
	# get the time units, start date and number of days per year - start date is 0, though
	sdate, timeu, n_days_py = ens_mems[0].get_t_dim().get_time_details()
	# get the old attributes and modify the source date
	time_attrs = tgt_dims[T_axis].get_attributes()
	# get the start_date as a daytime
	sd_dt = float_to_daytime(sd, n_days_py)
	time_attrs['units'] = "days since " + sd_dt.isoformat(" ")
	# minus sd from time_data after converting to numpy array
	time_data = numpy.array(time_data)-sd

	# create a new time dimension and overwrite the existing dimension
	tgt_t_dim = cpdn_boxdim(name=tgt_dims[T_axis].get_name(), 
							vals=time_data, attrs=time_attrs, axis="T", 
							start_date=sd, time_units=timeu, n_days_per_year=n_days_py)
	tgt_dims[T_axis] = tgt_t_dim

	#####

	# modify the variables attributes - add a cell method and history
	method_str = "ensemble mean"
	history_str = datetime.now().isoformat() + " altered by CPDN: ensemble mean."
	tgt_attrs = amend_attributes(src_box.get_attributes(), method_str, history_str)
	tgt_glob_attrs = src_box.get_global_attributes()

	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	# create the box and return
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=src_box.get_name(),
						glob_attrs=tgt_glob_attrs, off=tgt_off, sf=tgt_sf, data=out_data)

	return tgt_box
