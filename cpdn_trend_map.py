###############################################################################
# File         : cpdn_trend_map.py
# Author       : Neil Massey
# Created      : 31/07/12
# Purpose      : routines to derive trend maps from a box
# Changes      : 
###############################################################################

import numpy
from cpdn_box import *
from cpdn_smooth import *
from cpdn_means import *

###############################################################################

def derive_trend(x_data, y_data, period=-1):
	# fit a straight line (1 degree of freedom)
	pl = numpy.polyfit(x_data, y_data, 1)
	x0 = x_data[0]
	x1 = x_data[-1]
	if period != -1:
		ptime = period
	else:
		ptime = x1 - x0
	# get the start and end value of the line and subtact them
	v0 = numpy.polyval(pl, x0)
	v1 = numpy.polyval(pl, x1)
	out_v = (v1 - v0) / ptime
#	print out_v, v1,v0, ptime
	return out_v, v0, v1

###############################################################################

def create_trend_map(box, z_lev=0, mv=None):
	# do we have a Z dimension?
	z_dim = "Z" in box.get_dimension_axes()
	# get the lengths of the X and Y axis
	X_len = box.get_dimension("X").get_len()
	Y_len = box.get_dimension("Y").get_len()
	# create the return data
	if z_dim:
		ret_data = numpy.zeros([1,1,Y_len,X_len], 'f')
	else:
		ret_data = numpy.zeros([1,Y_len,X_len], 'f')
	# loop through each grid point and get the timeseries
	for y in range(0, Y_len):
		for x in range(0, X_len):
			if z_dim:
				gp_ts = box[:,z_lev,y,x].get_values().squeeze()
			else:
				gp_ts = box[:,y,x].get_values().squeeze()
			# create the x axis and data for for the polyfit 
			# lists to store x axis and data
			x_axis = []
			data = []
			for i in range(0, gp_ts.shape[0]):
				if not isinstance(mv, NoneType):
					if abs(gp_ts[i]) < abs(mv):
						x_axis.append(i)
						data.append(gp_ts[i])
				else:
					x_axis.append(i)
					data.append(gp_ts[i])
			# check that there is a fair distribution of data points over the time series
			if len(x_axis) > 0 and (x_axis[-1] - x_axis[0]) > 0.5 * gp_ts.shape[0]:
				out_v = derive_trend(x_axis, data)
			else:
				if mv != None:
					out_v = (mv, 0.0, 0.0)
				else:
					out_v = (0.0, 0.0, 0.0)
			if z_dim:
				ret_data[0,0,y,x] = out_v[0]
			else:
				ret_data[0,y,x] = out_v[0]
	return ret_data

###############################################################################

def create_seasonal_trend_map(box, z_lev=0, mv=None):
	# box must contain a timeseries of seasonal mean maps
	# get the lengths of the X and Y axis
	X_len = box.get_dimension("X").get_len()
	Y_len = box.get_dimension("Y").get_len()
	# 
	# create the output data
	out_data = numpy.zeros([4,1,Y_len,X_len], 'f')
	#
	for s in range(0,4):        # seasons count
		s_ts = box[s::4] # get a season's time series, e.g. DJF
		out_data[s] = create_trend_map(s_ts, z_lev, mv)
	return out_data

