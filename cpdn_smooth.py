###############################################################################
# File         : cpdn_smooth.py
# Author       : Neil Massey
# Created      : 26/06/12
# Purpose      : routines to smooth data contained in a cpdn_box
# Changes      : 
###############################################################################

from cpdn_box import *
from copy import *
import numpy
from scipy.signal import convolve2d
from cpdn_regrid import fill_missing_data_field

###############################################################################

def get_smoothing_window(j, i, j_len, i_len, w_x, w_y):
	# calculate the smoothing window indices
	idx_list = []
	for y in range(-w_y, w_y+1):
		for x in range(-w_x, w_x+1):
			y0 = j + y
			x0 = i + x
			# wrap around for x
			if x0 < 0:
				x0 += i_len
			if x0 >= i_len:
				x0 -= i_len
			# but only add if j >= 0 and j < j_len
			if y0 >= 0 and y0 < j_len:
				idx_list.append([y0,x0])
	return idx_list

###############################################################################

def cpdn_spatial_smooth(src_box, w_x=3, w_y=None, type="flat"):
	"""Spatially smooth (in X and Y dimensions) data contained within a box."""
	# create or copy the window
	if isinstance(w_x, numpy.ndarray):
		win = w_x / numpy.sum(w_x)
	elif isinstance(w_x, int):
		if isinstance(w_y, NoneType):
			w_y = w_x
		if type == "flat":
			win = numpy.ones([w_x*2+1,w_y*2+1], 'f') / (w_x*w_y*4)
			win = win / numpy.sum(win)
		elif type == "gauss":
			x,y = numpy.mgrid[-w_x:w_x+1, -w_y:w_y+1]
			win = numpy.exp(-(x**2/float(w_x)+y**2/float(w_y)))
			win = win / numpy.sum(win)
	else:
		raise Exception("Unknown window type")

	# create the target data - same shape as the source box
	tgt_data = numpy.zeros(src_box.get_dimension_lengths(),'f')

	# create an iterator over the source box
	it_dims = []
	for d in src_box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

    # create an iterator
	box_it = cpdn_box_iterator(src_box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)

	# get missing value - don't want to convolve this
	mv = src_box.get_missing_value()

	# loop through
	while not box_it.end():
		# get the source data
		src_data = src_box[src_idx].get_values().squeeze()
		# check whether the data can be smoothed quickly using convolve
		if not numpy.isinf(mv):
			mv_in_src = numpy.where(numpy.abs(src_data) > mv*0.9)[0].shape[0] > 0
		else:
			mv_in_src = False
		if not numpy.isinf(mv) and not mv_in_src:
			# convolve the data with the window - ensure it's the same size and wrap around
			# the date line (and the poles, unfortunately)
			c_data = convolve2d(win, src_data, mode="same", boundary="wrap")
			tgt_data[tgt_idx] = numpy.reshape(c_data, src_data.shape)
		else:
			# do the slow dance with missing data
			win_f = win.flatten()
			for j in range(0, src_data.shape[0]):
				for i in range(0, src_data.shape[1]):
					if abs(src_data[j,i]) > 1000:	#abs(mv*0.9):
						tgt_data[tgt_idx][j,i] = mv
					else:
						# get the list of indices
						idx_list = get_smoothing_window(j, i, 
														src_data.shape[0], 
														src_data.shape[1], 
														w_x, w_y)
						# average using the list of indices
						sum = 0.0
						n = 0
						for k in range(0, len(idx_list)):
							idx = idx_list[k]
							y = idx[0]
							x = idx[1]
							v = src_data[y,x]
							if abs(v) < 1000: #abs(mv*0.9):
								sum += v * win_f[k]
								n += win_f[k]
						if n == 0:
							tgt_data[tgt_idx][j,i] = mv
						else:
							tgt_data[tgt_idx][j,i] = sum/n

		# iterate to next time / z / etc.
		src_idx, tgt_idx = box_it.next(True)
	
	# smoothing does weird things to the poles (due to wrap around boundary conditions)
	# so reinstate the poles
	Y_axis = src_box.get_dimension_axes().index("Y")
	idxs = []
	for x in range(0, src_box.get_ndims()):
		idxs.append(slice(None, None, None))
	# north pole / top
	idxs[Y_axis] = 0
	tgt_data[idxs] = src_box.get_values()[idxs]
	# south pole / bottom
	idxs[Y_axis] = -1
	tgt_data[idxs] = src_box.get_values()[idxs]

	# create return box, all the same except for attributes and history
	method_str = "longitude: latitude: kernel smoothing"
	history_str = datetime.now().isoformat() + " altered by CPDN: area-weighted mean."
	tgt_attrs = amend_attributes(src_box.get_attributes(), method_str, history_str)
	
	# create return box
	tgt_box = cpdn_box(dims=src_box.get_dimensions(), var_attrs=tgt_attrs, 
					   name=src_box.get_name(), glob_attrs=src_box.get_global_attributes(),
                       off=0.0, sf=1.0, data=tgt_data)
	return tgt_box
