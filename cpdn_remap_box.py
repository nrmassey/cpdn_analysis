###############################################################################
# File         : cpdn_remap_box.py
# Author       : Neil Massey
# Created      : 28/08/12
# Purpose      : routine to remap longitude between boxes - i.e. convert from
#                0 to 360 to -180 to 180 and vice versa
#				 also remap longitude from 90 to -90 to -90 to 90 and vice versa
# Changes      : 
###############################################################################

from cpdn_box import *

###############################################################################

def cpdn_remap_longitude(box):
	# get the box longitude and decide what to do
	LON_dim = box.get_dimension("X")
	LON = LON_dim.get_values()
	
	# create an iterator over the source box
	it_dims = []
	for d in box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

	# create an iterator
	box_it = cpdn_box_iterator(box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)

	# create the target data - same shape as the source box
	tgt_data = numpy.zeros(box.get_dimension_lengths(),'f')

	# get where the longitude is less or greater than the date line / meridion
	if LON[0] < 0.0:
		lon_under_0 = numpy.where(LON < 0)
		s = lon_under_0[0][0]
		e = lon_under_0[0][-1] + 1
		LON[lon_under_0] += 360
		LON.sort()
		under_0 = True
		p = box.get_dimension("X").get_length() - e
	else:
		lon_over_180 = numpy.where(LON > 180)
		s = lon_over_180[0][0]
		e = lon_over_180[0][-1] + 1
		LON[lon_over_180] -= 360
		LON.sort()
		under_0 = False
		p = box.get_dimension("X").get_length() - s

	# loop through
	while not box_it.end():
		data = box[src_idx].get_values().squeeze()
		if under_0:
			# remap to 0 to 360
			tgt_full_idx = list(tgt_idx)
			tgt_full_idx.append(slice(None, None, None))
			tgt_full_idx.append(slice(p, None, None))
			tgt_data[tgt_full_idx] = data[:,s:e]
			tgt_full_idx = list(tgt_idx)
			tgt_full_idx.append(slice(None, None, None))
			tgt_full_idx.append(slice(0, p, None))
			tgt_data[tgt_full_idx]  = data[:,e:]
		else:
			# remap to -180 to 180
			tgt_full_idx = list(tgt_idx)
			tgt_full_idx.append(slice(None, None, None))
			tgt_full_idx.append(slice(0, p, None))
			tgt_data[tgt_full_idx] = data[:,s:e]
			tgt_full_idx = list(tgt_idx)
			tgt_full_idx.append(slice(None, None, None))
			tgt_full_idx.append(slice(p, None, None))
			tgt_data[tgt_full_idx]  = data[:,0:s]

		# move onto next
		src_idx, tgt_idx = box_it.next(True)

	# create the output dimensions
	out_dims = []
	for d in box.get_dimensions():
		axis = d.get_axis()
		if axis == "X":
			lon_dim = cpdn_boxdim(LON_dim.get_name(), LON, LON_dim.get_attributes(), "X")
			out_dims.append(lon_dim)
		else:
			out_dims.append(d)

	# return a box
	out_box = cpdn_box(name=box.get_name(), dims=out_dims, var_attrs=box.get_attributes(),
					   glob_attrs = box.get_global_attributes(), off=0.0, sf=1.0,
					   data=tgt_data)
	return out_box

###############################################################################

def cpdn_remap_latitude(box):
	# create an iterator over the source box
	it_dims = []
	for d in box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

	# create an iterator
	box_it = cpdn_box_iterator(box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)

	# create the target data - same shape as the source box
	tgt_data = numpy.zeros(box.get_dimension_lengths(),'f')

	# loop through
	while not box_it.end():
		data = box[src_idx].get_values().squeeze()
		# flip in the Y direction
		tgt_data[tgt_idx]  = data[::-1,:]
		# move onto next
		src_idx, tgt_idx = box_it.next(True)

	# get the lat dim and flip the lat values
	LAT_dim = box.get_dimension("Y")
	LAT_dim_vals = LAT_dim.get_values()[::-1]
	
	# create the output dimensions
	out_dims = []
	for d in box.get_dimensions():
		axis = d.get_axis()
		if axis == "Y":
			lat_dim = cpdn_boxdim(LAT_dim.get_name(), LAT_dim_vals, LAT_dim.get_attributes(), "Y")
			out_dims.append(lat_dim)
		else:
			out_dims.append(d)

	# return a box
	out_box = cpdn_box(name=box.get_name(), dims=out_dims, var_attrs=box.get_attributes(),
					   glob_attrs = box.get_global_attributes(), off=0.0, sf=1.0,
					   data=tgt_data)
	return out_box

###############################################################################

def cpdn_replace_mv(box, new_mv):
	# replace a missing_value with a new missing value
	# this is useful if the missing value isn't massively massive (e.g. -54 instead of 2e20)
	# get data values
	src_data = box.get_values()
	# get missing value and indices into data
	mv = box.get_missing_value()
	mv_idx = numpy.where(src_data == mv)
	# replace missing values
	src_data[mv_idx] = new_mv

	# create new box
	var_attrs = box.get_attributes()
	if not isinstance(var_attrs, NoneType):
		if "missing_value" in var_attrs.keys():
			var_attrs["missing_value"] = new_mv
		if "_FillValue" in var_attrs.keys():
			var_attrs["_FillValue"] = new_mv
	
	out_box = cpdn_box(name=box.get_name(), dims=box.get_dimensions(), 
					   var_attrs=var_attrs, glob_attrs = box.get_global_attributes(), 
					   off=0.0, sf=1.0, data=src_data)
	return out_box

###############################################################################

def cpdn_reshape(box, new_shape, new_dim_names=[], new_dim_axes=[]):
	"""Reshape the box to take on the new shape - will create dummy dimensions if
	   neccessary."""
	# get the values
	src_data = box.get_values("true_value")
	# check that the size of the arrays match
	if numpy.product(new_shape) != numpy.product(src_data.shape):
		raise Exception("Box cannot be refactored into new shape - size differs.")

	# match dimensions up by length
	dls = box.get_dimension_lengths()
	box_dims = box.get_dimensions()
	dimension_map = [[-1,-1] for s in new_shape]
	for s in range(0, len(new_shape)):
		if new_shape[s] in dls:
			i = dls.index(new_shape[s])
			if not i in dimension_map:			# only add once
				dimension_map[s] = dls.index(new_shape[s])		# create a map between dimension and new shape
			else:
				dimension_map[s] = -1
		else:
			dimension_map[s] = -1

	# now create the new dimensions
	new_dimensions = []
	n_new_dims = 0
	for d in dimension_map:
		if d == -1:
			# create a new dimension
			if new_dim_names != []:
				nd_name = new_dim_names[n_new_dims]
			else:
				nd_name = "null_dim_"+str(n_new_dims)
			if new_dim_axes != []:
				nd_axis = new_dim_axes[n_new_dims]
			else:
				nd_axis = "N"
			new_dim = cpdn_boxdim(name=nd_name, vals=[0.0], attrs={}, axis=nd_axis)
			new_dimensions.append(new_dim)
			n_new_dims += 1
		else:
			new_dimensions.append(box_dims[d])

	# reshape the data
	tgt_data = numpy.reshape(src_data, new_shape)

	# create the box
	tgt_box = cpdn_box(name = box.get_name(), dims = new_dimensions,
					   var_attrs = box.get_attributes(), 
					   glob_attrs = box.get_global_attributes(),
					   off = box.get_off(), sf = box.get_sf(), data = tgt_data)
	return tgt_box
