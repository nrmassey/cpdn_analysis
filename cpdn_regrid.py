###############################################################################
# File         : cpdn_regrid.py
# Author       : Neil Massey
# Created      : 01/06/12
# Purpose      : routines to regrid data from a cpdn_box to a new grid
#				 creates a new cpdn_box from the regridded data
#                Upscaling to a higher resolution grid will intepolate with a
#                user selected method - either nearest neighbour, linear or
#                bi-cubic
#                Downscaling to a lower resolution grid will use area-weighted
#                averaging
# Changes      : 
###############################################################################

from cpdn_box import *
from cpdn_means import calculate_grid_area
from copy import *
from types import *
from datetime import datetime
from scipy.ndimage.interpolation import map_coordinates as interp2d

###############################################################################

def create_tgt_box(src_box, tgt_data, tgt_X, tgt_Y, method_str, history_str, mv=numpy.inf):
	# helper function to create the tgt output box

	# get the length of the dimensions first
	src_X = src_box.get_dimension("X").get_values()
	src_Y = src_box.get_dimension("Y").get_values()

	# create the output box
	tgt_attrs = amend_attributes(src_box.get_attributes(), method_str, history_str)

	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0

	# copy the dimensions
	tgt_dims = []
	for d in src_box.get_dimensions():
		if d.get_axis() == "X":
			# create a new dimension with the tgt_X dimension values
			X_dim = cpdn_boxdim(name=d.get_name(), vals=tgt_X, attrs=d.get_attributes(),
								axis="X")
			tgt_dims.append(X_dim)
		elif d.get_axis() == "Y":
			# same as above for Y dim
			Y_dim = cpdn_boxdim(name=d.get_name(), vals=tgt_Y, attrs=d.get_attributes(),
								axis="Y")
			tgt_dims.append(Y_dim)
		else:
			tgt_dims.append(d)

	# amend the missing value
	if not numpy.isinf(mv):
		tgt_attrs["missing_value"] = mv

	# create the target box to return into
	tgt_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, name=src_box.get_name(),
                       off=tgt_off, sf=tgt_sf, data=tgt_data)
	return tgt_box

###############################################################################

def get_index_in_grid(lon, lat, X, Y):
	# convert to index
	idx_lon = int((lon - X[0]) / (X[1]-X[0]) + 0.5)
	idx_lat = int((lat - Y[0]) / (Y[1]-Y[0]) + 0.5)

	# tests on limits
	if idx_lat < 0:
		idx_lat = 0
	if idx_lat >= len(Y):
		idx_lat = len(Y)-1
	if idx_lon < 0:
		idx_lon = len(X)-1
	if idx_lon >= len(X):
		idx_lon = 0
	idx = [idx_lon, idx_lat]
	return idx

###############################################################################

def intersect(src_lon_L, src_lon_R, src_lat_T, src_lat_B,
			  tgt_lon_L, tgt_lon_R, tgt_lat_T, tgt_lat_B):
	# L = left, R = right, T = top, B = bottom
	intersect = True
	if src_lon_R < tgt_lon_L: intersect = False
	if src_lon_L > tgt_lon_R: intersect = False
	# remember that lat is upside down
	if src_lat_B > tgt_lat_T: intersect = False
	if src_lat_T < tgt_lat_B: intersect = False
	return intersect

###############################################################################

def calc_box_weight(src_lon, src_lat, src_DX, src_DY,
					tgt_lon, tgt_lat, tgt_DX, tgt_DY):
	# calculate what percentage of area the box from the source grid 
	# occupies in the target grid

	# half boxes
	src_D_lon = src_DX * 0.5
	src_D_lat = src_DY * 0.5
	tgt_D_lon = tgt_DX * 0.5
	tgt_D_lat = tgt_DY * 0.5

	# coordinates
	src_L1 = src_lon - src_D_lon
	src_P1 = src_lat - src_D_lat
	src_L2 = src_lon + src_D_lon
	src_P2 = src_lat + src_D_lat

	tgt_L1 = tgt_lon - tgt_D_lon
	tgt_P1 = tgt_lat - tgt_D_lat
	tgt_L2 = tgt_lon + tgt_D_lon
	tgt_P2 = tgt_lat + tgt_D_lat

	if intersect(src_L1, src_L2, src_P1, src_P2,
				 tgt_L1, tgt_L2, tgt_P1, tgt_P2):
		# check and modify extents
		if src_L1 < tgt_L1: src_L1 = tgt_L1
		if src_L2 > tgt_L2: src_L2 = tgt_L2
		if src_P1 > tgt_P1: src_P1 = tgt_P1
		if src_P2 < tgt_P2: src_P2 = tgt_P2

		src_area = abs(calculate_grid_area(src_P1, src_L1, src_P2, src_L2))
		tgt_area = abs(calculate_grid_area(tgt_P1, tgt_L1, tgt_P2, tgt_L2))

		# area will be 0.0 at the boxes which span the poles
		if tgt_area == 0.0:
			area_ratio = 1.0
		else:
			area_ratio = src_area / tgt_area
	else:
		area_ratio = 0.0
	return area_ratio

###############################################################################

def calculate_regrid_mapping_and_weights(src_X, src_Y, tgt_X, tgt_Y):
	# loop through the new grid (G2) and check whether grid boxes from the
	# old grid (G1) intersect with the boxes from G2
	
	# get the grid spacings
	tgt_DY = tgt_Y[1]-tgt_Y[0]
	tgt_DX = tgt_X[1]-tgt_X[0]
	src_DY = src_Y[1]-src_Y[0]
	src_DX = src_X[1]-src_X[0]

	# create an array to hold the weights - this is equal to the size of the
	# ratios of the grid spacing multiplied together, times the amount of data
	# stored per weight
	n_weights = abs((2+int(tgt_DX / src_DX)) *
					 (2+int(tgt_DY / src_DY)))
	weights = numpy.zeros([len(tgt_X), len(tgt_Y), 1+n_weights, 3], 'f')

	# get the starting latitude for the tgt grid
	tgt_lat = tgt_Y[0]
	# loop through to calculate weights
	while tgt_lat > tgt_Y[-1]+tgt_DY:
		tgt_lon = tgt_X[0]
		while tgt_lon < tgt_X[-1]+tgt_DX:
			# get the indices of the current latitude and longitude in the tgt grid
			idx_d = get_index_in_grid(tgt_lon, tgt_lat, tgt_X, tgt_Y)
			idx_w = 1	# position in weight array
			# calculate the source grid extents
			idx_e = get_index_in_grid(tgt_lon, tgt_lat, src_X, src_Y)
			# convert this back to the source grid latitude
			src_lat_mid = idx_e[1]*src_DY+src_Y[0]
			src_lat     = src_lat_mid - src_DY
			src_lat_end = src_lat_mid + src_DY

			# loop between these latitudes
			while src_lat >= src_lat_end:
				# calculate the longitude extents for the source grid
				src_lon = (idx_e[0]-1) * src_DX
				src_lon_end = src_lon + 2*src_DX
				# loop between these two extents
				while src_lon <= src_lon_end:
					# quick trap check
					if src_lon-0.5*src_DX > src_lon_end:
						continue
					# calculate the weight of the box
					w = calc_box_weight(src_lon, src_lat, src_DX, src_DY,
										tgt_lon, tgt_lat, tgt_DX, tgt_DY)
					# only add weights that are greater than 0
					if w > 0:
						# get the index in the target grid of the current source location
						idx_s = get_index_in_grid(src_lon, src_lat, src_X, src_Y)
						# assign the weights - these are lon index, lat index, weight
						weights[idx_d[0], idx_d[1], idx_w, 0] = idx_s[0]
						weights[idx_d[0], idx_d[1], idx_w, 1] = idx_s[1]
						weights[idx_d[0], idx_d[1], idx_w, 2] = w
						idx_w += 1
					# next source longitude grid box
					src_lon += src_DX
				# next source latitude grid box
				src_lat += src_DY
			# next target grid box
			tgt_lon += tgt_DX
			# set the number of grid boxes accounted in the box
			weights[idx_d[0], idx_d[1], 0, 0] = idx_w
		# next target grid box
		tgt_lat += tgt_DY
	return weights

###############################################################################

def get_grid_mapping(idx_lon, idx_lat, grid_map):
	n_weights = grid_map[idx_lon, idx_lat, 0, 0]
	if n_weights > 0:
		w = grid_map[idx_lon, idx_lat, 1:, :]
	else:
		w = numpy.zeros([1,3],'f')
	return w

###############################################################################

def cpdn_regrid_field(src_data, tgt_data, grid_map, mv=numpy.inf, accumulation="mean"):
	# regrid a field (X by Y) of data using the computed weights
	for x in range(0, tgt_data.shape[1]):
		for y in range(0, tgt_data.shape[0]):
			# get the indices and weights for this index
			w = get_grid_mapping(x, y, grid_map)
			# sum the weights and data
			sum_d = 0.0
			sum_w = 0.0
			for k in range(0, w.shape[0]):
				# get the data value for the index in the grid_mapping
				data_d = src_data[w[k,1], w[k,0]]
				# do not add missing value
				if abs(data_d) < abs(mv*0.9):		# inaccurate fp representation!
					sum_d += data_d * w[k,2]		# sum the value
					sum_w += w[k,2]					# sum the weight
			if sum_w > 1.0:
				print sum_w
			if sum_w != 0.0:
				# if the accumulation is a mean then divide through by the sum of the weights
				if accumulation == "mean":
					tgt_data[y,x] = sum_d / sum_w		# do the average
				elif accumulation == "sum":		# if accumulation == sum
					tgt_data[y,x] = sum_d
			else:
				tgt_data[y,x] = mv					# missing value
	return tgt_data

###############################################################################

def cpdn_regrid_downscale(src_box, tgt_X, tgt_Y, accumulation="mean"):
	# get the values of the src_box X and Y
	src_X = src_box.get_dimension("X").get_values()
	src_Y = src_box.get_dimension("Y").get_values()

	# get the source axis indices
	src_X_axis = src_box.get_dimension_axes().index("X")
	src_Y_axis = src_box.get_dimension_axes().index("Y")

	# get the missing value attribute
	mv = src_box.get_missing_value()

	# create the output data - size is the same for all dimensions except Y and X
	tgt_shape = list(src_box.get_values().shape)
	# change tgt_shape to reflect new shape
	tgt_shape[src_X_axis] = len(tgt_X)
	tgt_shape[src_Y_axis] = len(tgt_Y)

	# create the target array
	tgt_data = numpy.zeros(tuple(tgt_shape), src_box.get_values().dtype)

	# get the dimensions to iterate over
	it_dims = []
	for d in src_box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

	# calculate the mapping between the old grid and new grid
	grid_map = calculate_regrid_mapping_and_weights(src_X, src_Y, tgt_X, tgt_Y)
	
	# create an iterator
	box_it = cpdn_box_iterator(src_box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)

	# iterate over the box dimensions that are not X and Y
	while not box_it.end():
		src_data = src_box[src_idx].get_values().squeeze()
		cpdn_regrid_field(src_data, tgt_data[tgt_idx], grid_map, mv, accumulation)
		src_idx, tgt_idx = box_it.next(True)

	# create the box
	method_str = "downscale regridding from: X=" + str(len(src_X)) + ", Y=" + str(len(src_Y))
	history_str = datetime.now().isoformat() + " altered by CPDN: downscale regridding"
	tgt_box = create_tgt_box(src_box, tgt_data, tgt_X, tgt_Y, method_str, history_str)

	return tgt_box

###############################################################################

def cpdn_regrid_upscale(src_box, tgt_X, tgt_Y, interp_method="nearest"):
	"""Perform upscaling to higher resolution."""
	# interpret interpolation method
	if interp_method=="linear":
		interp_order = 1
	elif interp_method=="quadratic":
		interp_order = 2
	elif interp_method=="cubic":
		interp_order = 3
	else:
		interp_order = 0

	# get the length of the dimensions first
	src_X = src_box.get_dimension("X").get_values()
	src_Y = src_box.get_dimension("Y").get_values()

	# get the source axis values
	src_X_axis = src_box.get_dimension_axes().index("X")
	src_Y_axis = src_box.get_dimension_axes().index("Y")

	# create the output data - size is the same for all dimensions except Y and X
	tgt_shape = list(src_box.get_values().shape)
	# change tgt_shape to reflect new shape
	tgt_shape[src_X_axis] = len(tgt_X)
	tgt_shape[src_Y_axis] = len(tgt_Y)

	# create the target array
	tgt_data = numpy.zeros(tuple(tgt_shape), src_box.get_values().dtype)

	# create the mesh for the target array - these are "fractional" coordinates relating
	# to the source length and so should be between 0 and the source length
	ratio_X = float(len(src_X)) / len(tgt_X)		# calculate ratio
	ratio_Y = float(len(src_Y)) / len(tgt_Y)
	# construct fractional coordinates
	fractional_tgt_X = [x * ratio_X + 0.5*ratio_X for x in range(0, len(tgt_X))]
	fractional_tgt_Y = [y * ratio_Y + 0.5*ratio_Y for y in range(0, len(tgt_Y))]
	# construct the meshgrid - determine order
	if src_X_axis < src_Y_axis:
		out_grid = numpy.meshgrid(fractional_tgt_Y, fractional_tgt_X)
	else:
		out_grid = numpy.meshgrid(fractional_tgt_X, fractional_tgt_Y)

	# create the output grid
	out_grid = numpy.array(out_grid)
	out_grid = out_grid[[1,0]]		# these have to be swapped
	# get the dimensions to iterate over
	it_dims = []
	for d in src_box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

	# create an iterator
	box_it = cpdn_box_iterator(src_box, it_dims)
	src_idx, tgt_idx = box_it.begin(True)
	mv = src_box.get_missing_value()

	while not box_it.end():
		# get the data here
		src_data = src_box[src_idx].get_values().squeeze()
		# get the min / max value of the field to use in assigning mv later
		mv_in = len(numpy.where(numpy.abs(src_data) > abs(mv*0.9))[0]) != 0
		if mv_in:
			# first fill the missing data
			v0 = numpy.copy(src_data)
			v0 = fill_missing_data_field(v0, mv)
			# regrid the filled data using the interpolation order
			v1 = interp2d(v0, out_grid, order=interp_order, mode='nearest',
			 			  cval=numpy.NAN, prefilter=False)
			# regrid using just nearest neighbour to get missing data values
			v2 = interp2d(src_data, out_grid, order=0, mode='nearest',
			 			  cval=numpy.NAN, prefilter=False)
			# reset mv in v1
			v1[v2==mv] = mv
			v = v1
		else:
			v = interp2d(src_data, out_grid, order=interp_order, mode='nearest', 
						 cval=numpy.NAN, prefilter=False)

		tgt_data[tgt_idx] = v
		# onto the next time / z value
		src_idx, tgt_idx = box_it.next(True)

	# create the box
	method_str = "upscale regridding from: X=" + str(len(src_X)) + ", Y=" + str(len(src_Y))
	history_str = datetime.now().isoformat() + " altered by CPDN: upscale regridding"
	tgt_box = create_tgt_box(src_box, tgt_data, tgt_X, tgt_Y, method_str, history_str)

	return tgt_box

###############################################################################

def cpdn_regrid(src_box, tgt_X, tgt_Y, method="nearest", accumulation="mean"):
	"""Regrid data from a box to the target_X and target_Y coordinate-grid.
		Creates a new cpdn_box from the regridded data.
		Upscaling to a higher resolution grid will interpolate with a spline.
		Downscaling to a lower resolution grid will use area-weighted averaging
		method = interpolation method
		accumulation = sum | mean"""

	# first check if the source box has an X and Y axis
	dim_axes = src_box.get_dimension_axes()
	if not "X" in dim_axes or not "Y" in dim_axes:
		raise Exception("Box has no X or Y dimension, cannot regrid.")

	# get the length of the dimensions first
	src_X = src_box.get_dimension("X").get_values()
	src_Y = src_box.get_dimension("Y").get_values()

	# calculate the ratios between the target and the source
	dim_rat_X = float(len(tgt_X)) / len(src_X)
	dim_rat_Y = float(len(tgt_Y)-1) / (len(src_Y)-1)

	# first check if equivalent
	if dim_rat_X == 1.0 and dim_rat_Y == 1.0:
		# just return a copy of the box
		return copy(src_box)

	# otherwise if > 0 then downscale, < 0 upscale
	if dim_rat_X < 1.0 or dim_rat_Y < 1.0:
		return cpdn_regrid_downscale(src_box, tgt_X, tgt_Y, accumulation)
	elif dim_rat_X > 1.0 or dim_rat_Y > 1.0:
		return cpdn_regrid_upscale(src_box, tgt_X, tgt_Y, method)

###############################################################################

def __fill_missing_data(data, y, x, mv_o, scan_mv=True, scan_mv_v=0.0):
	# fill in the missing data that occurs when regridding with a land-sea mask
	# and the land sea masks do not match
	w_size_x = 2
	w_size_y = 1
	sum = 0.0
	n = 0
	# first try to take an average of a window 5x3
	value_str = ""
	for q in range(-w_size_y, w_size_y+1):
		c_y = y + q
		if c_y < 0: c_y = 0					# don't wrap for y / lat
		if c_y >= data.shape[0]: c_y = data.shape[0]-1
		for p in range(-w_size_x, w_size_x+1):
			c_x = x + p
			if c_x < 0: c_x = data.shape[1]-1	# wrap around for x / lon
			if c_x >= data.shape[1]: c_x = 0
			# average if not original or new missing value
			dc = abs(data[c_y, c_x])
			if dc < abs(mv_o):
				n += 1
				sum += data[c_y, c_x]
				value_str += str(data[c_y, c_x]) + ", "
			else:
				value_str += "_, "
		value_str += "\n"
	if n > 0:
		data[y, x] = sum / n
	elif not scan_mv:
		data[y, x] = scan_mv_v
	else:
		# if no data has been assigned then the box still has missing data in it
		# scan left and right to find values near to the box in the same latitude
		x0 = x
		x1 = x
		# scan left
		while abs(data[y,x0]) >= abs(mv_o)*0.99:
			x0 -= 1
			if x0 < 0: x0 = data.shape[1]-1
			if x0 == x: break
		# scan right
		while abs(data[y,x1]) >= abs(mv_o)*0.99:
			x1 += 1
			if x1 >= data.shape[1]: x1 = 0
			if x1 == x: break
		# interpolate on distance
		x0_d = abs(x0-x)
		x1_d = abs(x1-x)
		td = x0_d + x1_d
		if td != 0:
			data[y,x] = ((td-x0_d) * data[y,x0] + (td-x1_d) * data[y,x1]) / td

	return data

###############################################################################

def fill_missing_data_field(field_data, mv):
	for y in range(0, field_data.shape[0]):
		for x in range(0, field_data.shape[1]):
			__fill_missing_data(field_data, y, x, mv, mv)
	# sometimes top line not filled
	for x in range(0, field_data.shape[1]):
		if abs(field_data[0,x]) > abs(mv*0.9):	#inaccurate fp
			field_data = __fill_missing_data(field_data, 0, x, mv, mv)
	return field_data

###############################################################################

def cpdn_fill_missing_data(src_box):
	# fill the missing data over the box
	pass

###############################################################################

def cpdn_impose_lsm(src_box, lsm, mv_i=2e20, scan_mv=True, scan_mv_v=0.0):
	"""Impose an lsm onto the src_box data.  Any missing data outside the lsm 
		will be assigned a value from the other data"""
	# check that the lsm and data are compatible
	lsm_X_len = lsm.get_dimension("X").get_length()
	lsm_Y_len = lsm.get_dimension("Y").get_length()
	box_X_len = src_box.get_dimension("X").get_length()
	box_Y_len = src_box.get_dimension("Y").get_length()

	if lsm_X_len != box_X_len or lsm_X_len != box_X_len:
		raise Exception("Dimensions of LSM do not match dimensions of src_box")

	# create the target array - copy the data from the source array
	tgt_data = numpy.array(src_box.get_values(), src_box.get_values().dtype)

	# otherwise iterate through the data and impose the lsm
	it_dims = []
	for d in src_box.get_dimensions():
		axis = d.get_axis()
		if axis != "X" and axis != "Y":
			it_dims.append(axis)

	# get the original missing value
	mv_o = src_box.get_missing_value()
	# create the box iterator
	box_it = cpdn_box_iterator(src_box, it_dims)
	src_idx = box_it.begin()
	# get the indices of the land sea mask that are 1
	lsm_idx = numpy.where(lsm.get_values().squeeze() != 0)
	# reuse the iterator
	while not box_it.end():
		# get the location of the data where the value is mv_o and interpolate over it
		c_data = tgt_data[src_idx]
		mv_idx = numpy.where(numpy.abs(c_data) > abs(mv_o*0.9))	# inaccurate fp
		for i in range(0, mv_idx[0].squeeze().shape[0]):
			x = mv_idx[1][i]
			y = mv_idx[0][i]
			# fill the missing data
			__fill_missing_data(c_data, y, x, mv_o, scan_mv, scan_mv_v)
		# now impose the land sea mask
		c_data[lsm_idx] = mv_i
		src_idx = box_it.next()

	# create the target box
	method_str = "imposition of Land Sea Mask"
	history_str = datetime.now().isoformat() + " altered by CPDN: imposed LSM"
	tgt_box = create_tgt_box(src_box, tgt_data, 
							 src_box.get_dimension("X").get_values(), 
							 src_box.get_dimension("Y").get_values(), 
							 method_str, history_str, mv_i)

	return tgt_box
