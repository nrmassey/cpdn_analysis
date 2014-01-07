###############################################################################
# File         : cpdn_helpers.py
# Author       : Neil Massey
# Created      : 01/08/12
# Purpose      : miscellaneous helper functions for cpdn analysis framework
# Changes      : 
###############################################################################

from cpdn_box import *
import numpy
from copy import *

###############################################################################

def cpdn_clone_box(o_box, new_data=None, T_vals=None, T_bnds=None, 
				   method_string=None, history_string=None, calendar=None,
				   name=None):
	# clone a box but with a new T axis, data and additions to the attributes
	# commonly done when analysing data - use the original box to produce a
	# 2nd box which is then altered to reflect what the analysis has done.
    # create tgt_dims
	tgt_dims = []
	for d in o_box.get_dimensions():
		if d.get_axis() == "T":
			# get the start date, units and number of days in the year
			sd, un, nd = d.get_time_details()
			# create the values and the bounds
			if T_vals == None:
				vals = copy(d.get_values())
			else:
				vals = T_vals
			if T_bnds == None:
				bnds = copy(d.get_bounds())
			else:
				bnds = T_bnds
			# get the attributes - manipulate the calendar if required
			tgt_attrs = d.get_attributes()
			if calendar != None:
				tgt_attrs["calendar"] = calendar

			# create the time dimension and append
			td = cpdn_boxdim(name=d.get_name(), vals=vals, attrs=tgt_attrs,
                             axis="T", bounds=bnds, start_date=sd, time_units=un,
                             n_days_per_year=nd)
			tgt_dims.append(td)
		else:
			tgt_dims.append(d)

	# create the return box
	tgt_attrs = amend_attributes(o_box.get_attributes(), method_string, history_string)
	tgt_glob_attrs = o_box.get_global_attributes()
	# overwrite tgt_off and tgt_sf as get_values scales by these
	tgt_off = 0.0
	tgt_sf = 1.0
	if name != None:
		oname = name
	else:
		oname = o_box.get_name()
	if new_data == None:
		new_data = o_box.get_values()
	# copy rotated grid if exists
	if o_box.has_rotated_grid():
		rot_grid = o_box.get_rotated_grid()
	else:
		rot_grid = None
		
	clone_box = cpdn_box(dims=tgt_dims, var_attrs=tgt_attrs, glob_attrs=tgt_glob_attrs,
         	             name=oname, off=tgt_off, sf=tgt_sf, data=new_data,
         	             rotated_grid=rot_grid)
	return clone_box

