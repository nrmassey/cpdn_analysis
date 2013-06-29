###############################################################################
# File         : cpdn_region.py
# Author       : Neil Massey
# Created      : 18/05/12
# Purpose      : routine to subset a cpdn_box based on a given region
# Changes      : 
###############################################################################
# regions are: name: [min_lat, min_lon, max_lat, max_lon], abbreviation
###############################################################################

from cpdn_box import *

region_list = { "Australia"                 : ([-45, 110, -10, 155], "Aus"),
                "Amazon Basin"              : ([-20, 280,  10, 325], "Amz"),
                "Southern South America"    : ([-55, 285, -20, 320], "SSA"),
                "Central America"           : ([ 10, 245,  30, 275], "CAm"),
                "Western North America"     : ([ 30, 230,  60, 255], "WNA"),
                "Central North America"     : ([ 30, 255,  50, 275], "CNA"),
                "Eastern North America"     : ([ 25, 275,  50, 300], "ENA"),
                "Alaska"                    : ([ 60, 190,  70, 255], "AK" ),
                "Greenland"                 : ([ 50, 255,  85, 350], "Grn"),
                "Mediterranean Basin"       : ([ 30, 350,  50,  40], "Med"),
                "Northern Europe"           : ([ 50, 350,  75,  40], "NEu"),
                "Western Africa"            : ([-10, 340,  20,  20], "WAf"),
                "Eastern Africa"            : ([-10,  20,  20,  50], "EAf"),
                "Southern Africa"           : ([-35, 350, -10,  50], "SAf"),
                "Sahara"                    : ([ 20, 340,  30,  65], "Sah"),
                "Southeast Asia"            : ([-10,  95,  20, 155], "SEA"),
                "East Asia"                 : ([ 20, 100,  50, 145], "EAi"),
                "South Asia"                : ([  5,  65,  30, 100], "SAi"),
                "Central Asia"              : ([ 30,  40,  50,  75], "CAi"),
                "Tibet"                     : ([ 30,  75,  50, 100], "Tib"),
                "Northeast Asia"            : ([ 30,  40,  70, 180], "NEA"),
                "Antarctica"                : ([-90,   0, -65, 360], "Ant"),
                "Northern Hemisphere"       : ([  0,   0,  90, 360], "NH" ),
                "Southern Hemisphere"       : ([-90,   0,   0, 360], "SH" ),
                "NH Extratropics"           : ([ 30,   0,  90, 360], "NHE"),
                "SH Extratropics"           : ([-90,   0, -30, 360], "SHE"),
                "Global"                    : ([-90,   0,  90, 360], "Glo"),
                "United Kingdom"            : ([51.25, 350.625, 57.5, 1.875], "UK"),
                "United Kingdom 2"          : ([58.75,  -11.25, 50.0, 1.875], "UK2"),
                "North Atlantic"            : ([ 20, 285,  60, 360], "NA"),
        # HadAM3P Regions
                "HA3P Global"               : ([-90,   0,  90, 358.125], "HA3P_Glo"),
                "HA3P Northern Hemisphere"  : ([  0,   0,  90, 358.125], "HA3P_NH"),
                "HA3P Southern Hemisphere"  : ([-90,   0,   0, 358.125], "HA3P_SH"),
                "HA3P Eur Atl"              : ([36.25, 315, 67.5, 380.625], "HA3P_Eur_Atl"),
                "HA3P pa Eur Atl"           : ([36.25, 315, 67.5, 20.625], "HA3P_pa_Eur_Atl"),
                "HA3P Eur Atl 2 "           : ([67.5, -45.0, 36.25, 20.625], "HA3P_pa_Eur_Atl2"),
                "HA3P NW USA"               : ([ 30, 228.75, 55, 255], "HA3P_NW_USA"),
                "HA3P South Africa"         : ([-37.5, 11.25, -16.25, 37.5], "HA3P_SA"),
                "HA3P South Africa 2"       : ([-16.25, 11.25, -37.5, 37.5], "HA3P_SA2"),
                "HA3P India China"          : ([ 10, 69.375, 50, 131.25], "HA3P_IndChi"),
                "HA3P pa UK"                : ([51.25, 350.625, 58.0, 1.875], "HA3P_pa_UK"),
                "HA3P pb UK"                : ([51.25, 350.625, 58.0, 361.875], "HA3P_pb_UK"),

        # HadISST oceans
                "HadISST North Atlantic"    : ([   0, 280,  54,   3], "HSST_NA"),
                "HadISST South Atlantic"    : ([ -60, 290,   0,  20], "HSST_SA"),
                "HadISST North Pacific"     : ([   0, 130,  54, 260], "HSST_NP"),
                "HadISST South Pacific"     : ([ -60, 130,   0, 290], "HSST_SP"),
                "HadISST Arctic"            : ([  54,   0,  90, 359], "HSST_Art"),
                "HadISST Antarctic"         : ([ -90,   0, -60, 359], "HSST_Ant"),
                "HadISST Indian"            : ([ -60,  20,  15, 130], "HSST_Ind"),
                "HadISST Mediterranean"     : ([  32,   3,  45,  32], "HSST_Med"),
                }

###############################################################################

def get_region_from_box(box, region_name):
	"""Get a new box from a cpdn_box subset to the region.
	   returns : new box which is a subset of the original box over the region"""

	# first look for the region in the keys
	if region_name in region_list.keys():
		region = region_list[region_name][0]
	else:
		# otherwise look in the abbreviations
		found = False
		for k in region_list.keys():
			if region_list[k][1] == region_name:
				found = True
				region = region_list[k][0]
				break
		if found == False:
			raise Exception("Region name " + region_name + " not defined in region_list")

	# build the indexing list
	idx = []
	box_dims = box.get_dimensions()
	for d in range(0, len(box_dims)):
		if box_dims[d].get_axis() == "X":
			idx.append(slice(float(region[1]),float(region[3]),None))

		if box_dims[d].get_axis() == "Y":
			idx.append(slice(float(region[0]),float(region[2]),None))

		if box_dims[d].get_axis() == "T":
			idx.append(slice(None, None, None))

		if box_dims[d].get_axis() == "Z":
			idx.append(slice(None, None, None))
	return box.__getitem__(idx)
