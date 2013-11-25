import math

########################################################################################
### Convert global coords to rotated (regional) coords

def glob2rot(lon, lat, pole_lon, pole_lat):

	# Make sure rotlon is between 0 and 360
	while (lon  >= 360.0):
		lon  -= 360.0
	while (lon  <    0.0):
		lon  += 360.0

    # Make sure pole_lon is between 0 and 360
	while (pole_lon >= 360.0):
		pole_lon -= 360.0
	while (pole_lon <    0.0):
		pole_lon += 360.0

	# Convert inputs to radians
	lon_r = math.radians(lon)
	lat_r = math.radians(lat)
	pole_lon_r = math.radians(pole_lon)
	pole_lat_r = math.radians(pole_lat)

	# Amount rotated about 180E meridian
	if (pole_lon_r == 0.0):
		sock = 0.0;
	else:
		sock = pole_lon_r - math.pi

	# Need to get the screw in range -pi to pi
	screw = lon_r - sock;
	while (screw < -1.0 * math.pi):
		screw += 2.0 * math.pi
	while (screw > math.pi):
		screw -= 2.0 * math.pi
	bpart = math.cos(screw) * math.cos(lat_r)

	x = (-1.0 * math.cos(pole_lat_r) * bpart) + (math.sin(pole_lat_r) * math.sin(lat_r))
	if (x >=  1.0):
		x =  1.0
	if (x <= -1.0):
		x = -1.0
	lat2 = math.asin(x)

	t1 = math.cos(pole_lat_r) * math.sin(lat_r)
	t2 = math.sin(pole_lat_r) * bpart

	x = (t1 + t2) / math.cos(lat2)

	if (x >=  1.0):
		x =  1.0
	if (x <= -1.0):
		x = -1.0
	lon2 = -1.0 * math.acos(x)

	if (screw >= 0.0 and screw <= math.pi):
		lon2 *= -1.0

	# Convert back to degrees
	lon = math.degrees(lon2)
	lat = math.degrees(lat2)

	return lon, lat

########################################################################################
### Convert rotated (regional) to global grid

def rot2glob(lon, lat, pole_lon, pole_lat):
	# Make sure rotlon is between 0 and 360
	while (lon >= 360.0):
		lon -= 360.0
	while (lon <    0.0):
		lon += 360.0
	# Make sure pole_lon is between 0 and 360
	while (pole_lon >= 360.0):
		pole_lon -= 360.0
	while (pole_lon < 0.0): 
		pole_lon += 360.0

	# Convert inputs to radians
	lon = math.radians(lon)
	lat = math.radians(lat)
	pole_lon = math.radians(pole_lon)
	pole_lat = math.radians(pole_lat)

	# Amount rotated about 180E meridian
	if (pole_lon == 0.0):
		sock = 0.0
	else:
		sock = pole_lon - math.pi

	cpart = math.cos(lon) * math.cos(lat)
	x = math.cos(pole_lat) * cpart + math.sin(pole_lat) * math.sin(lat)

	if (x >=  1.0):
		x =  1.0
	if (x <= -1.0):
		x = -1.0
	lat_out = math.asin(x)

	t1 = -1.0 * math.cos(pole_lat) * math.sin(lat)
	t2 = math.sin(pole_lat) * cpart

	x = (t1 + t2) / math.cos(lat_out)
	if (x >=  1.0):
		x =  1.0
	if (x <= -1.0):
		x = -1.0
	lon_out = -1.0 * math.acos(x)
 
    # Make sure rotlon is between -PI and PI    
	while (lon < -1*math.pi):
		lon += 2.0*math.pi
	while (lon > math.pi):
		lon -= 2.0*math.pi

	if (lon >= 0.0 and lon <= math.pi):
		lon_out *= -1.0
	lon_out += sock;
	# Convert back to degrees
	lon_out = math.degrees(lon_out)
	lat_out = math.degrees(lat_out)

	return lon_out, lat_out

########################################################################################

rot2reg = rot2glob
reg2rot = glob2rot

if __name__ == "__main__":
	pole_lat = 39.25
	pole_lon = 198
	lon = 6.8
	lat = 38.5
	lon_r, lat_r = glob2rot(lon, lat, pole_lon, pole_lat)
	lon_g, lat_g = rot2glob(lon_r, lat_r, pole_lon, pole_lat)
	print lon, lat
	print lon_r, lat_r
	print lon_g, lat_g
