import matplotlib.pyplot as plt
import numpy
from os.path import expanduser
EU = expanduser
from scipy.io.netcdf import *
from map_plot import *

# load some CPDN data in
fname = EU("~/ORWELL/HadAM3P_Output/1959/hadam3p_m24k_1959_2_006104590_2/m24kma.paf9dec.nc")
nc_fh = netcdf_file(fname, 'r')
lon = nc_fh.variables["longitude0"][:]
lat = nc_fh.variables["latitude0"][:]
data = nc_fh.variables["field16"][:]
nc_fh.close()

mp = map_plot(projection="robinson", grid=True)
fig = mp.get()

#mp.draw_continents(fillcolor='blue')
mp.draw_colorbar("Temperature in K", precision=0, ptype="dec_pla")
mp.block_plot(lon, lat, data[0,0].squeeze())
mp.contour_plot(lon, lat, data[0,0].squeeze(), nc=6, ptype="integer")

#mp.set_xlim([-90,90])
#mp.set_ylim([-45,45])

mp.save("robin.png")
