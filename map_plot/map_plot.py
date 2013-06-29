############################################################################
# File         : map_plot.py
# Author       : Neil Massey
# Created      : 09/03/12
# Purpose      : General class to plot climate model fields on various 
#                projections, in various manners (block, contours, wind, etc.)
#				 LATitude coordinates have the range(90, -90) = (90N, 90S)
#				 LONgitude coordinates have the range(-180, 180) = (180W, 180E)
# Changes	   : 19/06/12 - adapted to work with cpdn_box
############################################################################

import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import os
import projections
import pickle
import numpy
import matplotlib.cm as cm
from copy import *
from types import *

############################################################################

class map_plot:

	def __init__(self, projection="cylindrical", grid=True):
		# currently enabled projections: cylinder, mercator, robinson, polar
		# record the projection
		self.projection = projection
		# set up the figure
		self.fig = plt.figure()
		self.sp = self.fig.add_subplot("211", projection=projection)
		# set the position of the plots
		self.cbh = 0.025
		self.cbb = 0.05
		self.sp.set_position([self.cbb, self.cbh+2*self.cbb, 
							  1.0-2*self.cbb, 1.0-self.cbh-3*self.cbb])
		self.sp.grid(grid)
		self.continents = []
		self.block_stage = "none"		# record what has and hasn't been
		self.cb_stage = "none"			# rendered yet
		self.z_max = -2e20				# max and min of data
		self.z_min = 2e20				# for the colour bar drawing
		self.x_limit = [-180,180]
		self.y_limit = [-90, 90]

	########################################################################

	def set_grid_spacing(self, lon_spacing, lat_spacing):
		"""
			Set the number of degrees between each longitude / latitude grid.
		"""
		self.sp.set_longitude_grid(lon_spacing)
		self.sp.set_latitude_grid(lat_spacing)

	########################################################################

	def prepare_data(self, LON, LAT, data, invert=False):
		# mirror LAT if necessary - imshow needs extents to match array order
		n_data = data
		if self.projection == "cylindrical" and invert:
			if LAT[0] > LAT[-1]:
				LAT = LAT[::-1]
			else:
				data = data[::-1,:]

		# change the LONgitude if necessary, so that it runs between -180 and 180
		lon_over_180 = numpy.where(LON > 180)
		if lon_over_180[0].shape[0] != 0:
			LON[lon_over_180] -= 360
			s = lon_over_180[0][0]
			e = lon_over_180[0][-1] + 1
			p = data.shape[1] - s
			n_data = numpy.zeros([data.shape[0], data.shape[1]], data.dtype)
			n_data[:,0:p] = data[:,s:e]
			n_data[:,p:]  = data[:,0:s]
			LON.sort()
		else:
			n_data = data
		return LON, LAT, n_data

	########################################################################

	def __get_rotated_lon_lat_corners(self, LON, LAT, io, jo):
		# get the four corners of the latitude / longitude box at position i,j
		if io == LON.shape[0]-1:
			i = LON.shape[0]-2
		elif io == 0:
			i = 1
		else:
			i = io
			
		if jo == LON.shape[1]-1:
			j = LON.shape[1]-2
		elif jo == 0:
			j = 1
		else:
			j = jo
		# box origin
		lon_o = LON[io,jo]
		lat_o = LAT[io,jo]
		
		# upper left
		lon_ul = lon_o - 0.5*(lon_o - LON[i-1,j-1])
		lat_ul = lat_o - 0.5*(lat_o - LAT[i-1,j-1])
		# upper right
		lon_ur = lon_o + 0.5*(LON[i-1,j+1] - lon_o)
		lat_ur = lat_o + 0.5*(LAT[i-1,j+1] - lat_o)
		# lower right
		lon_lr = lon_o + 0.5*(LON[i+1,j+1] - lon_o)
		lat_lr = lat_o + 0.5*(LAT[i+1,j+1] - lat_o)
		# lower left
		lon_ll = lon_o - 0.5*(lon_o - LON[i+1,j-1])
		lat_ll = lat_o - 0.5*(lat_o - LAT[i+1,j-1])
		
		x = [lon_ul, lon_ur, lon_lr, lon_ll]
		y = [lat_ul, lat_ur, lat_lr, lat_ll]
		return x,y

	########################################################################

	def block_plot(self, box, colmap=cm.RdBu_r, interpolate='nearest', 
				   z_max=-2e20, z_min=2e20, mv=None, sym_col=False, mvc="#FFFFFF"):
		"""plot the field as blocks using the colormap referenced by cmap"""
		# adapted to work with a cpdn_box
		# have to copy LON and LAT so they can be manipulated
		# check first for rotated grid
		if box.has_rotated_grid():
			rg = box.get_rotated_grid()
			LON = rg.get_global_longitude_values()
			LAT = rg.get_global_latitude_values()
		else:
			LON = copy(box.get_dimension("X").get_values())
			LAT = copy(box.get_dimension("Y").get_values())
		
		# get the data values from the box
		Z = numpy.array(box.get_values().squeeze())
		# check that Z has the correct shape
		if len(Z.shape) != 2:
			raise Exception("Z has dimensions other than 2 - cannot plot")
		# check Z dimensions match LON and LAT
		
		if box.has_rotated_grid():
			if Z.shape[0] != LON.shape[0] or Z.shape[1] != LON.shape[1]:
				raise Exception("Z has dimensions different to LON - cannot plot")
			if Z.shape[0] != LAT.shape[0] or Z.shape[1] != LAT.shape[1]:
				raise Exception("Z has dimensions different to LAT - cannot plot")				
		else:
			if Z.shape[1] != LON.shape[0]:
				raise Exception("Z has dimensions different to LON - cannot plot")
			if Z.shape[0] != LAT.shape[0]:
				raise Exception("Z has dimensions different to LAT - cannot plot")

		# if mv set then create masked array
		if isinstance(mv, NoneType):
			if box.get_missing_value() != numpy.inf:
				mv = box.get_missing_value()
				Zm = numpy.ma.masked_where(Z==mv, Z)
			else:
				Zm = Z
		else:
			Zm = numpy.ma.masked_where(Z==mv, Z)

		if z_max == -2e20:
			self.z_max = numpy.max(Zm)
		else:
			self.z_max = z_max

		if z_min == 2e20:
			self.z_min = numpy.min(Zm)
		else:
			self.z_min = z_min

		# should we have symmetrical colorbars?
		if sym_col:
			if abs(self.z_min) > abs(self.z_max):
				self.z_max = abs(self.z_min)
			else:
				self.z_min = -abs(self.z_max)

		# get the location of missing data
		# normalize the data
		z_nrm = (Zm - self.z_min) / (self.z_max-self.z_min)		
		self.colmap = colmap
		# prepare the data in case the LAT etc. are upside down
		# get the region extents
		if not box.has_rotated_grid:
			LON, LAT, z_nrm = self.prepare_data(LON, LAT, z_nrm, invert=True)
			region_extent = [LON[0], LON[-1], LAT[0], LAT[-1]]
		else:
			region_extent = [numpy.min(LON), numpy.min(LAT), numpy.max(LON), numpy.max(LAT)]

		# can use imshow for the cylindrical projection - it is much faster as well
		if box.has_rotated_grid:
			# special routine for rotated grid
			for i in range(0, LON.shape[0]):
				for j in range(0, LAT.shape[1]):
					# get the four corners of the box - this will be a quadrilateral - not a rectangle
					# origin
					x,y = self.__get_rotated_lon_lat_corners(LON, LAT, i, j)
					# get the colour
					V = z_nrm[i,j]
					if V != numpy.inf:
						C = colmap.__call__(int(V*colmap.N))
					else:
						C = mvc	# missing value colour
					self.sp.fill(x, y, fc=C, ec=C, zorder=1)
		elif self.projection == "cylindrical":
			colmap.set_over(mvc,1)
			colmap.set_under(mvc,1)
			# imshow is fast but only works with cylindrical projection and non-rotated grids :(
			self.sp.imshow(z_nrm, cmap=colmap, interpolation=interpolate, extent=region_extent, 
						   clip_on=False, vmin=0.0, vmax=1.0, zorder=0)
		else:
			# get the half width of the grid boxes
			x_h = abs(LON[1] - LON[0]) / 2.0
			y_h = abs(LAT[1] - LAT[0]) / 2.0
			# have to do a nested loop plot
			for i in range(0, len(LON)):
				for j in range(0, len(LAT)):
					# get the four corners of the box
					x = [LON[i]-x_h, LON[i]+x_h, LON[i]+x_h, LON[i]-x_h]
					y = [LAT[j]-y_h, LAT[j]-y_h, LAT[j]+y_h, LAT[j]+y_h]
					# use the color from the colormap
					V = z_nrm[i,j]
					if V != numpy.inf:
						C = colmap.__call__(int(V*colmap.N))
					else:
						C = mvc
					self.sp.fill(x, y, fc=C, ec=C, zorder=1)
		# update the state machines
		self.block_stage = "drawn"
		if self.cb_stage == "to_draw":
			self.draw_colorbar(self.cb_title, self.precision, self.ptype)
		
		# set the region - maximum / minimum lats / lons for rotated grid
		self.set_region(region_extent)

	########################################################################

	def get_format_string(self, precision, ptype):
		if ptype=="sig_figs":
			format = "%.0"+str(precision)+"g"
		elif ptype=="integer":
			format = "%i"
		else:
			format = "%.0"+str(precision)+"f"
		return format

	########################################################################

	def draw_colorbar(self, title="none", precision=5, ptype="sig_figs", n_ticks=11,
					  V=[], v_max=-2e20, v_min=2e20):
		if self.block_stage == "none":
			# don't draw yet, unless the block plot has been drawn
			self.cb_stage = "to_draw"
			self.cb_title = title
			self.precision = precision
			self.ptype = ptype
		else:
			# block plot has been drawn so add the colorbar
			self.cb = self.fig.add_subplot("212")
			self.cb.set_position([self.cbb, self.cbb*1.5, 1.0-2*self.cbb, self.cbh])
			# set no y ticks or labels and don't draw the x tick lines
			self.cb.xaxis.set_ticks_position('none')
			self.cb.yaxis.set_ticks_position('none')
			self.cb.yaxis.set_ticks([])
			# set the x tick values and labels
			if V == []:
				xt = numpy.linspace(self.z_min, self.z_max, n_ticks)
				self.cb.xaxis.set_ticks(xt)
				format = self.get_format_string(precision, ptype)
				xtl = [format % float(xt[i]) for i in range(0, len(xt))]   # limit to two d.p.
			else:
				self.cb.xaxis.set_ticks(V)
				format = self.get_format_string(precision, ptype)
				xtl = [format % float(V[i]) for i in range(0, len(V))]   # limit to two d.p.
			self.cb.xaxis.set_ticklabels(xtl)
			# finally! draw the colorbar
			sc = float(self.z_max - self.z_min) / self.colmap.N
			for c in range(0, self.colmap.N):
				x = self.z_min + c * sc
				col = self.colmap.__call__(c)
				self.cb.fill([x,x,x+sc,x+sc],[0,1,1,0],fc=col,ec=col)
			if V == []:
				self.cb.set_xlim([self.z_min, self.z_max])
			else:
				self.cb.set_xlim([V[0], V[-1]])
			if v_max != -2e20 and v_min != 2e20:
				self.cb.set_xlim(v_min, v_max)

			if title != "none":
				self.cb.set_title(title)
				self.cb_title = title
			self.cb_stage = "drawn"

	########################################################################

	def contour_plot(self, box, nc=10, z_max=-2e20, z_min=2e20,
					 color='k', labels=True, precision=5, ptype="sig_figs",
					 contours=[], neg_ls="dashed", pos_ls="solid", lw=1.0):
		"""plot the field as a contour"""
		# adapted to work with a cpdn_box
		# have to copy LON and LAT so they can be manipulated
		LON = copy(box.get_dimension("X").get_values())
		LAT = copy(box.get_dimension("Y").get_values())
		Z = numpy.array(box.get_values().squeeze())
		# translate the LON, LAT and Z to the coordinate system
		LON, LAT, Z = self.prepare_data(LON, LAT, Z)
		# create the contour values to draw
		if contours != []:
			V = contours
		elif z_max != -2e20 and z_min != 2e20:
			V = numpy.linspace(z_min, z_max, nc)
		else:
			V = nc
		# get the line styles for the contours
		ls = []
		for v in V:
			if v < 0.0:
				ls.append(neg_ls)
			else:
				ls.append(pos_ls)
		# draw the contours - all in the same color
		CS = self.sp.contour(LON, LAT, Z, V, colors=color, linewidths=lw, linestyles=ls)
		# do the labels
		format = self.get_format_string(precision, ptype)
		self.sp.clabel(CS, fontsize=8, inline=True, inline_spacing=-2, 
					   fmt=format)

	########################################################################

	def vector_plot(self, box_u, box_v):
		"""plot the U and V components of the vector"""
		pass

	########################################################################

	def path_plot(self, path):
		"""plot a path from a list of (lon, lat) pairs"""
		pass

	########################################################################

	def set_xlim(self, xlim):
		self.x_limit = [xlim[0], xlim[1]]
		self.sp.set_xlim(self.x_limit)

	########################################################################

	def set_ylim(self, ylim):
		self.y_limit = [ylim[0], ylim[1]]
		self.sp.set_ylim(self.y_limit)

	########################################################################

	def set_region(self, region):
		# region is [xmin, ymin, xmax, ymax]
		self.x_limit = [region[0], region[2]]
		self.y_limit = [region[1], region[3]]
		self.sp.set_xlim(self.x_limit)
		self.sp.set_ylim(self.y_limit)

	########################################################################

	def set_title(self, title):
		self.sp.set_title(title)

	########################################################################

	def draw_continents(self, fillcolor='none', cn=0, lw=0.5, alpha=1.0):
		"""draw the continents on the map"""
		res=1
		if fillcolor == 'none':
			ec = 'k'
		else:
			ec = fillcolor
		if self.continents == []:
			self.read_continents()
		cont = self.continents[cn]
		for cont in self.continents:
			cont_p = patches.PathPatch(cont, edgecolor=ec, facecolor=fillcolor, lw=lw, zorder=1,
									   alpha=alpha)
			self.sp.add_patch(cont_p)

	########################################################################

	def save(self, filename, dpi=300):
		"""save the map as a file - the filename extension of filename implies
			the filetype it will be saved as, e.g. .eps, .png, .jpg etc."""
		self.fig.savefig(filename, dpi=dpi, bbox_inches="tight", pad_inches=0.25)

	########################################################################

	def gcp(self):
		"""Get the map as a matplotlib subplot"""
		return self.sp

	########################################################################

	def gcf(self):
		return self.fig

	########################################################################

	def read_continents(self):
		# check first whether the continents have already been pickled
		path = "/Users/massey/Coding/cpdn_analysis/map_plot/data/"
		pik_name = path + "map_continents.pik"
		if os.path.exists(pik_name):
			fh = open(pik_name)
			self.continents = pickle.load(fh)
			fh.close()
		else:
			# read in the continents file
			fh = open(path + "map_continents_lowres.txt")
			lines = fh.readlines()			
			fh.close()
			# get the info line by line - each continent takes up 3 lines
			continent_points = []
			c_cont = 0
			for i in range(0, len(lines), 3):
				if (lines[i][0] != "#"):
					print "error in file"
				lons = lines[i+1].split(",")
				lats = lines[i+2].split(",")
				# build a list of the continents lat & lon coords
				cont = [[float(lons[0]), float(lats[0])]]
				phase = 0
				for j in range(1, len(lons)):
					pt = [float(lons[j]), float(lats[j])]
					# check for wrap around date line
					if abs(pt[0] - cont[-1][0]) < 5.0:
						cont.append(pt)
					else:
						# create new continent - special case for europe / africa and antartica
						if c_cont == 0:
							# draw a box to the new poitn
							cont.append([-180,-90])
							cont.append([180,-90])
							cont.append(pt)
						elif c_cont == 186:
							phase += 1
							if phase == 1:
								continent_points.append(cont)
								cont = [pt]
							elif phase == 2:
								continent_points.append(cont)
								cont = continent_points[-2]
								cont.append(pt)

				c_cont += 1
				continent_points.append(cont)
			
			# create the path commands for the continent
			for cont in continent_points:
				cont_codes = [Path.MOVETO]
				for pt_n in range(1, len(cont)):
					cont_codes.append(Path.LINETO)
						
				# create a path from cont and cont_path
				cont_path = Path(cont, cont_codes)				
				self.continents.append(cont_path)
			# write the data out as a pickle so we only have to construct it once
			fo = open(pik_name, "w")
			pickle.dump(self.continents, fo)
			fo.close()
