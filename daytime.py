###############################################################################
# File         : daytime.py
# Author       : Neil Massey
# Created      : 21/05/12
# Purpose      : functions to convert between daytime and a floating points
#                representation of date, using either a 360 or 365 day calendar
# Changes      : 
###############################################################################

n_minutes_per_day = 24*60
n_secs_per_day = 24*60*60

# number of days elapsed in the year before the month
days_elapsed = [0,31,59,90,120,151,181,212,243,273,304,334,365]
#				  31 28 31  30  31  30  31  31  30  31  30  31
#				  J  F  M   A   M   J   J   A   S   O   N   D

from datetime import *

###############################################################################

class daytime(object):
	"""Class to replace datetime, for dates which might have 30 days in Feb,
		due to 360 day calendar"""

	###########################################################################

	def __init__(self, years, months, days, hours=0, minutes=0, seconds=0):
		self.year   = years
		self.month  = months
		self.day    = days
		self.hour   = hours
		self.minute = minutes
		self.second = seconds

	###########################################################################

	def __eq__(self, other):
		# convert to floating point representation and compare
		return daytime_to_float(self) == daytime_to_float(other)

	###########################################################################

	def __lt__(self, other):
		# convert to floating point representation and compare
		return daytime_to_float(self) < daytime_to_float(other)

	###########################################################################

	def __le__(self, other):
		# convert to floating point representation and compare
		return daytime_to_float(self) <= daytime_to_float(other)

	###########################################################################

	def __gt__(self, other):
		# convert to floating point representation and compare
		return daytime_to_float(self) > daytime_to_float(other)

	###########################################################################

	def __ge__(self, other):
		# convert to floating point representation and compare
		return daytime_to_float(self) >= daytime_to_float(other)

	###########################################################################

	def isoformat(self, separator="T"):
		"""Return a string in isoformat date"""
		tstring = "%02i-%02i-%02i" % (self.year, self.month, self.day)
		tstring += separator
		tstring += "%02i:%02i:%5f" % (self.hour, self.minute, self.second)
		return tstring

###############################################################################

def daytime_to_float(dt, n_days_per_year=365.25):
	"""Convert a daytime to the number of days since 0, given an arbitrary 
	   number of days in a year"""

	c_time = float(dt.year-1) * int(n_days_per_year)
	if n_days_per_year == 365.25:
		# calculate number of leap years
		n_leaps = 0
		for y in range(0, dt.year):
			if y % 400 == 0:
				n_leaps += 1
			elif y % 100 == 0:
				n_leaps += 0
			elif y % 4 == 0:
				n_leaps += 1
		c_time += n_leaps
		c_time += days_elapsed[dt.month-1] + 1
		# check for leap year
		if dt.year % 4 == 0 and dt.month == 2:
			c_time += 1.0
	else:
		c_time += (dt.month-1) * n_days_per_year/12
	c_time += (dt.day-1)
	c_time += float(dt.hour) / 24
	c_time += float(dt.minute) / n_minutes_per_day
	c_time += float(dt.second) / n_secs_per_day
	return c_time

###############################################################################

def float_to_daytime(fl, n_days_per_year=365.25):
	if n_days_per_year == 365.25:
		fl_to_datetime = date.fromordinal(int(fl)-1)
		n_years = fl_to_datetime.year
		n_month = fl_to_datetime.month
		n_days  = fl_to_datetime.day
		days_left = fl - int(fl)
	elif n_days_per_year == 365.0:
		n_years = int(fl / n_days_per_year)
		days_left = fl - int(n_years * n_days_per_year)
		n_month =  int(days_left / (n_days_per_year / 12))
		days_left = days_left - days_elapsed[n_month]
		n_days = int(days_left)
		days_left = days_left - n_days
	else:
		n_years = int(fl / n_days_per_year)
		days_left = fl - int(n_years * n_days_per_year)
		n_month =  int(days_left / (n_days_per_year / 12))
		days_left = days_left - n_month * (n_days_per_year / 12)
		n_days = int(days_left)
		days_left = days_left - n_days
	# calculate the remaining hours / mins / secs from the fractional part
	n_hours = int(1e-4+days_left * 24)					# + 1e-4 to remove rounding errors
	days_left = days_left - float(n_hours)/24
	n_minutes = int(1e-4+days_left * n_minutes_per_day)
	days_left = days_left - float(n_minutes) / n_minutes_per_day
	n_seconds = int(1e-4+days_left*n_secs_per_day)
	# convert to daytime
	if n_days_per_year == 365.25:
		dt = daytime(n_years, n_month, n_days, n_hours, n_minutes, n_seconds)
	else:
		dt = daytime(n_years+1, n_month+1, n_days+1, n_hours, n_minutes, n_seconds)
	return dt

###############################################################################

def datestring_to_daytime(date_string):
	# strip-time is such a pain to use - I've rewritten this
	try:
		date, time = date_string.split(" ")
	except:
		date = date_string      # cannot split!
		time = None
	date_split = date.split("-")
	if time != None:
		time_split = time.split(":")
		# make sure there are enough time splits to cope with indexing below
		while len(time_split) < 3:
			time_split.append("0")
		# fun and games!
		if "." in time_split[2]:
			time_split[2] = str(int(float(time_split[2])))
		do = daytime(int(date_split[0]), int(date_split[1]), int(date_split[2]),
                     int(time_split[0]), int(time_split[1]), int(time_split[2]))
	else:
		do = daytime(int(date_split[0]), int(date_split[1]), int(date_split[2]),
                     0,0,0)
	return do

###############################################################################

if __name__ == "__main__":
	dt = daytime(1, 12, 1, 0, 0, 0)
	# test both 365 & 365.25 day calendars
	x = daytime_to_float(dt, 360)
	print float_to_daytime(x, 360)
	x += 390.5
	print float_to_daytime(x, 360)

	x = daytime_to_float(dt, 365.25)
	print float_to_daytime(x, 365.25)
	x += 365
	print float_to_daytime(x, 365.25)

