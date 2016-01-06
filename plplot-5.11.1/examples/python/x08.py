#!/usr/bin/env python
#
#	x08c.c
#
#	3-d plot demo.

from numpy import *
import math
import sys
import os

module_dir = "@MODULE_DIR@"

if module_dir[0] == '@':
	module_dir = os.getcwd ()

sys.path.insert (0, module_dir)

XPTS = 35		# Data points in x
YPTS = 46		# Data points in y

opt = [1, 2, 3, 3]

alt = [60.0, 20.0, 60.0, 60.0]

az = [30.0, 60.0, 120.0, 160.0]

title = ["#frPLplot Example 8 - Alt=60, Az=30, Opt=1",
	 "#frPLplot Example 8 - Alt=20, Az=60, Opt=2",
	 "#frPLplot Example 8 - Alt=60, Az=120, Opt=3",
	 "#frPLplot Example 8 - Alt=60, Az=160, Opt=3"]

# main
#
# Does a series of 3-d plots for a given data set, with different
# viewing options in each plot.

def main(w):

##    # Parse and process command line arguments
##
##    pl.ParseOpts(sys.argv, pl.PARSE_FULL)
##
##    # Initialize plplot
##
##    pl.init()

##    x = []
##    y = []
##    z = []

    x = zeros( XPTS, 'd' ); y = zeros( YPTS, 'd' )
    #z = zeros( (XPTS, YPTS), 'd' )
    z = reshape( zeros( XPTS*YPTS, 'd' ), (XPTS, YPTS) )

    for i in range(XPTS):
	x[i] = float(i - (XPTS / 2)) / float(XPTS / 2)

    for i in range(YPTS):
	y[i] = float(i - (YPTS / 2)) / float(YPTS / 2)

    for i in range(XPTS):
	xx = x[i]
##	zz = []
	for j in range(YPTS):
	    yy = y[j]
	    r = math.sqrt(xx * xx + yy * yy)
##	    zz.append(math.exp(-r * r) *
##		      math.cos(2.0 * math.pi * r))
	    z[i,j] = math.exp(-r * r) * math.cos(2.0 * math.pi * r)
##    z.append(zz)

    for k in range(4):
	w.pladv(0)
	w.plvpor(0.0, 1.0, 0.0, 0.9)
	w.plwind(-1.0, 1.0, -0.9, 1.1)
	w.plcol0(1)
	w.plw3d(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0,
		alt[k], az[k])
	w.plbox3("bnstu", "x axis", 0.0, 0,
		 "bnstu", "y axis", 0.0, 0,
		 "bcdmnstuv", "z axis", 0.0, 0)

	w.plcol0(2)
	w.plot3d(x, y, z, opt[k], 1)
	w.plcol0(3)
	w.plmtex("t", 1.0, 0.5, 0.5, title[k])

	w.pleop()
##    pl.end()
##
##main()
