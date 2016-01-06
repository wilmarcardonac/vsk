# xw08.py PLplot demo for Python
#
# Copyright (C) 2004-2014  Alan W. Irwin
#
# This file is part of PLplot.
#
# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA

from plplot_py_demos import *

# These values must be odd, for the middle
# of the index range to be an integer, and thus
# to correspond to the exact floating point centre
# of the sombrero.

XPTS = 35		# Data points in x
YPTS = 45		# Data points in y

alt = [60.0, 40.0]

az = [30.0, -30.0]

title = ["#frPLplot Example 8 - Alt=60, Az=30",
	 "#frPLplot Example 8 - Alt=40, Az=-30"]

# Routine for restoring colour map1 to default.
# See static void plcmap1_def(void) in plctrl.c for reference.
def restore_cmap1():
    # For center control points, pick black or white, whichever is closer to bg
    # Be careful to pick just short of top or bottom else hue info is lost
    vertex = sum(array(plgcolbg()))/(3.*255.)
    if vertex < 0.5:
	vertex = 0.01
	midpt = 0.10
    else:
	vertex = 0.99
	midpt = 0.90
    # Independent variable of control points.
    i = array((0., 0.44, 0.50, 0.50, 0.56, 1.))
    # Hue for control points.  Blue-violet to red
    h = array((260., 260., 260., 0., 0., 0.))
    # Lightness ranging from medium to vertex to medium
    l = array((0.5, midpt, vertex, vertex, midpt, 0.5))
    # Saturation is complete for default
    s = array((1., 1., 1., 1., 1., 1.))
    # Default number of cmap1 colours
    plscmap1n(128)
    # Interpolate between control points to set up default cmap1.
    plscmap1l(0, i, h, l, s)

# Routine for defining a specific color map 1 in HLS space.
# if gray is true, use basic grayscale variation from half-dark to light.
# otherwise use false color variation from blue (240 deg) to red (360 deg).
def cmap1_init(gray):
    # Independent variable of control points.
    i = array((0., 1.))
    if gray:
	# Hue for control points.  Doesn't matter since saturation is zero.
	h = array((0., 0.))
	# Lightness ranging from half-dark (for interest) to light.
	l = array((0.5, 1.))
	# Gray scale has zero saturation
	s = array((0., 0.))
    else:
	# Hue ranges from blue (240 deg) to red (0 or 360 deg)
	h = array((240., 0.))
	# Lightness and saturation are constant (values taken from C example).
	l = array((0.6, 0.6))
	s = array((0.8, 0.8))

    # number of cmap1 colours is 256 in this case.
    plscmap1n(256)
    # Interpolate between control points to set up cmap1.
    plscmap1l(0, i, h, l, s)
# main
#
# Does a series of 3-d plots for a given data set, with different
# viewing options in each plot.

def main():

    rosen = 0
    dx = 2. / float( XPTS - 1 )
    dy = 2. / float( YPTS - 1 )

    x = -1. + dx*arange(XPTS)
    y = -1. + dy*arange(YPTS)
    if rosen == 1:
	x = 1.5*x
	y = 0.5 + y
    x.shape = (-1,1)
    r2 = (x*x) + (y*y)
    if rosen == 1:
	z = (1. - x)*(1. - x) + 100 * (x*x - y)*(x*x - y)
	# The log argument might be zero for just the right grid.
	z = log(choose(greater(z,0.), (exp(-5.), z)))
    else:
	z = exp(-r2)*cos((2.0*pi)*sqrt(r2))

    x.shape = (-1,)
    zmin = min(z.flat)
    zmax = max(z.flat)
    nlevel = 10
    step = (zmax-zmin)/(nlevel+1)
    clevel = zmin + step + arange(nlevel)*step

    # Set up data and arrays for plsurf3dl call below.
    indexxmin = 0
    indexxmax = XPTS
    # Must be same shape as z, and a row of z.
    zlimited = empty(z.shape)
    indexymin = empty(z.shape[0], dtype=int)
    indexymax = empty(z.shape[0], dtype=int)
    # Parameters of ellipse that limits the data.
    x0 = 0.5*(XPTS - 1)
    a = 0.9*x0
    y0 = 0.5*(YPTS - 1)
    b = 0.7*y0
    for i in range(indexxmin, indexxmax):
         square_root = sqrt(1. - min(1., ((double(i) - x0)/a)**2))
         # Add 0.5 to find nearest integer and therefore preserve symmetry
         # with regard to lower and upper bound of y range.
         indexymin[i] = max(0, int(0.5 + y0 - b*square_root))
         # indexymax calculated with the convention that it is 1
         # greater than highest valid index.
         indexymax[i] = min(YPTS, 1 + int(0.5 + y0 + b*square_root))
         zlimited[i][indexymin[i]:indexymax[i]] = z[i][indexymin[i]:indexymax[i]]

    pllightsource(1., 1., 1.)

    for k in range(2):
	for ifshade in range(5):
	    pladv(0)
	    plvpor(0.0, 1.0, 0.0, 0.9)
	    plwind(-1.0, 1.0, -0.9, 1.1)
	    plcol0(3)
	    plmtex("t", 1.0, 0.5, 0.5, title[k])
	    plcol0(1)
	    if rosen == 1:
		plw3d(1.0, 1.0, 1.0, -1.5, 1.5, -0.5, 1.5, zmin, zmax,
		alt[k], az[k])
	    else:
		plw3d(1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0, zmin, zmax,
		alt[k], az[k])
	    plbox3("bnstu", "x axis", 0.0, 0,
	    "bnstu", "y axis", 0.0, 0,
	    "bcdmnstuv", "z axis", 0.0, 0)

	    plcol0(2)
	    if ifshade == 0:
		# diffuse light surface plot.
		# set up modified gray scale cmap1.
		cmap1_init(1)
		plsurf3d(x, y, z, 0, ())
	    elif ifshade == 1:
		# magnitude colored plot.
		cmap1_init(0)
		plsurf3d(x, y, z, MAG_COLOR, ())
	    elif ifshade == 2:
		# magnitude colored plot with faceted squares
		cmap1_init(0)
		plsurf3d(x, y, z, MAG_COLOR | FACETED, ())
	    elif ifshade == 3:
		# magnitude colored plot with contours
		cmap1_init(0)
		plsurf3d(x, y, z, MAG_COLOR | SURF_CONT | BASE_CONT, clevel)
	    elif ifshade == 4:
		# magnitude colored plot with contoursmagnitude colored plot and index limits
		cmap1_init(0)
		plsurf3dl(x, y, zlimited, MAG_COLOR | SURF_CONT | BASE_CONT, clevel, indexxmin, indexymin, indexymax)

    # Restore defaults
    #plcol0(1)
    restore_cmap1()

main()
