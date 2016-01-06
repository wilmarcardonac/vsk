#  Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 Alan W. Irwin

#  Multiple window and color map 0 demo.
#
#  This file is part of PLplot.
#
#  PLplot is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Library General Public License as published
#  by the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  PLplot is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Library General Public License for more details.
#
#  You should have received a copy of the GNU Library General Public License
#  along with PLplot; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

from plplot_py_demos import *

# Demonstrates multiple windows and default color map 0 palette.

# Draws a set of numbered boxes with colors according to cmap0 entry
def draw_windows(nw, cmap0_offset):
    plschr(0.0, 3.5)
    plfont(4)

    for i in range(nw):
	plcol0(i+cmap0_offset)
	pladv(0)
	vmin = 0.1
	vmax = 0.9
	for j in range(3):
	    plwidth(j + 1)
	    plvpor(vmin, vmax, vmin, vmax)
	    plwind(0.0, 1.0, 0.0, 1.0)
	    plbox("bc", 0.0, 0, "bc", 0.0, 0)
	    vmin = vmin + 0.1
	    vmax = vmax - 0.1
	plwidth(1)
	plptex(0.5, 0.5, 1.0, 0.0, 0.5, `i`)

# Demonstrate multiple windows and default color map 0 palette.
def demo1():

    plbop()
    # Divide screen into 16 regions
    plssub(4, 4)

    draw_windows(16,0)

    pleop()


def demo2():
    lmin = 0.15
    lmax = 0.85

    plbop()

    # Divide screen into 100 regions
    plssub(10,10)
    
    r = zeros(116,"int")
    g = zeros(116,"int")
    b = zeros(116,"int")

    for i in range(100):
        h = (360./10.) * (i%10)
        l = lmin + (lmax-lmin)*(i/10)/9
        s = 1.0

        rgb = plhlsrgb(h, l, s)

        r[i+16] = rgb[0]*255.001
        g[i+16] = rgb[1]*255.001
        b[i+16] = rgb[2]*255.001

    # Load default cmap0 colors into out custom set
    for i in range(16):
        rgb = plgcol0(i)
        r[i] = rgb[0]
        g[i] = rgb[1]
        b[i] = rgb[2]

    # Now set cmap0 all at once (faster, since fewer driver calls)
    plscmap0(r,g,b)

    draw_windows(100,16)

    pleop()

def main():

    # For starting from scratch this call to pladv increments cursub, but 
    # then the following plssub sets it to zero so the whole thing is 
    # essentially a nop.  However, for the case when other examples are run 
    # first, this call to pladv is absolutely essential to finish the 
    # preceding page.
    pladv(0)

    # Run demos
    demo1()
    demo2()

    # Restore defaults
    plssub(1, 1)
    plfont(1)
#    plcol0(1)
    pleop()

main()
