#  Copyright (C) 2007, 2008 Andrew Ross

#  Grid data demo.
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
import time

pts = 500
xp = 25
yp = 20
nl = 16
knn_order = 20
threshold = 1.001
wmin = -1e3
randn = 0
rosen = 0


def cmap1_init():
    i = [0, 1]
    h = [240, 0]
    l = [0.6, 0.6]
    s = [0.8, 0.8]

    plscmap1n(256)
    plscmap1l(0, i, h, l, s)

    

    
# main
#
#
def main():

    title = ["Cubic Spline Approximation",
		   "Delaunay Linear Interpolation",
		   "Natural Neighbors Interpolation",
		   "KNN Inv. Distance Weighted",
		   "3NN Linear Interpolation",
		   "4NN Around Inv. Dist. Weighted"]

    opt = [0., 0., 0., 0., 0., 0.]

    xm = -0.2
    ym = -0.2
    xM = 0.6
    yM = 0.6

    opt[2] = wmin
    opt[3] = knn_order
    opt[4] = threshold
 
    cmap1_init()
    plseed(5489)

    # Create the sampled data
    # For consistency across languages use plrandd to create the 
    # pseudo-random data that are required.
    xt = zeros(pts)
    yt = zeros(pts)
    for i in range(pts) :
	xt[i] = (xM-xm)*plrandd()
	yt[i] = (yM-ym)*plrandd()
    if randn == 0 :
        x = xt + xm
        y = yt + ym
    else:
        x = sqrt(-2.*log(xt))*cos(2.*pi*yt) + xm
        y = sqrt(-2.*log(xt))*sin(2.*pi*yt) + ym

    if rosen == 0:
        r = sqrt( x*x + y*y)
        z = exp(-r*r)*cos(2*pi*r)
    else:
        z = log(pow(1.-x,2) + 100.*pow(y-pow(x,2),2))

    zmin = min(z)
    zmax = max(z)

    # Create regular grid
    xg = xm + (xM-xm)*arange(xp)/(xp-1.)
    yg = ym + (yM-ym)*arange(yp)/(yp-1.)

    xx = zeros(1)
    yy = zeros(1)

    plcol0(1)
    plenv(xm,xM,ym,yM,2,0)
    plcol0(15)
    pllab('X','Y','The original data sampling')
    for i in range(pts):
        plcol1( ( z[i] - zmin )/( zmax - zmin ) )
        xx[0] = x[i]
        yy[0] = y[i]
        plstring( xx, yy, '#(727)' )
    pladv(0)

    plssub(3,2)

    for k in range(2):
        pladv(0)
        for alg in range(1,7):
            zg = plgriddata(x, y, z, xg, yg, alg, opt[alg-1])

            if alg == GRID_CSA or alg == GRID_DTLI or alg == GRID_NNLI or alg == GRID_NNI:
                for i in range(xp):
                    for j in range(yp):
                        # Average (IDW) over the 8 neighbours for Nan
                        if isnan(zg[i][j]):
                            zg[i][j] = 0.0
                            dist = 0.0
                            for ii in range(max(i-1,0),min(i+2,xp)):
                                for jj in range(max(j-1,0),min(j+2,yp)):
                                    if (not isnan(zg[ii][jj])):
                                        d = abs(ii-i) + abs(jj-j)
                                        if (d != 1.0) :
                                            d = 1.4142
                                        zg[i][j] += zg[ii][jj]/(d*d)
                                        dist += d
                                        
                            if dist != 0.0 :
                                zg[i][j] /= dist
                            else:
                                zg[i][j] = zmin
            
            lzM = max(zg.flat)
            lzm = min(zg.flat)

	    lzm = min(lzm,zmin)
            lzM = max(lzM,zmax)

            lzm = lzm - 0.01
            lzM = lzM + 0.01
            
            plcol0(1)
            pladv(alg)

            if (k == 0):
                clev = lzm + (lzM-lzm)*arange(nl)/(nl-1)
                plenv0(xm,xM,ym,yM,2,0)
                plcol0(15)
                pllab('X','Y',title[alg-1])
                plshades(zg, xm, xM, ym, yM, clev, 1., 1, None, None)
                plcol0(2)
            else:
                clev = lzm + (lzM-lzm)*arange(nl)/(nl-1)
                plvpor(0.0, 1.0, 0.0, 0.9)
                plwind(-1.1, 0.75, -0.65, 1.20)
                plw3d(1.0, 1.0, 1.0, xm, xM, ym, yM, lzm, lzM, 30.0, -40.0)
                plbox3('bntu', 'X', 0.0, 0,'bntu', 'Y', 0.0, 0,'bcdfntu', 'Z', 0.5, 0)
                plcol0(15)
                pllab('','',title[alg-1])
                plot3dc(xg, yg, zg, DRAW_LINEXY | MAG_COLOR | BASE_CONT, clev)
    
main()
