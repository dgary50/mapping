#+
# Project     : SMM-XRP
#
# Name        : PLOT_HELIO
#
# Purpose     : Plot solar heliographic grid and limb
#
# Category    : Plotting
#
# Explanation : Uses heliographic formulae from:
#               "Textbook on Spherical Astronomy",
#               by W.M. Smart (see page 174).
#
# Syntax      : plot_helio,date,grid=grid,glabel=glabel
#
# Inputs      : DATE = date/time
#
# Keywords    : GRID_SPACING = spacing (deg) of lat-long  [default = 10]
#               GLABEL = label lat-long grid with coordinate values (default = no labels)
#               OVER = control erasing of previously drawn plot. Setting this
#                        to 1 will cause grid to be overlaid on a drawn plot
#               GCOLOR  = color index for grid
#               GFONT   = gfont index for grid label
#               ROLL = solar roll (deg) clockwise from solar north
#               RCENTER = roll center [default = 0,0]
#               GSTYLE = grid linestyle [default=0]
#               GTHICK = grid thickness
#               XRANGE,YRANGE = data ranges in arcsecs
#               LIMB_PLOT =  plot limb
#               LCOLOR = color index for limb
#               LTHICK = thickness of limb
#               NO_ROLL_CORRECT = do not roll grid
#
# History     : 2020-May-26  DG
#                 Ported from Dominic Zarro's IDL code of the same name.
#
# Contact     : dgary@njit.edu
#-

from __future__ import print_function
from astropy.time import Time
import numpy as np
from mapping.map_util import *
from mapping.transform_solar import hel2arcmin
dtor = np.pi/180.0   # Global conversion of degrees to radians

def plot_helio(date=None, roll=False, rcenter=[0.,0.], 
               grid = 10., glabel=False, gcolor='k', gfont=None, gstyle='--', gthick=1, 
               over=False, xrange=[-1300,1300], yrange=[-1300,1300], 
               limb=False, lcolor='k', lthick=1, lstyle='-',
               tgrid=None, no_roll_correct=False, l0=0., axis=None):

    import matplotlib.pylab as plt
    from copy import deepcopy
    grid = abs(grid)

    #-- date

    try:
        tdate = Time(date).iso
    except ValueError:
        print('PLOT_HELIO: Error, could not interpret given date as ISO time.')
        return

    #-- need solar radius for grid

    if limb or roll:
        out = hel2arcmin(0, 0, date=tdate)#, arcsec=True)
        radius = out['angles']['rsun']*60.

    #-- check for non-zero roll

    roll_center = [0.,0.]
    do_roll = False

    if roll:
        do_roll = (roll % 360.) != 0.
        on_disk = (rcenter[0]**2 + rcenter[1]**2) < (radius)**2
        if do_roll and on_disk:
            roll_center = rcenter

    if no_roll_correct:
        do_roll = False

    if limb:
        ang = np.linspace(0,2*np.pi,361)
        xlimb = radius*np.cos(ang)
        ylimb = radius*np.sin(ang)
        if do_roll:
            xlimb, ylimb = roll_xy(xlimb,ylimb,roll,rcenter=roll_center)

    #-- define latitude-longitude matrices (in degrees) at desired grid spacing

    do_grid = grid > 0.

    xcor = ycor = None
    if do_grid:
        gv = -90. + np.arange(0,181.,grid)
        ng = len(gv)
        if tgrid is None:
            tgrid = 0.5
        v = -90. + np.arange(0,180.1,tgrid)
        npts = len(v)
        lon = v.repeat(npts).reshape(npts,npts)
        lat = deepcopy(np.transpose(lon))
        lon = lon + l0

        #-- compute cartesian coordinates of grid points
        #lat = np.transpose(lat)
        #lon = np.transpose(lon)
        out = hel2arcmin(lat, lon, date=tdate, l0=l0)#, arcsec=True)
        xcor = out['xx'].reshape(npts,npts)*60.
        ycor = out['yy'].reshape(npts,npts)*60.
        if do_roll:
             xcor, ycor = roll_xy(xcor, ycor, roll, rcenter=roll_center)

    if not over:
        f, axis = plt.subplots(1,1)
        axis.set_aspect('equal','box')
        axis.set(xlim=xrange, ylim=yrange)
    if axis is None:
        print('PLOT_HELIO: Error, no axis provided for overplot.')
        return

    if limb:
        axis.plot(xlimb, ylimb, color=lcolor, linewidth=lthick, linestyle=lstyle)

    #-- plot latitude-longitude lines

    if do_grid:
        for i in range(npts):
            ok, = np.where(v[i] == gv)
            if len(ok) > 0:
                axis.plot(xcor[i], ycor[i], linestyle=gstyle, color=gcolor, linewidth=gthick)
        for j in range(npts):
            ok, = np.where(v[j] == gv)
            if len(ok) > 0:
                axis.plot(xcor[:,j], ycor[:,j], linestyle=gstyle, color=gcolor, linewidth=gthick)

    #-- label grid coordinates within current viewport

        if glabel:
            glon = gv.repeat(ng).reshape(ng,ng)
            glat = np.transpose(glon)
            out = hel2arcmin(glat, glon, date=tdate)#, arcsec=True)
            xcor = out['xx'].reshape(ng,ng)*60.
            ycor = out['yy'].reshape(ng,ng)*60.
            if do_roll:
                xcor, ycor = roll_xy(xcor, ycor, roll, rcenter=roll_center)
            x0, x1, y0, y1 = axis.axis()
            within, = np.where( np.logical_and(np.logical_and(np.logical_and(np.logical_and(np.logical_and(abs(glon) < 90.,abs(glat) < 90.), x0 <= xcor),xcor <= x1), y0 <= ycor), ycor <= y1))
            if len(within) != 0:
                for idx in within:
                    wlon = glon[idx]
                    wlat = glat[idx]
                    clon = str(np.rint(wlon)).strip()
                    clat = str(np.rint(wlat)).strip()
                    axis.text(xcor[idx], ycor[idx], clon+','+clat) #,font=gfont
    return
