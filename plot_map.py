#+
# Project     : SOHO_CDS
#
# Name        : PLOT_MAP
#
# Purpose     : Plot an image map
#
# Category    : imaging
#
# Syntax      : plot_map,map
#
# Inputs      : MAP = image structure map created by MAKE_MAP
#
# Keywords    :
#     /OVERLAY = overlay on previous image
#     /CONT = contour the image
#     SMOOTH_WIDTH = smoothing width (> 1)
#     FOV = [fx,fy] = field of view to be plotted
#     GRID_SPACING = grid spacing (deg) for latitude-longitude grid [default= 0, no grid]
#     GLABEL = label grid with coordinate values [default = 0, no labels]
#     GSTYLE = grid linestyle [default =0]
#     CENTER = [xc,yc] = center coordinates (arcsec) of FOV [default = center of image]
#            (if center is a valid map, then use its center)
#     DMIN,DMAX = min, max data for plot [default = data min,max]
#     BORDER = draw border around image [default = no]
#     /DEBUG = turn on extra debugging
#     /TAIL = allows user to tailor contours
#     /LOG_SCALE  = log_10 scale image
#     WINDOW = window index to send plot to
#     /NOAXES = inhibit plotting axes
#     /NODATA = inhibit plotting data (don't plot anything)
#     /NO_DATA = inhibit plotting data (but include axes, grid, etc) 
#     /NOTITLE = inhibit printing title
#     /NOLABELS = inhibit axis labels
#     /NOXTICKS = inhibit X-tick labels
#     /NOYTICKS = inhibit Y-tick labels
#     /DROTATE  = solar rotate image contour
#     LEVELS  = user specified contour levels
#     NLEVELS = # of default levels [default=10]
#     /PLUS_ONLY = plot positive data
#     /MINUS_ONLY = plot negative data
#     XRANGE,YRANGE = cartesian plot limits
#     /INTERLACE = interlace two images when overlaying
#     /COMPOSITE = simultaneously plot two images when overlaying
#                = type of compositing:
#                  1: original , 2: latest(new), 3: max(new/old) 4: min(new/old)
#     /AVERAGE   = average two images when using /COMPOSITE
#     BOTTOM = lowermost color index to use in color scaling [def=0]
#     LAST_SCALE = use MIN/MAX from previous plot
#     LIMB_PLOT = overplot solar limb
#     BTHICK = border thickness
#     BCOLOR = border line color (defaults to white)
#     LCOLOR = limb line color
#     LTHICK = limb line thickness
#     MULTI = set for multiple plots per page, e.g. mult=[2,2] (or 2) for 4
#             plots per page (!p.multi remains at this value until cleared)
#     NOERASE = don't erase previous plot
#     SQUARE_SCALE = force equal aspect ratio (default)
#     ERR_MSG = String error message if any
#     STATUS = 0/1 means failure/success
#     CBAR = 0/1 means draw colorbar on image plots (only works in > IDL 5.2)
#     PERCENT = if levels are entered, they are in % of data max
#     MARK_POINT = if set to a 2-element array, it is the x,y data coords of a point to be marked
#        If point is not within plot, and arrow off edge of plot shows direction to point.
#     DURTITLE = If set, plot title will include duration of image.  Title will be
#        map.id +  dd-mmm-yy hh:mm:ss-hh:mm:ss
#     NO_BYTE_SCALE = set to not byte scale images
#     ALPHA = transparency factor between 0 and 1 to blend to images.
#     XSHIFT, YSHIFT = translation shifts (arcecs) to input map (positive in
#     X and Y). 
#     ROLL_ANGLE = angle (degrees) to rotate input map (positive clockwise).
#     NO_FILL = don't fill area outside of plot box with background.
#     RED, GREEN, BLUE = (R,G,B) color arrays
#     TRUE_COLORS = display using true colors
#     BLEND = blend overlayed image with base image
#     ABOUT_CENTER = rotate map about_center when overlaying
#
# Restrictions:
#      - do not set /OVERLAY unless a plot exists on the current device
#
# History     : 2020-May-26  DG
#                 Ported from Dominic Zarro's IDL code of the same name.
#
# Contact     : dgary@njit.edu
#-

from __future__ import print_function
import matplotlib.pylab as plt
import numpy as np
from mapping.map_util import *
from mapping.transform_map import rot_map, shift_map, drot_map
from mapping.transform_solar import pb0r

def last_wp():
    import pickle
    import os
    
    last = None
    # Implement common block functionality by use of a pickle file "pmap_shared.pkl"
    common_file = 'pmap_shared.pkl'
    if os.path.exists(common_file):
        # Read file if it exists
        fp = open(common_file, "rb")
        last = pickle.load(fp)
        fp.close()
    
    if last:
        window = last['window']
        panel = last['panel'] 
    else:
        window = None
        panel = None
    return window, panel

def wp_ok(window, panel):
    ''' Quick test to see if requested window, panel combination is available
    '''
    try:
        if window in plt.get_fignums():
            if panel > 0 and panel <= len(plt.figure(window).get_axes()):
                return True
    except:
        print('WP_OK: Syntax truth = wp_ok(window, panel)')
    return False

def get_axes(window=None):

    panel_axes = []
    cbar_axes = []
    if window is None:
        window, panel = last_wp()
    if window in plt.get_fignums():
        fig = plt.figure(window)
        axes = fig.get_axes()
        # Remove any axes from list that are not "Subplot" axes
        for ax in axes:
            if str(ax).find('Subplot') == -1:
                cbar_axes.append(ax)
            else:
                panel_axes.append(ax)
    return {'pax':panel_axes, 'cax':cbar_axes}

def get_cbar(window=None, panel=None):
    w, p = last_wp()
    if window is None:
        window = w
    if panel is None:
        panel = p

    if wp_ok(window, panel):
        out = get_axes(window)
        ax = out['pax'][panel-1]
        pleft, pbot, pwid, pht = ax.get_position().bounds
        ptop = pbot + pht
        pright = pleft + pwid
        for cax in out['cax']:
            cleft, cbot, cwid, cht = cax.get_position().bounds
            if cwid > cht:
                # This is a horizontal cbar, so compare with top and bottom
                if abs(cbot - ptop) < 0.01 and abs(pleft - cleft) < 0.01:
                    return cax
            if cwid < cht:
                # This is a vertical cbar, so compare with right and left
                if abs(cbot - pbot) < 0.01 and abs(pright - cleft) < 0.01:
                    return cax
    return None
            
    

def clr(window=None, panel=None, item='all'):
    ''' Helper routine to clear contours from the current window, or panel specified
        by window, panel combination.
    '''
    
    if wp_ok(window, panel):
        fig = plt.figure(window)
        axes = fig.get_axes()
        axis = axes[panel-1]
        if item == 'all':
            axis.cla()
        elif item == 'images':
            img_list = axis.images
            while(len(img_list) > 0):
                img_list[-1].remove()
        elif item == 'lines':
            line_list = axis.lines
            while(len(line_list) > 0):
                line_list[-1].remove()
        elif item == 'contours':
            coll_list = axis.collections
            while(len(coll_list) > 0):
                coll_list[-1].remove()
        elif item == 'cbar':
            cax = get_cbar(window, panel)
            cax.remove()
        else:
            print('CLR: Item must be images, lines, contours, cbar, or all')

def plot_map(map=None, cont=False, over=False, panel=None,
                   about_center=False, bcolor='white', border=True, bthick=1,
                   cbar=False, center=None, color=None, composite=False, date_only=False,
                   dmax=None, dmin=None, drange=None, drotate=False, duration=None,
                   durtitle=False, fov=None, grid=0., last_scale=False, levels=None,
                   limb=False, log=False, newwindow=False, #mark_point=mark_point,
                   minus_only=False, multi=[1,1], new=False, nlevels=10, no_data=False,
                   noaxes=False, nodata=False, nodate=False, noerase=False, saved_map=None,
                   nolabels=False, notitle=False, original_time=False, percent=False,
                   plus_only=False, roll_angle=0., smooth=False, title=None, 
                   window=None, xrange=None, xshift=0., yrange=None, yshift=0., 
                   **kwargs):

    #-- some variables saved in memory for overlay

    import pickle
    import os
    from copy import deepcopy
    from scipy.ndimage import gaussian_filter
    from time import sleep
    
    # Implement common block functionality by use of a pickle file "pmap_shared.pkl"
    common_file = 'pmap_shared.pkl'
    # Set default values if file does not exist
    last = {'window':None, 'panel':None, 'time':None, 'drange':None, 'scale': False,
            'xrange':None, 'yrange':None, 'multi':None, 'roll':0., 'rcenter':[0.,0.],
            'b0':None, 'l0':0., 'rsun':960., 'saved_map':None, 'roll_correct':False,
            'center':[0.,0.]}
    if os.path.exists(common_file):
        # Read file if it exists
        fp = open(common_file, "rb")
        last.update(pickle.load(fp))
        fp.close()

    shifting = False
    rolling = False

    # Strip any grid and limb keywords from kwargs, and set to value or default
    glabel = default(kwargs, 'glabel', False)
    gcolor = default(kwargs, 'gcolor', 'w')
    gstyle = default(kwargs, 'gstyle', ':')
    gthick = default(kwargs, 'gthick', 1)
    gfont  = default(kwargs, 'gfont', None)
    lcolor = default(kwargs, 'lcolor', 'w')
    lstyle = default(kwargs, 'lstyle', '-')
    lthick = default(kwargs, 'lthick', 1)
    
    #-- overlay limb and/or grid on previous plot

    if over and map is None:
        from mapping.plot_helio import plot_helio
        if last['time'] is None:
            print('PLOT_MAP: Error: No previous image on which to overlay limb/grid')
            return False
        if about_center:
            rcenter = last['center'] 
        else:
            rcenter = last['rcenter']
        if xrange is None:
            xrange = last['xrange']
        if yrange is None:
            yrange = last['yrange']
        axes = np.array(plt.figure(last['window']).get_axes()).flatten()
        if len(axes) <= last['panel']:
            axis = axes[last['panel']-1]
        else:
            axis = None
        plot_helio(last['time'], roll=last['roll'], axis=axis, grid=grid, 
                   over=True, rcenter=rcenter, l0=last['l0'], 
                   glabel=glabel, gcolor=gcolor, gfont=gfont, gstyle=gstyle, gthick=gthick, 
                   limb=limb, lcolor=lcolor, lthick=lthick, lstyle=lstyle,
                   xrange=xrange, yrange=yrange,
                   tgrid=0.5, no_roll_correct=last['roll_correct'])
        return

    #-- check input map

    if map is None:
        print('PLOT_MAP: Syntax plot_map,map')
        return
    if not valid_map(map):
        print('PLOT_MAP: Error, the provided argument is not a valid map dict.')
        return

    #-- check image scalings

    if log and 'log_scale' in map.keys():
        if map['log_scale']:
            print('PLOT_MAP: Input map already log-scaled.')
            log = False

    #-- always overlay as a contour unless composite = True

    if not cont and not composite:
        cont = over
    if not noerase: noerase = over

    rolling = (roll_angle % 360.) != 0.
    shifting = xshift != 0. or yshift != 0.

    #-- open a new figure if one doesn't exist
    #-- else get axis for previous plot

    nowindow = False
    nopanel = False
    if window is None:
        # No window specified in call, so default to the last-used window
        window = last['window']
    if not (window in plt.get_fignums()):
        # The requested window is not open
        nowindow = True
        nopanel = True
        last['panel'] = None
    if not nowindow:
        # There is an open window, so look for a viable axis
        if panel is None:
            # No panel was specified in call, so default to the last-used panel
            panel = last['panel']
        # Compare with available axes in the window
        fig = plt.figure(window)
        axes = fig.get_axes()
        # Remove any axes from list that are not "Subplot" axes
        good_axes = []
        for ax in axes:
            if str(ax).find('Subplot') != -1:
                good_axes.append(ax)
        axes = good_axes
        # See if the current axis is in the window
        naxes = len(axes)
        if panel <= naxes:
            axis = axes[panel-1]
            # Current axis was found, so use it for an overplot, or increment
            # to next available axis if not
            if not over:
                if panel == naxes:
                    # Need a new axis, but the current axis is already the last one.
                    nopanel = True
                    if newwindow:
                        nowindow = True
                    else:
                        panel = 1
                        axis = axes[panel-1]
                        axis.cla()
                else:
                   # Set axis to the next available axis.
                   axis = axes[panel]
                   panel += 1
        else:
            panel = None
            nopanel = True
            
    if over:
        if nowindow:
            print('PLOT_MAP: Overlay base window unavailable')
            return
        if nopanel:
            print('PLOT_MAP: Overlay panel unavailable')
            return
    else:
        if nowindow:
            # Open a new window for plotting
            fig, axes = plt.subplots(multi[1],multi[0])
            last['multi'] = multi
            window = fig.number
            axes = fig.get_axes()        
            panel = 1
            axis = axes[panel-1]
                
    # fig, window and panel now point to the correct figure, window, and panel for the next plot.
    last['window'] = window
    last['panel'] = panel
    
    #-- translating or rolling map

    saved_icenter = [map['xc'],map['yc']]
    saved_rcenter = get_map_prop(map,'roll_center',default=saved_icenter)

    if shifting or rolling:
        saved_map = deepcopy(map)
        if rolling: map = rot_map(map, roll_angle)
        if shifting: map = shift_map(map, xshift, yshift)

    odmin = float(np.min(map['data']))
    odmax = float(np.nanmax(map['data']))
    if odmin == 0 and odmax == 0:
        print('PLOT_MAP: Error, All data are zero')
        return

    #-- filter NaN's

    off_scale = odmax*100.
    if map['data'].dtype == np.byte:
        pic = np.float(map['data'])
    else:
        pic = deepcopy(map['data'])
    inan = np.where(np.isnan(pic))
    if cont:
        pic[inan] = off_scale
    else:
        pic[inan] = 0.

    #-- smoothing?

    if smooth:
        pic = gaussian_filter(pic, sigma=smooth)
        odmin = float(np.nanmin(pic))
        odmax = float(np.nanmax(pic))
    odrange = [odmin,odmax]

    if log:
        ok = np.where(pic > 0.)
        if len(ok[0]) == 0:
           print('PLOT_MAP: Warning, all data are negative. Cannot plot on a log scale. Using linear scale.')
           log = False
    if log:
        pmin = np.nanmin(pic[ok])
        pmax = np.nanmax(pic[ok])
        odrange = [np.float(pmin),np.float(pmax)]

    #-- establish plot labels

    units = get_map_prop(map,'units',default='arcsec')
    units = units
    xunits = units
    yunits = units
    xunits = get_map_prop(map,'xunits',xunits)
    yunits = get_map_prop(map,'yunits',yunits)
    xtitle='X ('+xunits+')'
    ytitle='Y ('+yunits+')'
    if nolabels:
        xtitle=''
        ytitle=''

    #-- if solar rotating, check that we are rotating relative to last time image
    #   rotated

    otime = get_map_prop(map, 'time')
    rtime = get_map_prop(map, 'rtime', default=otime)
    mtitle = get_map_prop(map, 'id', default='')
    if not over:
        mtime = rtime
        if original_time:
            mtime = otime
        if durtitle:
            mtitle += ' '+mtime
            s_time = Time(mtime).mjd
            e_time = Time(s_time + get_map_prop(map, 'dur', default=0.)/86400.,format='mjd').iso
            mtitle += '-'+e_time
        else:
            if date_only:
                date_obs = mtime[:10]
            else:
                date_obs = mtime+' UT'
            mtitle = mtitle+' '+date_obs
    mtitle = mtitle.strip()
    if title:
        mtitle = title
    if notitle:
        mtitle = ''

    #-- get some map properties

    oxrange = get_map_xrange(map, edge=True)
    oyrange = get_map_yrange(map, edge=True)
    dx = map['dx']
    dy = map['dy']
    dx2 = dx/2. 
    dy2 = dy/2.
    icenter = [map['xc'], map['yc']]
    curr_roll = get_map_prop(map, 'roll_angle', default=0.)
    curr_rcenter = get_map_prop(map, 'roll_center', default=icenter)

    #-- retrieve coordinate transformation angles for plotting limb and
    #   overlaying

    ang_error = ''
    #if map['id'].find('SOHO'):
    # Calculate the angles for this date in case they are not already in the map
    angles = pb0r(otime)
    # Override with the map values if they exist
    b0 = get_map_prop(map, 'b0', default=angles['b0'])
    l0 = get_map_prop(map, 'l0', default=angles['l0'])
    angles['rsun'] *= 60.
    rsun = get_map_prop(map, 'rsun', default=angles['rsun'])

    #-- establish plot ranges
    #   (start with image, then FOV, then XRANGE/YRANGE keywords)

    dcenter = icenter
    if valid_map(fov):
        dcenter = get_map_center(fov)
    if center:
        if valid_map(center):
            dcenter = get_map_center(center)
        else:
            dcenter = center

    dxrange = oxrange
    dyrange = oyrange
    if fov:
        if valid_map(fov):
            dfov = get_map_fov(fov, edge=True)
        else:
            dfov = 60.*np.array([fov[0],fov[-1]])
        half_fov = dfov/2.
        dxrange = [dcenter[0] - half_fov[0], dcenter[0] + half_fov[0]]
        dyrange = [dcenter[1] - half_fov[1], dcenter[1] + half_fov[1]]

    if center and not fov:
        dxrange[0] += dcenter[0] - icenter[0]
        dxrange[1] += dcenter[0] - icenter[0]
        dyrange[0] += dcenter[1] - icenter[1]
        dyrange[1] += dcenter[1] - icenter[1]

    #-- if overlaying, match with previous viewport

    if over:
        if not last['xrange'] is None: dxrange = last['xrange']
        if not last['yrange'] is None: dyrange = last['yrange']

    #-- overide with user input ranges

    if not xrange is None: dxrange = xrange
    if not yrange is None: dyrange = yrange
    dxrange.sort()
    if min(dxrange) == max(dxrange): dxrange = oxrange
    dyrange.sort()
    if min(dyrange) == max(dyrange): dyrange = oyrange

    #-- bail out if trying to display at the sub-pixel level

    diff_x = max(dxrange) - min(dxrange)
    diff_y = max(dyrange) - min(dyrange)
    if diff_x < dx2 or diff_y < dy2:
        print('PLOT_MAP: Error, Cannot display below half pixel resolution limit')
        return

    #-- define viewport

    xmin, xmax = dxrange
    ymin, ymax = dyrange

    if xmin == xmax or ymin == ymax:
        print('PLOT_MAP: Error, Plot scale MIN/MAX must differ')
        return

    #-- don't extract sub-region if contouring, since contour procedure
    #   takes care of it via drange

    if not cont:
        ranges = get_map_sub_ranges(map, xrange=dxrange, yrange=dyrange, xcor=None, ycor=None)
        irange = ranges['irange']
        if irange is None:
            return

        # NB: data array x and y dimensions of numpy arrays are swapped
        pic = pic[irange[0]:irange[1]+1, irange[2]:irange[3]+1]
        xmin, xmax, ymin, ymax = ranges['arange']

    #-- plot axes & viewport
    #-- try to preserve aspect ratio (won't work if multi is set)
    #-- if contouring and not overlaying, then check for roll

    no_drotate = not drotate
    no_project = no_drotate

    ilimb = False 
    olimb = False
    no_roll_correct = rolling or last['roll_correct']

    #-- get data plot limits

    if cont:
        if over:
            trans = [0,0]
            if about_center: 
                rcenter = last['center'] 
            else: 
                rcenter = last['rcenter']
            xp, yp = drot_map(map, time=last['time'], trans=trans, b0=last['b0'], l0=last['l0'],
                           rsun=last['rsun'], roll=last['roll'], rcenter=rcenter,
                           no_data = True, no_drotate=no_drotate, 
                           no_project=no_project, no_roll_correct=no_roll_correct, 
                           about_center=about_center)

#            ranges = get_map_sub_ranges(map, xrange=dxrange, yrange=dyrange, xcor=None, ycor=None)
#            xp = xp[irange[0]:irange[1]+1,irange[2]:irange[3]+1]
#            yp = yp[irange[0]:irange[1]+1,irange[2]:irange[3]+1]

            #-- send off-limb points to outside fov

            olimb = np.isnan(xp)
            if np.sum(olimb) != 0:
                pic[np.where(olimb)] = off_scale
            ilimb = (xp == -9999.)
            if np.sum(ilimb) != 0:
                pic[np.where(ilimb)] = off_scale
        else:
            xp = get_map_xp(map)
            yp = get_map_yp(map)
#            ranges = get_map_sub_ranges(map, xrange=dxrange, yrange=dyrange, xcor=None, ycor=None)
#            xp = xp[irange[0]:irange[1]+1,irange[2]:irange[3]+1]
#            yp = yp[irange[0]:irange[1]+1,irange[2]:irange[3]+1]

    prange = odrange

    #-- override with user keywords

    if dmin: prange[0] = dmin
    if dmax: prange[1] = dmax

    if drange:
        if valid_map(drange):
            prange = [np.nanmin(drange['data']),np.nanmax(drange['data'])]
        else:
            prange = np.float(drange)

    prange.sort()
    if min(prange) == max(prange): 
        prange = odrange

    if plus_only or minus_only:
        if plus_only:
            ok = np.logical_and(pic >= 0., pic != off_scale)
            nok = np.logical_not(ok)
        else:
            ok = np.logical_and(pic <= 0., pic != off_scale)
            nok = np.logical_not(ok)
        if np.sum(ok) == 0:
            if plus_only:
                print('PLOT_MAP: Error, All data are negative') 
            else:
                print('PLOT_MAP: Error, All data are positive') 
            return
        if np.sum(nok) > 0: 
            pic[np.where(nok)] = off_scale

    #-- log scale?

    if log:
        ok = np.logical_and(pic > 0., pic != off_scale)
        nok = np.logical_not(ok)
        if np.sum(ok) == 0:
            print('PLOT_MAP: Warning, all data are negative. Cannot plot on a log scale. Using linear scale.')
            log = False
        else:
            okidx = np.where(ok)
            pmin = np.nanmin(pic[okidx])
            pmax = np.nanmax(pic[okidx])
            if np.sum(nok) > 0:
                nokidx = np.where(nok)
                if cont:
                    pic[nokidx] = off_scale 
                else:
                    pic[nokidx] = pmin
            pic = np.log10(pic)
            if prange[0] <= 0: prange[0] = odrange[0]
            if prange[1] <= 0: prange[1] = odrange[1] 
            prange = np.log10(prange)

    #-- override with last scaling

    if last['scale']:
        if last['drange']:
            prange = last['drange']
        if log:
           prange = np.log10(last['drange'])

    #-- make an empty plot to establish scaling

    axis.set_xlim(xmin, xmax)
    axis.set_ylim(ymin, ymax)
    axis.set_aspect('equal','box')
    if no_data: return

    if cont:
        #-- plot contours
        nlevels = max(nlevels,2)
        dlevels = (prange[1] - prange[0])/(nlevels - 1.)
        def_levels = prange[0] + np.arange(nlevels)*dlevels
        plevels = def_levels
        if not levels is None:
            plevels = np.unique(levels)
            if percent:
                plevels = plevels*prange[1]/100. 
            elif log:
                ok, = np.where(plevels > 0)
                if len(ok) == 0:
                    print('PLOT_MAP: Warning, Contour levels must be greater than zero for log scale - using default set')
                    plevels = def_levels
                else: 
                    plevels = np.log10(plevels[ok])

        if not color is None:
            # For compatibility with IDL version, setting color keyword should be
            # interpreted as plotting all contours the same color
            kwargs.update({'colors':color})

        if 'fill' in kwargs.keys():
            kwargs.pop('fill')
            axis.contourf(xp, yp, pic, plevels, origin='lower', **kwargs)
        else:
            kwargs['linewidths'] = default(kwargs,'linewidths',0.5)
            kwargs['colors'] = default(kwargs,'colors','k')
            axis.contour(xp, yp, pic, plevels, origin='lower', **kwargs)
    else:
        #-- plot image
        if cbar:
            cbpos = default(kwargs,'cbpos','right')
            if cbpos == 'top' or cbpos == 'bottom':
                cbor = 'horizontal'
            else:
                cbor = 'vertical'

        im = axis.imshow(pic, origin='lower', extent=[xmin, xmax, ymin, ymax], vmin=prange[0], vmax=prange[1], **kwargs)

        # Add colorbar
        if cbar:
            from mpl_toolkits.axes_grid1 import make_axes_locatable
            divider = make_axes_locatable(axis)
            # if 'cbpos' in kwargs.keys():
                # cbpos = kwargs.pop('cbpos')
            # else:
                # cbpos = 'right'
            cax = divider.append_axes(cbpos, size=0.1, pad = 0.05)
            fig.colorbar(im, cax=cax, orientation=cbor, **kwargs)

    #-- plot axes and labels

    if not over:
        axis.set_xlabel(xtitle)
        axis.set_ylabel(ytitle)
        axis.set_title(mtitle)

    #-- overlay a solar latitude-longitude grid

    if (grid != 0. or limb):
        from mapping.plot_helio import plot_helio
        if shifting:
            rcenter = saved_rcenter 
        else: 
            rcenter = curr_rcenter
        plot_helio(otime, roll=curr_roll, axis=axis,
            over=True, rcenter=rcenter, no_roll_correct=no_roll_correct,
            grid=grid, glabel=glabel, gcolor=gcolor, gfont=gfont, gstyle=gstyle, gthick=gthick, 
            limb=limb, lcolor=lcolor, lthick=lthick, lstyle=lstyle,
            l0=l0)

    #-- mark point

#    mark_point,mark_point

    #-- plot border edges

#    if border then begin
#     !p.multi[0]=sp
#     oplot,xedge,yedge,thick=bthick,color=bcolor
#    endif

    #-- save last settings

    if not over:
        last.update({'xrange':dxrange})
        last.update({'yrange':dyrange})
        if not last['scale']:
            if log:
                prange = 10**prange
            last.update({'drange':prange})
        if not rtime is None: 
            last.update({'time':rtime}) 
        else: 
            last.update({'time':otime})
        if not curr_roll is None: last.update({'roll':curr_roll})
        if not curr_rcenter is None: last.update({'rcenter':curr_rcenter})
        if b0:        last.update({'b0':b0})
        if l0:        last.update({'l0':l0})
        if rsun:      last.update({'rsun':rsun})
        if not icenter is None:   last.update({'center':icenter})
        if rolling:   last.update({'roll_correct':True})

    if shifting or rolling:
        if valid_map(saved_map):
            map = saved_map

    # Write any changed information for the next call.
    fp = open(common_file,"wb")
    pickle.dump(last,fp)
    fp.close()
    sleep(0.1)

    return

def aia_lct(wave=None, load=False):
    ''' Functions like the parallel IDL version.  Returns the colormap for the given 
        wavelength snd addes the name 'aia'+str(wave) to the colortable list. 
        
        Inputs:
          wave    integer or string specifying the wavelength--must be one of 
                    'all', '1600', '1700', '4500', '94', '131', '171', '193', '211', '304', '335'
                    If 'all', then all of the aia colormaps are added to the colortable list.
          load    bool.  If True, it makes the colormap the default for subsequent plots. 
                    Note: If wave is 'all', load=True generates a warning and is ignored.
          
        Returns:
          cmap    A handle to the colormap, so that cmap=cmap in a plot command will use that
                    colormap.  But it also adds 'aia'+str(wave) to the colortable list, so
                    an alternative is to ignore the return value and just use, e.g. cmap='aia211'
                    Note: Returns None on error, or if wave is 'all'
    '''
    from matplotlib.pyplot import register_cmap, colormaps
    from matplotlib.colors import ListedColormap
    from matplotlib import rcParams
    # load in standard color table from which to start
    if wave is None:
        print('AIA_LCT: Syntax  cmap = aia_lct(wave[, load=True)')
        return None
        
    # allowed values of wave:
    wavelist = np.array(['1600', '1700', '4500', '94', '131', '171', '193', '211', '304', '335', 'all'])
    select, = np.where(str(wave) == wavelist)
    if len(select) == 0:
        print('AIA_LCT: Error, invalid wavelength/channel')
        return None

    # The below creates the standard IDL red colortable (3) 
    r0 = np.concatenate((np.linspace(0, 255, 256 - 80),np.ones(80)*255))/256.
    g0 = np.concatenate((np.zeros(120),np.linspace(0,255,256-120)))/256.
    b0 = np.concatenate((np.zeros(256-66),np.linspace(0,255,66)))/256.

    c0 = np.arange(256)/256.
    c1 = (np.sqrt(np.arange(256))*np.sqrt(255.))/256.
    c2 = (np.arange(256)**2/255.)/256.
    c3 = ((c1+c2/2.)*255./(np.max(c1)+np.max(c2)/2.))/256.
    a = np.ones(256)  # Alpha value is 1 for all

    rgb = np.ones((256, 4))
    if select == 10: select = np.arange(10)
    for sel in select:
        if sel == 0:   # 1600
            rgb = np.array([c3, c3,    c2, a]).T
        elif sel == 1: # 1700
            rgb = np.array([c1, c0,    c0, a]).T
        elif sel == 2: # 4500
            rgb = np.array([c0, c0, b0/2., a]).T
        elif sel == 3: # 94
            rgb = np.array([c2, c3,    c0, a]).T
        elif sel == 4: # 131
            rgb = np.array([g0, r0,    r0, a]).T
        elif sel == 5: # 171
            rgb = np.array([r0, c0,    b0, a]).T
        elif sel == 6: # 193
            rgb = np.array([c1, c0,    c2, a]).T
        elif sel == 7: # 211
            rgb = np.array([c1, c0,    c3, a]).T
        elif sel == 8: # 304
            rgb = np.array([r0, g0,    b0, a]).T
        elif sel == 9: # 355
            rgb = np.array([c2, c0,    c1, a]).T

        name = 'aia'+wavelist[sel]
        cmap = ListedColormap(rgb)
        if not name in colormaps():
            register_cmap(name, cmap)
            
    if wave == 'all':
        if load:
            print('AIA_LCT: Warning, no default color table loaded when wave="all"')
        return None
        
    if load:
        rcParams['image.cmap'] = name
        
    return cmap