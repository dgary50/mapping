from __future__ import print_function
from map_util import *
from transform_solar import arcmin2hel, hel2arcmin

from copy import deepcopy
import numpy as np
from astropy.time import Time
dtor = np.pi/180.

def rot_map(map=None, angle=None, rcenter=None, full_size=False,
                roll_angle=None, about_center=False, align_roll=False, no_copy=False):
    
    if not valid_map(map):
        print('ROT_MAP: Syntax rmap=rot_map(map,angle) OR rmap=rot_map(map,roll_angle=roll_angle)')
        print('% ANGLE = angle (deg clockwise to roll) or ROLL = new map ROLL angle')
        if not map is None:
            return map
        else:
            return False

    if angle is None and roll_angle is None:
        print('ROT_MAP: Error, Enter rotation angle in degrees clockwise from North')
        if not map is None:
            return map
        else:
            return False

    #-- don't rotate if multiple of 360.

    if angle:
        if (angle % 360.) == 0.: return map

    #-- read image and pixel arrays

    dx = map['dx']
    dy = map['dy']
    xc = map['xc']
    yc = map['yc']
    xp = get_map_xp(map)
    yp = get_map_yp(map)

    icenter = [xc, yc]
    roll_center = get_map_prop(map, 'roll_center', default = icenter)
    curr_roll = get_map_prop(map, 'roll_angle', default=0.)
    if valid_map(rcenter): 
        roll_center = get_map_center(rcenter)
    else:
        if rcenter: roll_center = rcenter
    if angle is None:
        angle = roll_angle - curr_roll

    if not no_copy:
        nmap = deepcopy(map) 
    else:
        nmap = map
        
    if about_center:
        roll_center = icenter
    new_roll = (angle + curr_roll) % 360

    nmap.update({'roll_angle':new_roll})
    nmap.update({'roll_center':roll_center})

    apply_roll = (angle % 360.) != 0.

    if apply_roll:
        rx, ry = roll_xy(xp,yp,angle,rcenter=roll_center)
        nmap = repack_map(nmap, rx, ry, no_copy=True)

        #-- rebin image
        #-- do this by regridding rotated coordinates and
        #   computing image data value in pre-rotated image by interpolation

        gxy = grid_xy(rx, ry, gspace=[dx,dy], adjust_resolution=full_size)
        gx = gxy['gx']
        gy = gxy['gy']
        rx, ry = roll_xy(gx, gy, -angle, rcenter=roll_center)
        xmin = np.nanmin(xp)
        xmax = np.nanmax(xp)
        ymin = np.nanmin(yp)
        ymax = np.nanmax(yp)
        out = np.logical_or(np.logical_or(np.logical_or(rx < xmin, rx > xmax), 
                                                        ry < ymin), ry > ymax)
        if np.sum(out) == np.prod(rx.shape):
            print('ROT_MAP: Error, No data in rotated image')
        else:
            find = np.isnan(nmap['data'])
            if np.sum(find) > 0:
                bad = np.where(find)
                nmap['data'][bad] = 0.
                rx[bad] = -999999.
                ry[bad] = -999999.
            from scipy.interpolate import RectBivariateSpline as interp2d
            sdata = deepcopy(nmap['data'])
            # zero-out a one-pixel boundary around the array to be rotated
            # so that extrapolations outside the boundary will be zero
            sdata[0] = 0.
            sdata[-1] = 0.
            sdata[:,0] = 0.
            sdata[:,-1] = 0.
            spl = interp2d(xp[0], yp[:,0], sdata.T)
            rdata = spl.ev(rx,ry)
            nmap['data'] = rdata
            nmap = repack_map(nmap, gx, gy, no_copy=True)

    if align_roll:
        nmap['roll_angle'] = 0.
        nmap['roll_center'] = map['roll_center']
        
    return nmap

def shift_map(map=None, xshift=0., yshift=0., xc=None, yc=None, no_copy=False, 
                        reverse=False, keep_roll_center=False):

    #-- check inputs 

    if not valid_map(map):
        print('SHIFT_MAP: Syntax smap = shift_map(map, xshift, yshift, [xc=xc, yc=yc])')
        if map is None:
            return False
        else:
            return map

    if no_copy:
        tmap = map
    else: 
        tmap = deepcopy(map)
        
    pxc = tmap['xc']
    pyc = tmap['yc']

    pxshift = 0.
    pyshift = 0.

    if xc is None:
        pxshift = xshift
    else:
        pxshift = xc - pxc

    if yc is None:
        pyshift = yshift
    else:
        pyshift = yc - pyc

    if reverse:
        pxshift = -pxshift
        pyshift = -pyshift

    if pxshift != 0. or pyshift != 0.:
        xrange = get_map_xrange(tmap)
        yrange = get_map_yrange(tmap)
        rcenter = tmap['roll_center']
        inside = np.logical_and(np.logical_and(np.logical_and(rcenter[0] <= max(xrange), 
                                                              rcenter[0] >= min(xrange)), 
                                                              rcenter[1] <= max(yrange)), 
                                                              rcenter[1] >= min(yrange))
        if np.sum(inside) > 0:
            tmap['roll_center'][0] += pxshift
            tmap['roll_center'][1] += pyshift
        tmap['xc'] += pxshift
        tmap['yc'] += pyshift

    return tmap

def sub_map(map=None, xrange=None, yrange=None, ref_map=None, preserve=False,
            noplot=False, irange=None, fov=None, pixel=None, 
            dimensions=None, project=False):

    if not valid_map(map):
        print('SUB_MAP: Syntax, smap = sub_map(map,xrange=xrange,yrange=yrange)')
        if map is None:
            return None
        else:
            return map

    if not valid_map(ref_map) and valid_map(fov):
        ref_map = fov

    do_plot = xrange is None and yrange is None and pixel is None

    #-- get data coordinates to extract
    #-- get arcsec ranges of sub map

    if do_plot:
        print('SUB_MAP: Error, Manual selection of region is not yet available')
        return map
        #plot_map(map)
        #region = get_sub_region(pixel=pixel)
        # xrange=[region[0],region[1]]
        # yrange=[region[2],region[3]]

    #-- get pixel indicies (REF_MAP overrides input or plot-derived x/y ranges)

    ny, nx = map['data'].shape
     
    if pixel is None:
        if valid_map(ref_map):
            xrange = get_map_prop(ref_map, 'xrange')
            yrange = get_map_prop(ref_map, 'yrange')

            if project:

                #-- compute heliographic center of map in coordinate frame of
                #   reference map, allowing for differential rotation.
       
                h2 = proj_xy(ref_map, map)

                #-- compute sub FOV around this center

                dxs = (max(xrange) - min(xrange))/2.
                dys = (max(yrange) - min(yrange))/2.
                xrange = [h2[0] - dxs, h2[0] + dxs]
                yrange = [h2[1] - dys, h2[1] + dys]

        ranges = get_map_sub_ranges(map, xrange=xrange, yrange=yrange)
        
        irange = ranges['irange']
        if irange[1]-irange[0] <= 2 and irange[3]-irange[2] <=2:
            print('SUB_MAP: Error, Insufficient points in sub-region. Skipping map.')
            return map

        x1, x2, y1, y2 = irange

    else:
        # Ranges are in pixel units
        if xrange is None and yrange is None:
            print('SUB_MAP: Error, Non-zero XRANGE or YRANGE must be entered for pixel selection.')
            return map
        # Swap the order of the pixel ranges (x becomes y and vice versa), 
        # due to Python array order
        y1, y2 = np.clip(xrange, 0, nx-1)
        x1, x2 = np.clip(yrange, 0, ny-1)

    #-- keep same pixel dimensions as REF_MAP map

    if valid_map(ref_map) and not preserve is None:
        nrx, nry = get_map_size(ref_map)
        fspace = [get_map_prop(ref_map,'dx'),get_map_prop(ref_map,'dy')]
        dspace = [get_map_prop(map,'dx'),get_map_prop(map,'dy')]
        d1 = dspace[0]-fspace[0]
        d2 = dspace[1]-fspace[1]
        if abs(d1) >= dspace[0] or abs(d2) >= dspace[1]:
            print('SUB_MAP: Warning, input image and reference image have different pixel spacings')
        x2 = x1 + nrx - 1
        y2 = y1 + nry - 1

    #-- if user entered DIMENSIONS keyword, then they must really want to control
    #   the output dimensions of the sub map.
     
    if dimensions:
        x2 = x1 + dimensions[0] - 1
        y2 = y1 + dimensions[1] - 1
        if x2 >= nx or y2 >= ny:
            print('SUB_MAP: Error, Specified dimensions of sub map exceed input map size. Skipping map')
            return map

    smap = deepcopy(map)
    # Definitions of x, y are swapped for Python
    smap.update({'data':map['data'][x1:x2+1,y1:y2+1]})
    xp = get_map_xp(map)
    yp = get_map_yp(map)
    cenx = get_arr_center(xp[x1:x2+1,y1:y2+1])
    ceny = get_arr_center(yp[x1:x2+1,y1:y2+1])
    smap.update({'xc':cenx['center'],'yc':ceny['center'],'dx':cenx['dx'],'dy':ceny['dy']})

    return smap
    
def diff_rot(ddays, latitude, howard=True, allen=False, synodic=False, 
             sidereal=True, rigid=False, rate=None, carrington=False, snodgrass=False):

    try:
        # Convert latitude list to an array if possible (this does nothing if already an array)
        npts = len(latitude)
        latitude = np.array(latitude)
    except:
        # Looks like inputs are not arrays
        pass

#-- check if rotating as rigid body

    if rigid:
        sin2l = np.zeros_like(latitude)
        if rate and rate > 0: return ddays*rate + sin2l
        else:
            print('DIFF_ROT: If rigid is True, must supply a rotation rate [deg/day]')
            return None
    else:
        sin2l = np.sin(latitude*dtor)**2
        sin4l = sin2l*sin2l

    if allen:
        #  Allen, Astrophysical Quantities
        rotation = ddays*(14.44 - 3.*sin2l)
    elif snodgrass:
        #  Magnetic features as used by the Solar Orbiter project for planning
        #  (Snodgrass and Ulrich, 1990, Ap. J., 351, 309-316)
        rotation = ddays*(14.252 - 1.678*sin2l - 2.401*sin4l)
    else:
        #  Small magnetic features 
        #  (Howard, Harvey, and Forgach, Solar Physics, 130, 295, 1990)
        rotation = (1.e-6)*ddays*(2.894 - 0.428*sin2l - 0.37*sin4l)*86400./dtor

    if synodic:
        rotation -= 0.9856*ddays
    elif carrington:
        rotation -= (360./25.38)*ddays

    return rotation


def rot_xy(xx=None, yy=None, interval=None, tstart=None, tend=None,
            keep=False, offlimb=False, index=None,
            back_index=None, vstart=None, vend=None,
            return_backside=False, **kwarg):
    ''' Perform differential solar rotation of point(s) (xx,yy) and return the rotated
        coordinates plus information about visibility of the points as a result of the
        rotation.
        
        Inputs:
            xx       initial heliocentric X position (arcsec) or vector of such positions
            yy       initial heliocentric Y position (arcsec) or vector of such positions
            tstart   time (as Time() object or string in ISO format) for which xx and yy
                        apply
            tend     time (as Time() object or string in ISO format) for which xx and yy
                        are differentially rotated (tend may be omitted if "interval" 
                        is supplied)
            interval optional duration of rotation, in s, used if tend is not supplied
        
        Keyword Parameters:
            keep     if True, keep the same epoch date for P, B0, and RSUN for both
                        start and end of rotation
        Returns:
            xydict   dictionary of outputs with keys:
                        'rx'    = rotated positions of xx
                        'ry'    = rotated positions of yy
                        'vis'   = indicates whether points are still visible after
                                    rotation (True) or are on the backside (False)
    '''
    from astropy.time import Time
    
    if xx is None or yy is None:
        print('ROT_XY: Error, Syntax: a = rot_xy(x, y, interval)')
        return None

    badinputs = type(xx) != type(yy)
    if not badinputs:
        if type(xx) is list:
            badinputs = len(xx) != len(yy)
        else:
            try:
                ans = xx*yy
            except:
                bandinputs = True
    if badinputs:
        print('ROT_XY: Error, XX and YY must be scalar or vector with the same elements')
        return None

    #-- validate time inputs

    if tstart:
        try:
            if type(tstart) is str: tstart.replace('/','-')
            dstart = Time(tstart).tai
        except:
            print('ROT_XY: Error, could not interpret tstart time as ISO time.')
            return None
    if tend:
        try:
            if type(tend) is str: tend.replace('/','-')
            dend = Time(tend).tai
        except:
            print('ROT_XY: Error, could not interpret tend time as ISO time.')
            return None
    else:
        if interval is None:
            print('ROT_XY: Error, either interval or tend must be supplied.')
            return None
        else:
            dend = dstart + interval/86400.

    # dstart and dend are now valid TAI Time() objects

    # Declare output variables
    rx = np.zeros_like(xx) + np.nan
    ry = np.zeros_like(yy) + np.nan
    
    #---------------------------------------------------------------------------
    #  Compute heliographic coordinates of initial time
    #---------------------------------------------------------------------------
    toarcmin = False
    if vstart['rsun'] > 500:
        toarcmin = True
        vstart['rsun'] /= 60.  # Convert rsun to arcmin for this call
    llout = arcmin2hel(xx/60., yy/60., date=dstart.iso, vdict=vstart, **kwarg)
    if toarcmin: vstart['rsun'] *= 60.  # Convert back to arcsec
    #---------------------------------------------------------------------------
    #  If time interval is zero, bail out unless Solar angles passed via VSTART/VEND
    #---------------------------------------------------------------------------

    ang_in = type(vstart) is dict or type(vend) is dict
    if (dend - dstart) == 0. and not ang_in:
        print('ROT_XY: Zero time interval. No need to rotate.')
        return {'rx':xx, 'ry':yy, 'vis':llout['offlimb']}

    if sum(llout['offlimb']) == len(llout['offlimb']):
        print('ROT_XY: Error, All initial points are off the limb; cannot rotate!')
        return {'rx':xx, 'ry':yy, 'vis':llout['offlimb']}

    #-- differentially rotate if interval is non-zero

    if (dend - dstart) != 0.:
        ddays = (dend - dstart).value
        llout['lon'] += diff_rot(ddays, llout['lat'], synodic=True)

    if keep:
        datstr = dstart.iso 
    else:
        datstr = dend.iso
    back_index = True
    off = llout['offlimb']
    if np.sum(off) == 0:
        # All points are on the solar disk
        toarcmin = False
        if vend['rsun'] > 500:
            toarcmin = True
            vend['rsun'] /= 60.  # Convert rsun to arcmin for this call
        xyout = hel2arcmin(llout['lat'], llout['lon'], date=datstr, vdict=vend)
        if toarcmin: vend['rsun'] *= 60.  # Convert back to arcsec
        rx = xyout['xx'] * 60.
        ry = xyout['yy'] * 60.
        visible = xyout['visible']
        back = np.logical_not(xyout['visible'])
        if np.sum(back) > 0:
            if return_backside:
                # Flag frontside points
                idx, = np.where(not back)
            else:
                # Flag backside points
                idx, = np.where(back)
            rx[idx] = -9999.0
            ry[idx] = -9999.0
    else:
        on = np.logical_not(off)   # True for points on visible disk of Sun
        if np.sum(on) > 0:
            # Some points are on the solar disk
            # Only call this for points on the visible disk
            lat = llout['lat'][on]
            toarcmin = False
            if vend['rsun'] > 500:
                toarcmin = True
                vend['rsun'] /= 60.  # Convert rsun to arcmin for this call
            xyout = hel2arcmin(llout['lat'][on], llout['lon'][on], date=datstr, vdict=vend)
            if toarcmin: vend['rsun'] *= 60.  # Convert back to arcsec
            rx[on] = xyout['xx'] * 60.
            ry[on] = xyout['yy'] * 60.
            back = np.logical_not(xyout['visible'])
            if np.sum(back) > 0:
                bad, = np.where(on)
                if return_backside:
                    # Flag frontside points
                    really_bad, = np.where(not back)
                    idx = bad[really_bad]
                else:
                    # Flag backside points
                    really_bad, = np.where(back)
                    idx = bad[really_bad]
                rx[idx] = -9999.0
                ry[idx] = -9999.0
            visible = np.logical_not(llout['offlimb'])
            visible[on] = np.logical_and(visible[on],xyout['visible'])
                        
    return {'rx':rx, 'ry':ry, 'vis':visible}


def proj_xy(map=None, rmap=None, drotate=True):

    if not valid_map(map) and not valid_map(rmap):
        print('PROJ_XY: Syntax pcor = proj_xy(map, rmap, drotate=True)')
        return None, None

    xc = map['xc']
    yc = map['yc']

    vstart = get_map_angles(map)
    vend = get_map_angles(rmap)
    tstart = map['time']
    tend = rmap['time']

    #-- if not on disk then bail

    on_disk = np.sqrt(xc**2 + yc**2) < vstart['rsun']
    if not on_disk:
        print('PROJ_XY: Error, Center of input map is off limb. Cannot project.')
        return xc, yc

    # #-- project

    if not drotate:
        tend = tstart

    rx, ry = roll_xy(xc, yc, -map['roll_angle'], rcenter=map['roll_center'])
    rcor = rot_xy(rx, ry, tstart=tstart, tend=tend, vstart=vstart, vend=vend, sphere=True)
    xc = rcor['rx']
    yc = rcor['ry']

    rx, ry = roll_xy(xc, yc, rmap['roll_angle'], rcenter=rmap['roll_center'])
    on_disk = np.sqrt(xc**2 + yc**2) < vend['rsun']

    if not on_disk:
        print('PROJ_XY: Error, Center of projected map is outside field-of-view.')

    return xc, yc

def drot_map(inmap=None, duration=None, proj_rcenter=False, time=None, 
                  trans=None, resolution=None, roll=None, no_rtime=False,
                  rcenter=None, center=None, keep_limb=False, unrotated=False, 
                  no_drotate=False, ref_map=None, dimensions=None, outsize=None,
                  same_center=False, degrees=False, track_center=False,
                  preserve_area=False, adjust_resolution=False, xrange=None,
                  yrange=None, b0=None, l0=None, rsun=None, xp=None, yp=None, 
                  no_project=False, no_data=False, no_roll_correct=False,
                  about_center=False, align_center=False):

    if not valid_map(inmap):
        print('DROT_MAP: Syntax rmap = drot_map(map, duration, [time=time])')
        return None

    deg_per_day = diff_rot(1, 0, synodic=True)
    sec_per_day = 24.*3600.
    sec_per_deg = sec_per_day/deg_per_day

    # Make a copy to modify for output
    map = deepcopy(inmap)
    
    tdur = None
    if degrees:
        if duration is None:
            print('DROT_MAP: Error, Must provide an amount of rotation in hours(or degrees if degrees=True)')
            return map
        if degrees:
            tdur = 24.*duration/deg_per_day
    elif duration:
        tdur = duration

    #-- check keywords

    track_center = track_center or align_center
    if preserve_area or adjust_resolution: track_center = False

    #-- if REF_MAP entered then use it's TIME, SPACING, CENTER, ROLL, and
    #    DIMENSIONS

    etime = None
    if valid_map(ref_map): 
        etime = Time(ref_map['time'].replace('/','-')).tai
        xc = ref_map['xc']
        yc = ref_map['yc']
        dx = ref_map['dx']
        dy = ref_map['dy']
        ny, nx = ref_map['data'].shape
        droll_center = ref_map['roll_center']
        droll = ref_map['roll_angle']
        dspace = [dx,dy]
        dcenter = [xc,yc]
        dsize = [nx,ny]
    if valid_map(time): 
        etime = Time(time['time'].replace('/','-')).tai
    elif time: 
        etime = Time(time.replace('/','-')).tai
            
    #-- translate after rotation?

    xs = 0.
    ys = 0.
    do_trans = False
    if trans:
        xs, ys = trans
        do_trans = xs != 0. or ys != 0.

    sub_range = False
    if not xrange is None and not yrange is None: 
        sub_range = True

    #-- input data type is less than float, then make it float for better
    #   precision

    #-- get solar rotation duration

    cur_time = Time(map['time'].replace('/','-')).tai
    if etime is None:
        if tdur is None:
            print('DROT_MAP: Error, Must supply one of duration, reference map, or time.')
            return map
        etime = cur_time + tdur/24.
    cdur = etime - cur_time

    # #-- check if differentially rotating

    if no_drotate: cdur *= 0.
    new_time = cur_time + cdur
    do_drot = new_time != cur_time

    if do_drot and (cdur.value*sec_per_day > 180*sec_per_deg):
        print('DROT_MAP: Warning, most of Sun will rotate over limb')

    #-- get start and end projection angles
    #-- override with b0, l0, or rsun entered as keywords

    start_angles = get_map_angles(map)
    end_angles = start_angles
    if do_drot:
        end_angles = get_map_angles(map, new_time.iso)
    if b0: end_angles['b0'] = b0
    if l0: end_angles['l0'] = l0
    if rsun: end_angles['rsun'] = rsun
    if no_project: 
        do_proj = False
    else:
        do_proj = False
        for key in ['b0','l0','rsun']:
            if start_angles[key] != end_angles[key]:
                do_proj = True

    # #--- if not projecting, ensure that different RSUN perspectives match

    do_rad = False
    if not do_proj and 'rsun' in map.keys():
        ratio = end_angles['rsun']/start_angles['rsun']
        if ratio != 1.:
            print('DROT_MAP: Info: Adjusting map to match solar distance in reference map.')
        map['dx'] *= ratio
        map['dy'] *= ratio
        map['xc'] *= ratio
        map['yc'] *= ratio
        map['roll_center'] = [i*ratio for i in map['roll_center']]
        map['rsun'] = end_angles['rsun']
        do_rad = True

    #-- extract the map data

    if sub_range:
        map = sub_map(map, xrange=xrange, yrange=yrange)

    if valid_map(map): 
        xc = map['xc']
        yc = map['yc']
        dx = map['dx']
        dy = map['dy']
        ny, nx = map['data'].shape
        curr_rcenter = map['roll_center']
        curr_roll = map['roll_angle']
        xp = get_map_xp(map)
        yp = get_map_yp(map)
        pxrange = deepcopy(xrange)
        pyrange = deepcopy(yrange)

    if about_center: curr_rcenter = [xc, yc]
     
    #-- check if rolling

    have_roll = (curr_roll % 360.0) != 0. 
    new_roll = curr_roll
    try:
        new_roll = droll
    except:
        pass
    if roll: new_roll = roll
    roll_diff = new_roll - curr_roll

    do_roll = not no_roll_correct and (roll_diff % 360.) != 0.

    #-- check if new roll center

    new_rcenter = curr_rcenter
    try:
        new_rcenter = droll_center
    except:
        pass
    if not rcenter is None: new_rcenter = rcenter

    do_rcenter = False
    if do_roll and have_roll:
        do_rcenter = (new_rcenter[0] != curr_rcenter[0] or
                      new_rcenter[1] != curr_rcenter[1])

    #-- check if recentering 

    do_center = False
    new_center = None
    try:
        if dcenter: new_center = dcenter
    except:
        pass
    if center: new_center = center
    if new_center: do_center = (new_center[0] != xc) or (new_center[1] != yc)

    #-- check if resizing

    new_size = [nx,ny]
    try:
        if dsize: new_size = dsize
    except:
        pass
    if outsize: new_size = outsize
    if dimensions: new_size = dimensions
    do_resize = (new_size[0] != nx) or (new_size[1] != ny) or preserve_area
      
    #-- check if rebinning

    new_space = [dx,dy]
    if do_resize:
        new_space[0] = dx*nx/new_size[0]
        new_space[1] = dy*ny/new_size[1]
    try:
        if dspace: new_space = dspace
    except:
        pass
    if resolution:
        try:
            new_space = [resolution[0], resolution[1]]
        except:
            new_space = [resolution, resolution]
    do_rebin = (new_space[0] != dx) or (new_space[1] != dy) or adjust_resolution

    onx, ony = new_size
    map.update({'roll_angle':0})
    map.update({'roll_center': [xc,yc]})

    if not do_proj and not do_drot and not do_rcenter and not do_trans and \
       not do_rebin and not do_roll and not do_center and not do_resize and not do_rad: 
        print('DROT_MAP: Nothing to do!')
        if no_data:  
            return xp, yp
        else:
            return map

    #-- get the before and after solar radii since we will need these to
    #   flag offlimb points when rotating or projecting.

    if do_proj or do_drot:
        sol_rad1 = start_angles['rsun']
        sol_rad2 = end_angles['rsun']
        sol_rat = sol_rad1 / sol_rad2

    if do_roll:
        map['roll_angle'] = new_roll
        map['roll_center'] = new_rcenter
    else:
        if have_roll: 
            map['roll_angle'] = curr_roll

    #-- if image is rolled and roll-center is within image, then
    #   project roll-center if /proj_rcenter

            roll_in_image = ((curr_rcenter[0] <= max(pxrange) and
                             (curr_rcenter[0] >= min(pxrange))) or
                            ((curr_rcenter[1] <= max(pyrange)) and
                             (curr_rcenter[1] >= min(pyrange))))
            if roll_in_image and (do_proj or do_drot) and proj_rcenter:
                out = rot_xy(curr_rcenter[0], curr_rcenter[1], tstart=cur_time,
                             tend=new_time, sphere=True, vstart=start_angles, 
                             vend=end_angles)
                drot_rcenter = [out['rx'],out['ry']]
            else: drot_rcenter = curr_rcenter
            map['roll_center'] = drot_rcenter

    #-- correct current roll before projecting or drotating

    xr = deepcopy(xp)
    yr = deepcopy(yp)

    if have_roll: xr, yr = roll_xy(xp, yp, -curr_roll, rcenter=curr_rcenter)
     
    #-- flag offlimb pixels

     # icount=0 & ocount=0 & olimb=-1 & 
    fsize = nx*ny
    fov = None
    if do_proj or do_drot:
        rad1 = np.sqrt(xr**2 + yr**2)
        outside = rad1 > sol_rad1
        if np.sum(outside) == fsize: 
            print('DROT_MAP: Error, All points off limb, cannot project')
            return map
        #-- apply solar rotation/projection

        xyshape = xr.shape
        xr = xr.flatten()
        yr = yr.flatten()
        out = rot_xy(xr, yr, tstart=cur_time, tend=new_time, sphere=True,
                      vstart = start_angles, vend=end_angles)
        if out is None:
            print('DROT_MAP: Error in ROT_XY.  Cannot continue.')
            return map
        xr = out['rx'].reshape(xyshape)
        yr = out['ry'].reshape(xyshape)
        
        #-- flag pixels that projected over limb

        rad2 = np.sqrt(xr**2 + yr**2)
        rad2[np.isnan(rad2)] = -9999.**2
        outside = rad2 > sol_rad2
        if np.sum(outside) == fsize: 
            print('DROT_MAP: Error, All points off limb at end time.')
            return map

    #-- determine valid pixels still on disk

        fov = np.logical_and(rad1 <= sol_rad1, rad2 <= sol_rad2)
        if np.sum(fov) == 0:
            print('DROT_MAP: Error, All points projected outside original FOV')
            return map

    #-- apply translation

    xr += xs
    yr += ys

    #-- apply roll
      
    if do_roll: 
        xr, yr = roll_xy(xr, yr, new_roll, rcenter=new_rcenter)
    elif have_roll:  
        xr, yr = roll_xy(xr, yr, curr_roll, rcenter=drot_rcenter)

    #-- return if just need coordinates
     
    if no_data:
        return xr, yr

    #-- update map properties 

    bad = np.where(xr == -9999.)
    xr[bad] = np.nan
    bad = np.where(yr == -9999.)
    yr[bad] = np.nan    
    map = repack_map(map, xr, yr, no_copy=True)
    bad = np.where(np.isnan(xr))
    xr[bad] = -9999.
    bad = np.where(np.isnan(yr))
    yr[bad] = -9999
      
    #-- remap image

    #-- first make a regularized grid using only pixels that are still in fov
    #   (i.e. limb pixels and disk pixels that haven't projected over limb)

    if same_center: new_center = [xc, yc]

    #-- track FOV center 

    if track_center and (do_proj or do_drot):
        xcen = xc
        ycen = yc
        if have_roll: xcen, ycen = roll_xy(xcen, ycen, -curr_roll, rcenter=curr_rcenter)
        out = rot_xy(xcen, ycen, tstart=cur_time, tend=new_time, sphere=True,
                         vstart=start_angles, vend=end_angles)
        xcen, ycen = [out['rx'], out['ry']]
        if do_roll: 
            xcen, ycen = roll_xy(xcen, ycen, new_roll, rcenter=new_rcenter)
        elif have_roll:
            xcen, ycen = roll_xy(xcen, ycen, curr_roll, rcenter=drot_rcenter)
        new_center = [xcen,ycen]

    if not fov is None:
        xr = xr[fov]
        yr = yr[fov]

    out = grid_xy(xr, yr, gspace=new_space, center=new_center, gsize=new_size, 
                  preserve_area=preserve_area, adjust_resolution=adjust_resolution)

    ony, onx = out['gsize']
    do_resize = onx != nx or ony != ny
    gx = out['gx']
    gy = out['gy']
    
    #-- project grid points back to find where each point came from
    map = repack_map(map, gx, gy, no_copy=True)
    xr = gx
    yr = gy
    
    #-- roll back 

    if do_roll: 
        xr, yr = roll_xy(xr, yr, -new_roll, rcenter=new_rcenter)
    elif have_roll:
        xr, yr = roll_xy(xr, yr, -curr_roll, rcenter=drot_rcenter)

    #-- shift back

    xr -= xs
    yr -= ys

    #-- project backwards

    if do_proj or do_drot:

    #-- flag projected limb pixels 

        rad2 = np.sqrt(xr**2 + yr**2)
        outside = rad2 > sol_rad2
        if keep_limb and np.sum(outside > 0):
            xlimb = deepcopy(xr[outside])
            ylimb = deepcopy(yr[outside])

        xr = xr.flatten()
        yr = yr.flatten()
        out = rot_xy(xr, yr, tstart=new_time, tend=cur_time, sphere=True,
                     vstart=end_angles, vend=start_angles)
        if out is None:
            print('DROT_MAP: Error in ROT_XY.  Cannot continue.')
            return map
        xr = out['rx'].reshape(ony, onx)
        yr = out['ry'].reshape(ony, onx)

        if keep_limb and np.sum(outside > 0):
            xr[outside] = xlimb*sol_rat
            yr[outside] = ylimb*sol_rat

    #-- roll back to initial roll

    if have_roll: xr, yr = roll_xy(xr, yr, curr_roll, rcenter=curr_rcenter)

    nans = np.isnan(map['data'])
    if np.sum(nans) > 0:
        map['data'][nans] = 0.

    from scipy.interpolate import  RectBivariateSpline as interp2d
    sdata = deepcopy(map['data'])
    # zero-out a one-pixel boundary around the array to be transformed
    # so that extrapolations outside the boundary will be zero
    sdata[0] = 0.
    sdata[-1] = 0.
    sdata[:,0] = 0.
    sdata[:,-1] = 0.
    spl = interp2d(xp[0], yp[:,0], sdata.T)
    bad = np.where(np.isnan(xr))
    xr[bad] = -9999.
    bad = np.where(np.isnan(yr))
    yr[bad] = -9999.
    rdata = spl.ev(xr,yr)

    map['data'] = rdata
    bad = np.where(xr == -9999.)
    xr[bad] = np.nan
    bad = np.where(yr == -9999.)
    yr[bad] = np.nan

    if do_proj or do_drot:
        dmin = np.nanmin(map['data'])
        dmax = np.nanmax(map['data'])
        if dmin == dmax  and dmin == 0.:
            print('DROT_MAP: Warning, Image projected out of field of view')

    rtime = map['time']
    if do_drot: rtime = new_time.iso
    if do_proj:
        map.update(end_angles)
    if not no_rtime and (do_drot or do_proj): map.update({'rtime':rtime})    

    return map