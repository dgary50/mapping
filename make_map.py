#+
# Project     : SOHO-CDS
#
# Name        : MAKE_MAP
#
# Purpose     : Make an image map 
#
# Category    : imaging
#
# Syntax      : map=make_map(data)
#
# Inputs      : DATA = 2d data array
#
# Outputs     : 
#               MAP ={'data':data,'xc':xc,'yc':yc,'dx':dx,'dy':dy,'time':time}
#               where,
#               DATA  = 2d image array
#
# Opt. Outputs: None
#
# Keywords    : XC,YC = center of image [arsces]
#               DX,DY = pixel spacing in X and Y directions [arcsecs]
#               TIME = image time [UT format]
#               ID = unique string identifier
#               SUB   = [x1,x2,y1,y2] = indicies of sub-array to extract
#               DATA_UNITS = units of subarray in data units [def=pixel number]
#               ROLL_ANGLE = image roll (deg clockwise from N) [optional]
#               ROLL_CENTER = roll center [optional]
#               FOV = same as SUB, but with /DATA_UNITS
#               NO_COPY= set to not make new copy of data
#               BYTE_SCALE = bytescale data
#
# History     : 2020-May-21  DG
#                 Modified from Dominic Zarro's IDL code of the same name
#
# Contact     : dgary@njit.edu
#-

from __future__ import print_function
from astropy.time import Time
import numpy as np
from map_util import *

def filt_data(data=None, min=None, max=None, inverse=False, nan=False, 
              missing=None, positive=False):
    ''' Returns array of boolean flags (True indicates bad data if inverse=True or good data 
        if inverse=False, defined by min, max, nan, missing, or positive arguments.
        
    '''
    above = False
    if data is None:
        print('FILT_DATA: Syntax: fdata = filt_data(data,[/pos,missing=missing,dmin=dmin,dmax=dmax])')
        return False

    if inverse:
        # Start with all zeros
        flag = np.zeros_like(data, dtype='bool')
        if positive:
            flag = np.logical_or(flag, data <= 0)
        if missing:
            flag = np.logical_or(flag, data == missing)
        if min:
            flag = np.logical_or(flag, data < min)
        if max:
            flag = np.logical_or(flag, data > max)
        if nan:
            flag = np.logical_or(flag, np.isfinite(data))
    else:
        # Start with all ones
        flag = np.ones_like(data, dtype='bool')
        if positive:
            flag = np.logical_and(flag, data > 0)
        if missing:
            flag = np.logical_and(flag, data != missing)
        if min:
            flag = np.logical_and(flag, data >= min)
        if max:
            flag = np.logical_and(flag, data <= max)
        if nan:
            flag = np.logical_and(flag, np.isfinite(data))

    return flag

def dscale(array=None, min=None, max=None, missing=None, log=False, err='', 
           no_copy=False, nan=False, drange=None):
    ''' Routine to scale a 2d array
    '''
    if array is None:
        print('data[,nok] = dscale(data, [/log, missing=missing][, nok=True])')
        err = 'invalid input'
        return [False,False,False,None]

    #-- set data limits
    amax = np.nanmax(array)
    amin = np.nanmin(array)
    if min: 
        cmin = min 
    else: 
        cmin = amin
    if max:
        cmax = max 
    else: 
        cmax = amax
        
    if drange:
        try:
            cmin, cmax = drange
            temp = np.array([cmin,cmax]).astype(np.float)
            temp.sort()
            cmin, cmax = temp
        except:
            print('DSCALE: Data range (drange) invalid.  Ignored')
            pass

    nx, ny = array.shape

    if (amax == amin) and (amax == 0):
        print('DSCALE: Error -- No data in specified range.')
        return [None,None,None,None]
    
    if log:
        log = True

    #-- flag missing or negative/zero data values for log case

    above = False
    if cmax != amax:
        above = array > cmax
        acount = np.sum(above)
        
    do_filt = missing or log or (cmin != amin) or (cmax != amax) or nan

    nok = False
    if do_filt:
        nok = filt_data(array, miss=missing, positive=log, min=cmin, max=cmax,
                               nan=nan, inverse=True)
        count = np.sum(nok)
        npts = np.prod(array.shape)
        if count == npts:
            print('DSCALE: Error -- No data in specified range.')
            return [None,None,None,None]

    #-- create output array
                  
    if no_copy:
        darray = array
    else:
        from copy import copy
        darray = copy(array)

    if nok:
        darray[nok] = abs(amax) + 1

    if log:
        darray = np.log10(darray)
        if cmin > 0: 
            cmin = np.log10(cmin)
        else:
            cmin = np.nanmin(darray)
        if cmax > 0:
            cmax = np.log10(cmax)
        else: 
            cmax = np.nanmax(darray)

    #dprint,'% DSCALE: cmin, cmax: ',cmin,cmax

    crange = np.array([cmin, cmax])

    if cmax == cmin:
        print('DSCALE: Error -- No data in specified range.')
        return [None,None,None,crange]

    return [darray, nok, above, crange]

def cscale(array, top=False, bottom=False, reverse=False, no_copy=False):

    carray, nok, above, crange = dscale(array, no_copy=no_copy)

    if carray is None: 
        return [carray, None]

    #-- set color limits

    cmax = np.max(crange)
    cmin = np.min(crange)
    lmax = 255
    if top:
        ctop = np.min(top,lmax)
    else:
        ctop = lmax
    if bottom:
        cbot = np.max(bottom,0.)
    else:
        cbot = 0. 
    ctop = float(ctop)
    cbot = float(cbot)

    if reverse:
        ctop, cbot = [cbot, ctop]

    slope = (ctop - cbot)/(float(cmax) - float(cmin))
    off = ctop - slope*float(cmax)
    carray = slope*carray.astype(float) + off
    carray = np.rint(carray)
    if nok:
        carray[nok] = cbot
    chk = np.logical_or(carray < cbot, carray > ctop)
    if np.sum(chk) > 0:
        carray[chk] = cbot
    if above:
        carray[above] = ctop
    carray = carray.astype('uint8')

    obottom = np.uint8(cbot)
    otop = np.uint8(ctop)
    return [carray, [obottom, otop]]

def make_map(data, xcen=0.0, ycen=0.0, dx=1.0, dy=1.0, time=None,
             sub=None, roll_angle=0.0, id='Generic Map', data_units=False,
             roll_center=None, fov=None, xunits='arcsec', extra=None,
             yunits='arcsec', dur=0.0, no_copy=False, byte_scale=False,
             xrange=None, yrange=None):

    #-- check input image
    map=0
    sz=data.shape
    if len(sz) != 2:
        print('Input image must be 2-d')
        return {}

    #-- check time
    try:
        # Try to interpret time as either an iso string or a Time() object
        t = Time(time)
        stime = t.iso
    except ValueError:
        # Error interpreting time, so set to current time.
        stime = Time.now().iso

    nx = sz[0]; ny = sz[1]
    dx = np.float64(dx); dy = np.float64(dy)
    dx = np.abs(dx); dy = np.abs(dy)

    #-- extract subarray?
    do_sub = False
    if not fov is None:
        if len(fov) == 4:
            subreg = fov
            data_units = True
            do_sub = True
    elif not sub is None:
        if len(sub) == 4:
            subreg = sub
            do_sub = True

    xc = np.float64(xcen); yc = np.float64(ycen)
    if do_sub:
        xp = mk_map_xp(xc,dx,nx,ny)
        yp = mk_map_yp(yc,dy,nx,ny)

    #-- start building map structure

    #bscale=keyword_set(byte_scale)
    if byte_scale:
        bdata, brange = cscale(data, no_copy=no_copy)
    else:
        bdata = None

    if do_sub:
        if data_units:
            sub = get_map_region(xp, yp, subreg)
            print(sub)
            if sub[0]:
                if (sub[1]-sub[0]) < 2 or (sub[3]-sub[2]) < 2:
                    print('MAKE_MAP: Insufficient data points to produce map')
                    return {}

                xp = xp[sub[0]:sub[1]]
                yp = yp[sub[2]:sub[3]]
            else:
                print('MAKE_MAP: Zero data points within selected subregion')
                return {}
            if byte_scale:
                #sub_data = bdata[sub[0]:sub[1],sub[2]:sub[3]]
                sub_data = data[sub[0]:sub[1],sub[2]:sub[3]]
            else:
                sub_data = data[sub[0]:sub[1],sub[2]:sub[3]]
        else:
            if np.min(subreg) < 0:
                print('MAKE_MAP: Negative index values -> use data_units=True for such values')
                return {}
            x1 = min(subreg[0],(nx-1)) 
            x2 = min(subreg[1],(nx-1))
            y1 = min(subreg[2],(ny-1))
            y2 = min(subreg[3],(ny-1))
            xp = xp[x1:x2]
            yp = yp[y1:y2]
            if byte_scale:
                #sub_data = bdata[x1:x2,y1:y2]
                sub_data = data[x1:x2,y1:y2]
            else:
                sub_data = data[x1:x2,y1:y2]
        xc = get_arr_center(xp)['center']
        yc = get_arr_center(yp)['center']

    #-- treat roll

    if roll_center is None:
        roll_center = [xc,yc]  

    #-- create final map dictionary 

    map = {'time':stime,'id':id,'dur':dur,'xunits':xunits,'yunits':yunits,
           'roll_angle': np.float64(roll_angle % 360.),
           'roll_center':np.array(roll_center, dtype=np.float64)}

    if type(extra) is dict:
        map.update(extra)
    map.update({'xc':xc, 'yc':yc, 'dx':dx, 'dy':dy})

    if do_sub:
        map.update({'data':sub_data})
    elif not bdata is None:
        map.update({'data':bdata})
    else:
        map.update({'data':data})

    return map
