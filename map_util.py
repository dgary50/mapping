#
# Collection of utilities for working with maps
#
# Based on routines written in IDL by D. Zarro
#
# History:  Written 2020-May-31  DG
#
# Contact     : dgary@njit.edu
#-

from __future__ import print_function
from copy import deepcopy
import numpy as np
import mapping.transform_solar as solar
dtor = np.pi/180.0   # Global conversion of degrees to radians

def default(src, kw, default, pop=True):
    try:
        if pop:
            val = src.pop(kw)
        else:
            val = src[kw]
    except KeyError:
        val = default
    return val
        
def valid_map(map=None):

    if map is None:
        return False

    #-- check for required keys
    if type(map) is dict:
        if not 'data' in map.keys(): return False
        if not 'time' in map.keys(): return False
        if not 'id' in map.keys(): return False
        if not 'xc' in map.keys(): return False
        if not 'yc' in map.keys(): return False
        if not 'dx' in map.keys(): return False
        if not 'dy' in map.keys(): return False
        return True

    #-- LIST object?
    if type(map) is list:
        if valid_map(map[0]):
            return True

    #-- MAP object? (no map objects yet!)
    #if obj_valid(map[0]) then begin
    # if obj_isa(map[0],'map') then begin
    #  pmap=map[0]->get(/map,/pointer)
    #  if ptr_exist(pmap) then begin
    #   chk=valid_map(*pmap,err=err,old_format=old_format)
    #   if chk then type=1
    #   return,chk
    #  endif
    # endif
    #endif

    return False

def get_map_prop(map, prop, default=None):

    if prop in map.keys():
        return map[prop]
    else:
        return default

def get_map_xrange(map, edge=False):

    if not valid_map(map):
        return [0.,0.]

    use_edge = 0.
    if edge:
        use_edge = 1.
    nx, ny = get_data_size(map['data'])
    dx = map['dx']
    xc = map['xc']
    xmin = xc - dx*(nx - 1.)/2. - use_edge*dx/2.
    xmax = xc + dx*(nx - 1.)/2. + use_edge*dx/2.
    return [xmin,xmax]

def get_map_yrange(map, edge=False):

    if not valid_map(map):
        return [0.,0.]

    use_edge = 0.
    if edge:
        use_edge = 1.
    nx, ny = get_data_size(map['data'])
    dy = map['dy']
    yc = map['yc']
    ymin = yc - dy*(ny - 1.)/2. - use_edge*dy/2.
    ymax = yc + dy*(ny - 1.)/2. + use_edge*dy/2.
    return [ymin,ymax]
    
def get_ij(v=None, nx=None):

    if v is None or nx is None:
        print('GET_IJ syntax ij = get_ij(v,nx)')
        return np.array([None,None])

    return np.array([v % nx, np.int(v/nx)])

def get_map_region(xp=None, yp=None, region=None, fast=False):
    ''' Compute pixel indicies of region [xmin,xmax,ymin,ymax] given in data 
        coordinates.
    '''
    badret = [False,False,False,False]
    if len(region) != 4 or xp is None or yp is None:
        print('GET_MAP_REGION syntax index = get_map_region(xp,yp,region)')
        return badret

    xmin, xmax = np.sort(region[:2])
    ymin, ymax = np.sort(region[2:])

    flags = np.logical_and(np.logical_and(np.logical_and(xp >= xmin, xp <= xmax),
                                          yp >= ymin), 
                           yp <= ymax)
    if (np.sum(flags) == 0):
        return badret

    vfind, = np.where(flags.flatten())                 

    #-- fast option

    if fast:
        ox, = np.where(np.logical_and(xp >= xmin, xp <= xmax))
        oy, = np.where(np.logical_and(yp >= ymin, yp <= ymax))
        if len(ox > 0) and len(oy > 0):
            return [np.min(ox),np.max(ox),np.min(oy),np.max(oy)] 
        else:
            return badret

    nx,ny = xp.shape
    
    xy1 = get_ij(vfind[0], nx)
    xy2 = get_ij(vfind[-1], nx)
    imin = [0,0]
    imax = [0,0]
    imin[0] = min([xy1[0],xy2[0]])
    imin[1] = min([xy1[1],xy2[1]])
    imax[0] = max([xy1[0],xy2[0]])
    imax[1] = max([xy1[1],xy2[1]])
    imin[0] = np.clip(imin[0], 0, nx-1)
    imax[0] = np.clip(imax[0], 0, nx-1)
    imin[1] = np.clip(imin[1], 0, ny-1)
    imax[1] = np.clip(imax[1], 0, ny-1)

    # Return in swapped order, so it can be used directly
    return [imin[1],imax[1],imin[0],imax[0]]

def get_map_center(map):

    if not valid_map(map):
        return False
    center = [map.xc, map.yc]

    return center

def get_map_fov(map, arcmin=False, round=False, **kwargs):

    if not valid_map(map):
        print('GET_MAP_FOV: syntax fov=get_map_fov(map,[arcmin=arcmin])')
        return False

    xrange = get_map_xrange(map, **kwargs)
    yrange = get_map_yrange(map, **kwargs)

    fovx = max(xrange) - min(xrange)
    fovy = max(yrange) - min(yrange)

    fov = [fovx, fovy]

    if arcmin: fov = [i/60. for i in fov]
    if round:  fov = [np.float(np.rint(i)) for i in fov]
    if arcmin:
        if fov[0] == 0: fov[0] = 1.
        if fov[1] == 0: fov[1] = 1.
    else:
        if fov[0] == 0: fov[0] = 60.
        if fov[1] == 0: fov[1] = 60.

    return fov

def get_map_angles(map, time=None):

    otime = default(map,'rtime',map['time'], pop=False)
    if time:
        otime = time
    # Calculate the angles for this date in case they are not already in the map
    angles = solar.pb0r(otime)
    # Override with the map values if they exist
    b0 = get_map_prop(map, 'b0', default=angles['b0'])
    l0 = get_map_prop(map, 'l0', default=angles['l0'])
    angles['rsun'] *= 60.
    rsun = get_map_prop(map, 'rsun', default=angles['rsun'])
    angles.update({'b0':b0, 'l0':l0, 'rsun':rsun})
    
    return angles

def get_map_sub_ranges(map, xrange=None, yrange=None, no_overlap=False, 
                            xcor=None, ycor=None):

    reterr = {'arange':None, 'irange':None}

    if not valid_map(map):
        print('GET_MAP_SUB_RANGES: syntax,region=get_map_sub(map,[xrange=xrange,yrange=yrange])')
        return reterr

    arange = 0.
    irange = 0.
    xenter = False
    if not xrange is None:
        dxmin = min(xrange)
        dxmax = max(xrange)
        xenter = dxmin < dxmax

    yenter = False
    if not yrange is None:
        dymin = min(yrange)
        dymax = max(yrange)
        yenter = dymin < dymax

    nx, ny = get_data_size(map['data'])

    #-- extract 1-d coordinate arrays

    if xcor is None:
        xarr = mk_map_xp(map['xc'], map['dx'], nx=nx, ny=1) 
    else:
        nx = xcor.shape
        xmin = np.nanmin(xcor)
        xmax = np.nanmax(xcor)
        xarr = np.linspace(xmin, xmax, nx)

    if ycor is None:
        yarr = mk_map_yp(map['yc'], map['dy'], nx=1, ny=ny) 
    else:
        ny = ycor.shape
        ymin = np.nanmin(ycor)
        ymax = np.nanmax(ycor)
        yarr = np.linspace(ymin, ymax, ny)

    if not xenter: dxmin, dxmax = get_map_xrange(map)
    if not yenter: dymin, dymax = get_map_yrange(map)

    dx = map['dx']
    dy = map['dy']
    dx2 = dx/2.
    dy2 = dy/2.

    if not no_overlap:
      dx2 = -dx2
      dy2 = -dy2 

    if not xenter:
        xstart = 0
        xend = nx - 1 
    else:
        xwhere, = np.where(np.logical_and((xarr + dx2) <= dxmax, (xarr - dx2) >= dxmin))
        if len(xwhere > 0):
            xstart = np.nanmin(xwhere)
            xend = np.nanmax(xwhere) 
        else:
            print('GET_MAP_SUB: Error, No data in specified X-range.')
            return reterr

    if not yenter:
        ystart = 0
        yend = ny - 1 
    else:
        ywhere, = np.where(np.logical_and((yarr + dy2) <= dymax, (yarr - dy2) >= dymin))
        if len(ywhere > 0):
            ystart = np.nanmin(ywhere)
            yend = np.nanmax(ywhere) 
        else:
            print('GET_MAP_SUB: Error, No data in specified X-range.')
            return reterr

    arange = [xarr[xstart],xarr[xend],yarr[ystart],yarr[yend]]
    irange = [ystart, yend, xstart, xend]
    
    return {'arange':arange, 'irange':irange}

def get_map_sub(map, xrange=None, yrange=None, no_overlap=False,
                     xcor=None, ycor=None):

    ranges = get_map_sub_ranges(map, xrange=xrange, yrange=yrange, xcor=xcor, ycor=ycor, 
                                     no_overlap=no_overlap)
    if ranges['irange'] is None:
        return False
        
    xstart, xend, ystart, yend = ranges['irange']
    nx = xend - xstart + 1
    ny = yend - ystart + 1
    if nx < 2 or ny < 2:
        print('GET_MAP_SUB: Error, Extracted data is not 2-D.')

    data = get_sub_data(map['data'], irange=ranges['irange'])

    #-- return coordinate subarrays  (disabled for now...)

    #if arg_present(xp) then begin
    # xp=get_map_xp(map)
    # xp=xp[xstart:xend,ystart:yend]
    #endif

    #if arg_present(yp) then begin
    # yp=get_map_yp(map)
    # yp=yp[xstart:xend,ystart:yend]
    #endif

    return data

def get_data_size(data):

    shape = np.array(data.shape)
    # Check for true-color data and find true-color plane
    if len(shape) == 3:
        true_idx, = np.where(shape <= 4)
        if len(true_idx) != 1:
            print('GET_SUB_DATA: Could not interpret true-color planes')
            return None, None
        else:
            nxy, = np.where(shape > 4)
            nx, ny = shape[nxy]
    else:
        nx, ny = shape
    return ny, nx
    
def get_sub_data(data, irange=None):
    
    if irange is None: return data
    
    ny, nx = get_data_size(data)
    
    if nx is None:
        return data
        
    shape = np.array(data.shape)
    if len(shape) == 3:
        true_idx, = np.where(shape <= 4)
    else:
        true_idx = None
        
    ystart = max(irange[0], 0)
    yend = min(irange[1], ny-1)
    xstart = max(irange[2], 0)
    xend = min(irange[3], nx-1)

    n1 = xend - xstart + 1
    n2 = yend - ystart + 1

    sdata = deepcopy(data[ystart:yend,xstart:xend])
    if (nx == n1) and (ny == n2): return data
    if true_idx:
        if true_idx == 0: sdata = deepcopy(data[:,ystart:yend,xstart:xend])
        if true_idx == 1: sdata = deepcopy(data[ystart:yend,:,xstart:xend])

    return sdata

def get_map_xp(map, oned=False):
    if not valid_map(map):
        return False

    nx, ny = get_data_size(map['data'])
    dx = map['dx']
    xc = map['xc']

    if oned: ny = 1
    xp = mk_map_xp(xc, dx, nx, ny)

    return xp

def get_map_yp(map, oned=False):
    if not valid_map(map):
        return False

    nx, ny = get_data_size(map['data'])
    dy = map['dy']
    yc = map['yc']

    if oned: nx = 1
    yp = mk_map_yp(yc, dy, nx, ny)

    return yp

def get_arr_center(array=None, dx=0., dy=0.):

    sz = len(array.shape)
    if sz == 2:
        ny, nx = array.shape
    else:
        nx, = array.shape 
        ny = 1

    if array is None:
        print('GET_ARR_CENTER syntax: center = get_arr_center(array)')
        return None

    amin = np.nanmin(array)
    amax = np.nanmax(array)

    dx = (amax - amin)/(nx - 1.)
    dy = dx
    if ny > 1:
        dy = (amax - amin)/(ny - 1.)
        xmin = np.nanmin(array)#[0])
        xmax = np.nanmax(array)#[0])
        ymin = np.nanmin(array)#[:,0])
        ymax = np.nanmax(array)#[:,0])
        if xmin == xmax:
            dx = None
        if ymin == ymax:
            dy = None

    center = (amin + amax)/2.

    return {'center':center, 'dx':dx, 'dy':dy}

def repack_map(map=None, xp=None, yp=None, no_copy=False):

    if not valid_map(map):
        print('REPACK_MAP: Syntax, rmap = repack_map(map, xp, yp)')
        if map is None:
            return None
        else:
            return map 
    
    if xp is None or yp is None:
        return None

    if no_copy:
       rmap = map
    else:
       rmap = deepcopy(map)

    xcen = get_arr_center(xp)
    rmap['xc'] = xcen['center']
    rmap['dx'] = xcen['dx']
    ycen = get_arr_center(yp)
    rmap['yc'] = ycen['center']
    rmap['dy'] = ycen['dy']

    return rmap


def mk_map_xp(xc, dx, nx, ny):
    ''' Compute X-coordinate arrays from center and spacing
        Inputs      : xc = x-coord image center (arcsecs)
                      dx = pixel spacing in x-direc (arcsecs)
                      nx,ny = output dimensions

        Outputs     : xp = 2d X-coordinate array;
    '''
    dumx = nx*dx/2.
    xp = (np.arange(nx)+0.5)*dx - dumx + xc
    if ny > 1:
        return np.repeat(xp,ny).reshape(nx,ny).T
    return xp.T

# def mk_map_xp(xc, dx, nx, ny):
    # ''' Compute X-coordinate arrays from center and spacing
        # Inputs      : xc = x-coord image center (arcsecs)
                      # dx = pixel spacing in x-direc (arcsecs)
                      # nx,ny = output dimensions

        # Outputs     : xp = 2d X-coordinate array;
    # '''
    # dumx = nx*dx/2.
    # xp = (np.arange(nx)+0.5)*dx - dumx + xc
    # if ny > 1:
        # return np.repeat(xp,ny).reshape(nx,ny)
    # return xp

def mk_map_yp(yc, dy, nx, ny):
    ''' Compute Y-coordinate arrays from center and spacing
        Inputs      : yc = y-coord image center (arcsecs)
                      dy = pixel spacing in y-direc (arcsecs)
                      nx,ny = output dimensions

        Outputs     : yp = 2d Y-coordinate array;
    '''
    dumy = ny*dy/2.
    yp = (np.arange(ny)+0.5)*dy - dumy + yc
    if nx > 1:
        return np.repeat(yp,nx).reshape(ny,nx)
    else:
        yp.shape = (ny,)
        return yp
     
def roll_xy(xarr,yarr,angle,rcenter=None):

    #if n_elements(angle) eq 0 then begin
    # repeat begin
    #  angle='' & read,'* enter angle [deg] by which to rotate image [+ clockwise]: ',angle
    # endrep until angle ne ''
    # angle=float(angle)
    #endif

    if (angle % 360.) == 0:
        return xarr, yarr

    theta = angle * dtor
    costh = np.cos(theta)
    sinth = np.sin(theta)

    #-- rotate pixel arrays about requested center 
    #-- (if input coords are 2d and center not specified, then use image center) 

    xc = 0.
    yc = 0.
    try:
        if len(rcenter) == 2:
            xc, yc = rcenter
    except TypeError:
        if len(xarr.shape) == 2:
            min_x = np.min(xarr)
            max_x = np.max(xarr)
            min_y = np.min(yarr)
            max_y = np.max(yarr)
            xc = (min_x + max_x)/2.
            yc = (min_y + max_y)/2.
            
    rx = xc + costh*(xarr - xc) + sinth*(yarr - yc)
    ry = yc - sinth*(xarr - xc) + costh*(yarr - yc)

    return rx, ry

def grid_xy(x=None, y=None, gsize=None, gspace=None, preserve_area=False, center=None,
            adjust_resolution=False):
    ''' Create 2-d coordinate grids corresponding to an existing map's grids
        and other inputs to modify the size and/or spacing of the grids.

        Inputs:
            x        The source map's xp grid
            y        The source map's yp grid            
            gsize    If specified, this is the desired x,y size of the updated
                       map, not the j,i size, hence the returned grids will have 
                       the complementary shape to gsize due to Python ordering. If
                       not specified, the returned grids will have the same shape
                       as x, y, unless overridden by preserve_area.
            gspace   If specified, this is the desired x,y scaling of the returned
                       grids.
            center   The output grid's center coordinates.  This just shifts the
                       center without changing the size or resolution.
            preserve_area
                     If True, gsize is ignored and the output grids have an altered
                     spacing (if needed) to match the input grid world coordinates.
            adjust_resolution
                     If True, gspace is ignored and the output grids have an altered
                     size (if needed) to match the input grid resolution

        Returns a dictionary with keys:
            'gx'     The new map's xp grid
            'gy'     The new map's yp grid
            'gsize'  The new map's pixel size in j,i order.  Note that 'gsize' is the
                       complementary shape to any input gsize!
            'gspace' The new map's dx, dy spacing in world units.
    '''
    try:
        if x.shape != y.shape:
            print('GRID_XY: Error, Input arrays do not match in size')
            return None
    except:
        print('GRID_XY: Syntax, outdict = grid_xy(x,y)')
        return None

    min_x = np.nanmin(x)
    max_x = np.nanmax(x)
    min_y = np.nanmin(y)
    max_y = np.nanmax(y)
    xc = (min_x + max_x)/2.
    yc = (min_y + max_y)/2.
    if not center is None:
        xc, yc = center

    xside = abs(max_x - min_x)
    yside = abs(max_y - min_y)
    try:
        ny, nx = x.shape
    except:
        # Arrays are 1-d, apparently
        nx = len(x)
        ny = len(y)
    dx = xside/(nx-1.) 
    dy = yside/(ny-1.)

    mspace = [dx,dy]
    msize = [nx,ny]
    
    # Override natural space and size?
    if gspace: mspace = [gspace[0], gspace[-1]]
    if gsize:  msize  = [gsize[0],  gsize[-1] ]

    # 
    if preserve_area:
        msize = [np.rint(i) for i in [xside/mspace[0] + 1., yside/mspace[1] + 1.]]
    elif adjust_resolution:
        mspace = [xside/(msize[0]-1.), yside/(msize[1]-1.)]
    xg = mk_map_xp(xc=xc, dx=mspace[0], nx=int(msize[0]) ,ny=int(msize[1]))
    yg = mk_map_yp(yc=yc, dy=mspace[1], nx=int(msize[0]) ,ny=int(msize[1]))

    gsize  = [i for i in xg.shape]
    gspace = [mspace[0],mspace[1]]

    return {'gx':xg, 'gy':yg, 'gsize':gsize, 'gspace':gspace}

