from __future__ import print_function
from astropy.time import Time
import numpy as np
dtor = np.pi/180.
wcs_au   = 1.49597870691e11  # Astronomical unit in m
wcs_rsun = 6.95508e8 # Solar radius in m

def sun_pos(dd):
    ''' Inputs  : date - fractional number of days since JD 2415020.0 
        Returns : a dictionary with:
                 longmed    -  Longitude of sun for mean equinox of date (degs)
                 ra         -  Apparent RA for true equinox of date (degs)
                 dec        -  Apparent declination for true equinox of date (degs)
                 l          -  Apparent longitude (degs)
                 oblt       -  True obliquity (degs)
    '''
    #
    #  This routine is a truncated version of Newcomb's Sun and
    #  is designed to give apparent angular coordinates (T.E.D) to a
    #  precision of one second of time

    #  time in Julian centuries from 1900.0
    t = dd/36525.0

    #  form sun's mean longitude
    l = (279.696678+((36000.768925*t) % 360.0))*3600.0

    #  allow for ellipticity of the orbit (equation of centre)
    #  using the Earth's mean anomoly ME
    me = 358.475844 + ((35999.049750*t) % 360.0)
    ellcor  = (6910.1 - 17.2*t)*np.sin(me*dtor) + 72.3*np.sin(2.0*me*dtor)
    l = l + ellcor

    # allow for the Venus perturbations using the mean anomaly of Venus MV
    mv = 212.603219 + ((58517.803875*t) % 360.0) 
    vencorr = (4.8 * np.cos((299.1017 + mv - me)*dtor) + 
               5.5 * np.cos((148.3133 +  2.0 * mv  -  2.0 * me )*dtor) + 
               2.5 * np.cos((315.9433 +  2.0 * mv  -  3.0 * me )*dtor) + 
               1.6 * np.cos((345.2533 +  3.0 * mv  -  4.0 * me )*dtor) + 
               1.0 * np.cos((318.15   +  3.0 * mv  -  5.0 * me )*dtor))
    l += vencorr

    #  Allow for the Mars perturbations using the mean anomaly of Mars MM
    mm = 319.529425  +  (( 19139.858500 * t)  %  360.0 )
    marscorr = (2.0 * np.cos((343.8883 -  2.0 * mm  +  2.0 * me)*dtor ) +
                1.8 * np.cos((200.4017 -  2.0 * mm  + me) * dtor))
    l += marscorr

    # Allow for the Jupiter perturbations using the mean anomaly of Jupiter MJ
    mj = 225.328328  +  (( 3034.6920239 * t)  %  360.0 )
    jupcorr = (7.2 * np.cos((179.5317 - mj + me )*dtor) +
               2.6 * np.cos((263.2167  -  mj ) *dtor) +
               2.7 * np.cos(( 87.1450  -  2.0 * mj  +  2.0 * me ) * dtor) +
               1.6 * np.cos((109.4933  -  2.0 * mj  +  me ) * dtor))
    l += jupcorr

    # Allow for the Moons perturbations using the mean elongation of
    # the Moon from the Sun D
    d = 350.7376814  + (( 445267.11422 * t)  %  360.0 )
    mooncorr  = 6.5 * np.sin(d*dtor)
    l += mooncorr

    # Allow for long period terms
    longterm  = + 6.4 * np.sin(( 231.19  +  20.20 * t )*dtor)
    l  += longterm
    l  =  ( l + 2592000.0)  %  1296000.0 
    longmed = l/3600.0

    # Allow for Aberration
    l  -=  20.5

    # Allow for Nutation using the longitude of the Moons mean node OMEGA
    omega = 259.183275 - (( 1934.142008 * t ) % 360.0 )
    l  -=  17.2 * np.sin(omega*dtor)

    # Form the True Obliquity
    oblt  = 23.452294 - 0.0130125*t + (9.2*np.cos(omega*dtor))/3600.0

    # Form Right Ascension and Declination
    l = l/3600.0
    ra  = np.arctan2( np.sin(l*dtor) * np.cos(oblt*dtor) , np.cos(l*dtor) ) / dtor

    if (ra < 0.0):
        ra += 360.0

    dec = np.arcsin(np.sin(l*dtor) * np.sin(oblt*dtor)) / dtor
    return {'longmed':longmed, 'ra':ra, 'dec':dec, 'l':l, 'oblt':oblt}

def pb0r(date=None, soho=False, arcsec=False, stereo=False, roll_angle=0.):

    defret = {'p':0.,'b0':0.,'rsun':16.,'l0':0.}
    if arcsec:
        defret['rsun'] *= 60.
    if stereo and soho:
        print('PB0R: Error, Cannot set stereo=True and soho=True simultaneously')
        return defret

    #-- STEREO by-pass

    if stereo:
       #return,pb0r_stereo(date,arcsec=arcsec,l0=l0,error=error, stereo=stereo,roll_angle=roll_angle,_extra=extra)
       print('PB0R: Error, cannot handle STEREO data (yet)')
       return defret

    #  common pb0r, prev_output, prev_soho, prev_utc,prev_l0, sd_const
    #  if exist(prev_output) then output = prev_output 
    #  if exist(prev_l0) then l0=prev_l0
     
    #---------------------------------------------------------------------------
    #  date supplied?
    #---------------------------------------------------------------------------

    if date:
        try:
            date_ = date.replace('/','-')
            time = Time(date_)
        except ValueError:
            print('PB0R: Error, could not interpret given date',date_,'as ISO time.')
            return defret
    else:
        time = Time.now()

    # Number of Julian days since 2415020.0

    de = time.jd - 2415020

    #---------------------------------------------------------------------------
    #  get the longitude of the sun etc.
    #---------------------------------------------------------------------------
    sp = sun_pos(de)
    longmed = sp['longmed']
    ra = sp['ra']
    dec = sp['dec'] 
    appl = sp['l']
    oblt = sp['oblt']

    #---------------------------------------------------------------------------
    #  form aberrated longitude
    #---------------------------------------------------------------------------
    lamda = longmed - (20.5/3600.0)

    #---------------------------------------------------------------------------
    #  form longitude of ascending node of sun's equator on ecliptic
    #---------------------------------------------------------------------------
    node = 73.666666 + (50.25/3600.0)*( (de/365.25) + 50.0 )
    arg = lamda - node

    #---------------------------------------------------------------------------
    #  calculate P, the position angle of the pole
    #---------------------------------------------------------------------------
    p = ((np.arctan(-np.tan(oblt*dtor) * np.cos(appl*dtor)) + 
            np.arctan( -0.12722 * np.cos(arg*dtor))) / dtor)

    #---------------------------------------------------------------------------
    #  ... and B0 the tilt of the axis
    #---------------------------------------------------------------------------
    b = np.arcsin( 0.12620 * np.sin(arg*dtor) ) / dtor

    #---------------------------------------------------------------------------
    #  ... and the semi-diameter
    #
    #
    #  Form the mean anomalies of Venus(MV),Earth(ME),Mars(MM),Jupiter(MJ)
    #  and the mean elongation of the Moon from the Sun(D).
    #
    #---------------------------------------------------------------------------
    t = de/36525.0

    mv = 212.6   + ( (58517.80   * t) % 360.0 )
    me = 358.476 + ( (35999.0498 * t) % 360.0 )
    mm = 319.5   + ( (19139.86   * t) % 360.0 )
    mj = 225.3   + ( ( 3034.69   * t) % 360.0 )
    d = 350.7    + ( (445267.11  * t) % 360.0 )

    #---------------------------------------------------------------------------
    #  Form the geocentric distance(r) and semi-diameter(sd)
    #---------------------------------------------------------------------------
    r = (1.000141 - (0.016748 - 0.0000418*t)*np.cos(me*dtor) 
          - 0.000140 * np.cos(2.0*me*dtor)                       
          + 0.000016 * np.cos((58.3 + 2.0*mv - 2.0*me)*dtor) 
          + 0.000005 * np.cos((209.1 + mv - me)*dtor)            
          + 0.000005 * np.cos((253.8 - 2.0*mm + 2.0*me)*dtor)
          + 0.000016 * np.cos(( 89.5 - mj + me)*dtor)            
          + 0.000009 * np.cos((357.1 - 2.0*mj + 2.0*me)*dtor) 
          + 0.000031 * np.cos(d*dtor))

    ## sd = (0.2665685/r)*60.0

    sd_const = wcs_rsun / wcs_au
    sd = np.arcsin(sd_const/r)*10800./np.pi
    if soho:
        soho_sd = sd*1.01
        l0=0.
        output = {'p':p, 'b0':b, 'rsun':soho_sd, 'l0':l0}
    else:
        l0=0.
        output = {'p':p, 'b0':b, 'rsun':sd, 'l0':l0}
     
    if arcsec:
        output['rsun'] *= 60.
    return output

def arcmin2hel(xx_in, yy_in, date=None, p=None, b0=None, r0=None, l0=None, rsun=None, 
                sphere=False, backside=False, vdict=None):

    from copy import deepcopy
    xx = deepcopy(xx_in)
    yy = deepcopy(yy_in)

    try:
        test = len(xx)  # Triggers an erro if xx is a scalar
        xx = np.array(xx).flatten()
        npts = len(xx)
        yy = np.array(yy).flatten()
        if len(yy) != npts:
            print('ARCMIN2HEL: Error, the two input parameters are not of same length')
            return None
    except:
        npts = 1
        xx = np.array([xx])
        yy = np.array([yy])

    offlimb = np.zeros_like(xx, dtype=bool)

#-- check for overriding keywords

    sun_angles = {'p':p, 'b0':b0, 'rsun':rsun, 'l0':l0}
    need_b0  = b0 is None
    need_rad = r0 is None and rsun is None
    need_l0  = l0 is None
    if need_rad or need_b0 or need_l0:
       sun_angles = pb0r(date)

#-- allow keywords to override individual angles

    if p: sun_angles['p'] = p
    if b0: sun_angles['b0'] = b0
    if r0: sun_angles['rsun'] = np.arctan(1./r0)*60./dtor # radians-to-arcmin conversion
    if rsun: sun_angles['rsun'] = rsun/60.
    if l0: sun_angles['l0'] = l0
    
#-- override with vdict, if supplied
    if type(vdict) is dict: sun_angles.update(vdict)

    b0_r = sun_angles['b0'] * dtor
    radius = sun_angles['rsun']
    robs = 1./np.tan(sun_angles['rsun']*dtor/60.)

    xxt = np.tan(xx*dtor/60.) #(Convert to radians & tanify)
    yyt = np.tan(yy*dtor/60.) #(Convert to radians & tanify)

# Convert to cylindrical angular coordinates and azimuth -- makes
# the final transformation easier.  Here, ra is the angle out from
# centerline; phi is the azimuth.  This reduces the problem to 2-D
# geometry in the observer -- Sun-center -- viewpoint plane.

# Load phi with atan(xxt,yyt)

    rat2 = xxt**2 + yyt**2
    phi = np.zeros_like(rat2)
    w_rat2 = np.where(rat2 != 0)
    if len(w_rat2[0]) > 0:
        phi[w_rat2] = np.arctan2(xxt[w_rat2],yyt[w_rat2])

    max_ra = np.arcsin(1./robs)
    max_rat2 = np.tan(max_ra)**2

    outside = rat2 > max_rat2
    if np.sum(outside) > 0:
        rat2[outside] = max_rat2
        offlimb[outside] = True

# Solving for the intersection of the line of sight with the sphere
# gives a z-coordinate (toward the observer) of
#   z = R * (sin(ra))^2 +/- sqrt( Ro^2 - (sin(ra))^2 * R^2 )
# with Ro = the solar radius, ra the angular displacement from disk 
# center, and R the viewpoint distance from Sun center.
#
# We normally want the positive branch, which represents the front
# side of the Sun; but who knows? Someone may want the opposite.

    ras2 = np.zeros_like(rat2)
    if len(w_rat2) > 0:
        ras2[w_rat2] = 1./(1. + 1./rat2[w_rat2])
    d1 = np.clip(1. - ras2, 0., None)
    d2 = np.clip(1. - (robs**2)*ras2, 0., None)
    if not backside:
        x = ras2*robs + np.sqrt(d1)*np.sqrt(d2)
    else:  # This branch is for the far side of the sun
        x = ras2*robs - np.sqrt(d1)*np.sqrt(d2) 

    rr = np.sqrt(np.clip(rat2,0.,None)) * (robs - x) 

# Now we can finally convert back to xyz coords and do the 
# helioraphic conversion.  x: towards obs., y: west, z: North

    t1 = np.sin(phi)*rr
    t2 = np.cos(phi)*rr
    xyz = np.array([[x], [t1], [t2]]).reshape(3,npts)
#---------------------------------------------------------------------------
#  rotate around y axis to correct for B0 angle (B0: hel. lat. of diskcenter)
#---------------------------------------------------------------------------
    rotmx = np.array([[np.cos(b0_r), 0., -np.sin(b0_r)], 
                      [          0., 1.,           0.], 
                      [np.sin(b0_r), 0.,  np.cos(b0_r)]])
    xyz = rotmx.dot(xyz)
#---------------------------------------------------------------------------
#  calc. latitude and longitude.
#---------------------------------------------------------------------------

    latitude = np.arcsin(xyz[2])
    latitude = np.clip(latitude, -89.99*dtor, 89.99*dtor) # force lat. between -pi/2 and pi/2
    longitude = np.arctan2(xyz[1], xyz[0]) # longitude
#---------------------------------------------------------------------------
#  longitude may be larger than 90 degrees due to nonzero B0: get proper value
#---------------------------------------------------------------------------
    if not sphere:
        ii = xyz[0] < 0.0       # where values are negative, ii is True
        for i,truth in enumerate(ii):
            if truth:
                if xyz[0,i] >= 0:
                    longitude[i] = np.pi - longitude[i]
                else:
                    longitude[i] = -np.pi - longitude[i]

#---------------------------------------------------------------------------
#  convert to degrees.  If we have an L0 offset, make sure the output is 
#  between +/- 180 degrees!
#---------------------------------------------------------------------------

    if l0 and l0 != 0:
        longitude = (((5.*np.pi + longitude + (l0 % 360.)*dtor) % 2.*np.pi) - np.pi)/dtor

    return {'lat':latitude/dtor, 'lon':longitude/dtor, 'offlimb':offlimb, 'angles':sun_angles} 

def hel2arcmin(ns=None, ew=None, date=None, arcsec=False,
               b0=None, l0=None, p=None, r0=None, soho=False,
               rsun=None, vdict=None):
    ''' Compute position relative to Sun center from heliographic coordinates
        Inputs    : ns  -   array of Heliographic latitudes in degrees (can be a
                            single string with N/S first character instead of sign).
                    ew  -   array of Heliographic longitudes in degrees (can be a
                            single string with E/W first character instead of sign).
        Returns   : xy      2 x npts array of heliocentric coordinates in 
                            arcmin or arcsec if arcsec=True
                  : visible npts boolean array indicating which of the coordinates
                            are on the visible side of the Sun (True) or not (False)
        optionally: sun_angles dictionary if angles=True, with keys: 
                    'p'   : p-angle(deg), 
                    'b0'  : b-angle (deg), 
                    'rsun': solar radius (arcmin) (or arcsec if arcsec=True), 
                    'l0'  : apparent longitude
        Parameters: 
                    date  -  The date (iso format string) to use in the calculation
          
            The following parameters are not normally needed, but are provided
            to override the values returned from the pb0r() call
                    b0    -  The b-angle, in degrees
                    p     -  The p-angle, in degrees
                    r0    -  The distance of the observer from the Sun, in solar radii
                    l0    -  The longitude of the observer, relative to Earth, in degrees
                    rsun  -  solar radius (input in arcsecs), overrides r0
                    vdict -  optional sun_angles dictionary with keys defined above, instead
                             of individual angles
    '''
    if type(ns) is str:
        try:
            n = np.float64(ns[1:])
            if ns[0].upper() == 'S':
                n = -n
            w = np.float64(ew[1:])
            if ew[0].upper() == 'E':
                w = -w
        except:
            print('HEL2ARCMIN: Error, could not interpret inputs ',ns,ew)
            return [None, None]
    else:
        n = ns
        w = ew
        try:
           test = len(n)
        except TypeError:
           n = np.array(n)
        try:
           test = len(w)
        except TypeError:
           w = np.array(w)
  
    #-- check for overriding keywords

    sun_angles = {'p':p, 'b0':b0, 'rsun':rsun, 'l0':l0}
    need_b0  = b0 is None
    need_rad = r0 is None and rsun is None
    need_l0  = l0 is None
    if need_rad or need_b0 or need_l0:
       sun_angles = pb0r(date, soho=soho, arcsec=arcsec)

#-- allow keywords to override individual angles

    if p: sun_angles['p'] = p
    if b0: sun_angles['b0'] = b0
    if r0: sun_angles['rsun'] = np.arctan(1./r0)*60./dtor # radians-to-arcmin conversion
    if rsun: sun_angles['rsun'] = rsun/60.
    if l0: sun_angles['l0'] = l0
    
#-- override with vdict, if supplied
    if type(vdict) is dict: sun_angles.update(vdict)

#-- convert to radians, and use L0 if present
    lon = (w - sun_angles['l0'])*dtor

# vect is the (x,y,z) location of the point for b0 = 0, where x is in the
# direction of Texas, y is west, and z is north. vect1 is rotated by b0. 
    radius = sun_angles['rsun']                       
    b0_r = sun_angles['b0'] * dtor
    colat = (90. - n) * dtor
    robs = 1./np.tan(radius/60. * dtor)
    maxdim = max([len(colat.flatten()), len(lon.flatten()), 1])
    answer = np.zeros((2, maxdim), dtype=np.float64)
    xcoord = np.zeros(maxdim, dtype=np.float64)

#  calculate the result
    scl = np.sin(colat)
    ccl = np.cos(colat)
    cb0 = np.cos(b0_r)
    sb0 = np.sin(b0_r)
    sl = np.sin(lon)
    cl = np.cos(lon)
    
    xcoord = sb0 * ccl + cb0 * cl * scl
    answer[0] = (np.arctan2( scl*sl, robs - xcoord) / dtor * 60.).flatten()
    answer[1] = (np.arctan2( -scl*cl*sb0 + ccl*cb0, robs - xcoord) / dtor * 60.).flatten()
    visible = xcoord > 0
  
    return {'xx':answer[0], 'yy':answer[1], 'visible':visible, 'angles':sun_angles}

