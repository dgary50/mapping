#+
# Project     : SOHO-CDS
#
# Name        : FITS2MAP
#
# Purpose     : Make an image map from a FITS file
#
# Category    : imaging
#
# Syntax      : fits2map,file,map
#
# Inputs      : FILE = FITS file name (or FITS data + HEADER)
#
# Outputs     : MAP = map structure
#
# Keywords    : HEADER = FITS header (output of last file read)
#               OBJECT = return map as an object
#               INDEX = FITS index matching HEADER
#               LIST = return multiple maps in LIST
#
# History     : 2020-May-26  DG
#                 Inspired by Dominic Zarro's IDL code of the same name,
#                 but this version only works for a small range of fits
#                 file types.
#               2020-Jun-07  DG
#                 Added check for file in url_copy(), and skips download
#                 if it exists.
#               2021-Jun-09  Jiajia Liu
#                 Added two new parameters sample and provider to vso_search()
#
# Contact     : dgary@njit.edu
#-

from __future__ import print_function
from astropy.io import fits
from astropy.time import Time
import os
from make_map import make_map
from map_util import *

def vso_search(trange=None, instrument=None, wavelength=None, sample=None,
               provider=None, source=None):
    try:
        from sunpy.net import Fido, attrs as a
    except:
        print('VSO_SEARCH: Error, requires sunpy. Sunpy not installed?')
        return None
    import astropy.units as u
    if not trange or not instrument:
        print('VSO_SEARCH: Error, must specify a valid time range and instrument')
        return None

    try:
        times = a.Time(Time(trange[0]).iso,Time(trange[1]).iso)
    except:
        print('VSO_SEARCH: Error, invalid time range')
        return None
    inst = a.Instrument(instrument)
    paras = [times, inst]
    if wavelength:
        wave = a.Wavelength(wavelength*u.angstrom)
        paras = paras + [wave]
    if sample is not None:
        sample = a.Sample(sample * u.second) # in units of seconds
        paras = paras + [sample]
    if provider is not None:
        provider = a.Provider(provider)
        paras = paras + [provider]
    if source is not None:
        source = a.Source(source)
        paras = paras + [source]

    qr = Fido.search(*paras)

    return qr

def vso_get(qr, outpath=None, names_only=True):
    if outpath is None:
        outpath = './'
    try:
        from sunpy.net import Fido
        files = Fido.fetch(qr, path=outpath)
    except:
        print('VSO_GET: Error, download failed.')
        return None
    if names_only:
        return files.data
    return files

def is_url(url, scheme=False):

    try:
        from urllib.parse import urlparse
    except:
        from urlparse import urlparse
    res = urlparse(url)
    if scheme:
        if res.scheme == '':
            return False
    if res.netloc == '':
        return False
    return True

def get_aia_daily(time=None, wave=None, out_path='./', map=True, **kwargs):

    import os

    wave_list = ['1600', '1700', '4500', '94', '131', '171', '193', '211', '304', '335', 'blos']

    if time is None:
        time = Time.now()
    try:
        t = Time(time).iso
        datstr = t[:10].replace('-','/')
    except:
        print('GET_AIA_DAILY: Error, Could not interpret time.')
        return None
    if str(wave) in wave_list:
        if wave == 'blos':
            fname = 'fblos.fits'
        else:
            fname = 'f%4.4d' % int(wave) +'.fits'
        #if os.path.exists(out_path+fname):
        #    os.remove(out_path+fname)
        url = 'http://suntoday.lmsal.com/sdomedia/SunInTime/'+datstr+'/'+fname
        replace = False
        if 'replace' in kwargs.keys():
            replace = kwargs['replace']
        filename = url_copy(url, out_path=out_path, replace=replace)
        if map:
            aiamap = fits2map(filename)
            return aiamap
        else:
            return filename
    else:
        print('GET_AIA_DAILY: Error, Could not interpret wavelength')
        return None

def url_copy(url, out_path='./', replace=True):
    ''' Copy a file from the given url to the given path.  If file
        already exists, download is skipped.
    '''
    import requests
    import os

    if out_path[-1] != '/':
        out_path += '/'
    # See if the file already exists
    filename = out_path+os.path.basename(url)
    print('URL_COPY: Downloading file',url,'to',out_path+os.path.basename(url))
    if os.path.exists(filename):
        if not replace:
            print('URL_COPY: Found file on disk.  Will skip download.')
            return filename
    try:
        r = requests.get(url, allow_redirects=True)
        filename = os.path.basename(url)
        try:
            if not os.path.exists(out_path):
                os.makedirs(out_path)
        except:
            print('URL_COPY: Error creating',out_path)
            return None
        open(out_path+filename, 'wb').write(r.content)
    except:
        print('URL_COPY: Error reading from',url)
        return None
    return out_path+filename

def fits2map(file=None, ext=0, header=False, **kwargs):

    #-- check inputs

    if file is None:
        print('FITS2MAP syntax: map = fits2map(file[,header=True])')
        return None
    if is_url(file):
        replace = False
        if 'replace' in kwargs.keys():
            replace = kwargs['replace']
        outfile = url_copy(file, replace=replace)
        if outfile:
            file = outfile
    if not os.path.exists(file):
        print('FITS2MAP Error: File not found',file)
        return None

    #-- return map structure

    hdu_list = fits.open(file)
    headr = hdu_list[ext].header
    if headr['NAXIS'] == 0:
        # Try extension 1
        print('FITS2MAP Warning: Extension',ext,'is not a 2-D image.')
        ext = 1
        print('Trying Extension',ext)
        headr = hdu_list[ext].header
    try:
        data = hdu_list[ext].data
    except:
        # Try fixing any header problem
        hdu_list[ext].verify('silentfix')
        data = hdu_list[ext].data
    if headr['NAXIS'] != 2:
        print('FITS2MAP Error: Extension',ext,'is not a 2-D image.')
        return None
    hdu_list.close(file)
    nx = default(headr, 'NAXIS1', 1024)
    ny = default(headr, 'NAXIS2', 1024)
    dx = default(headr, 'CDELT1', 1.0)
    dy = default(headr, 'CDELT2', 1.0)
    xc = default(headr, 'CRVAL1', 0.0)
    yc = default(headr, 'CRVAL2', 0.0)
    xc += dx*((nx+1)/2. - default(headr,'CRPIX1', 0.0))
    yc += dy*((ny+1)/2. - default(headr,'CRPIX2', 0.0))
    time = default(headr, 'DATE-OBS', Time.now())
    xunits = default(headr, 'CUNIT1', 'arcsec')
    yunits = default(headr, 'CUNIT2', 'arcsec')
    dur = default(headr, 'EXPTIME', 0.0)
    sub = default(kwargs, 'sub', None)
    fov = default(kwargs, 'fov', None)
    data_units = default(kwargs, 'data_units', None)
    roll_angle = default(kwargs, 'roll_angle', 0.0)
    byte_scale = default(kwargs, 'byte_scale', 0.0)
    id = default(headr, 'TELESCOP', '')
    wavelnth = str(default(headr, 'WAVELNTH', '')) + ' ' + default(headr, 'WAVEUNIT', '')
    if wavelnth != ' ':
        id += ' ' + wavelnth.replace('angstrom','A')
    frequency = default(headr, 'RESTFRQ',0)
    if frequency != 0:
        funit = default(headr, 'CUNIT3', 'Hz').strip()
        if funit == 'Hz':
            frequency = frequency/1e9
            funit = 'GHz'
        id += ' %5.2f ' % frequency + funit



    map = make_map(data, xcen=xc, ycen=yc, dx=dx, dy=dy, time=time,
             sub=sub, roll_angle=roll_angle, id=id, data_units=data_units,
             roll_center=None, fov=fov, xunits=xunits, extra=None,
             yunits=yunits, dur=dur, no_copy=True, byte_scale=byte_scale,
             xrange=None, yrange=None)
    if header:
        return [map, headr]
    else:
        return map
