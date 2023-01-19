#!/usr/bin/env python

"""utils.py: utility module for data parsing and helping analsysis or plots."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

sys.path.extend(["agw/", "agw/tid/"])

import configparser
import datetime as dt
import glob

import numpy as np
from numpy import ma
import pandas as pd
from loguru import logger


def get_config(key, section="LWS"):
    config = configparser.ConfigParser()
    config.read("config/conf.ini")
    val = config[section][key]
    return val


def read_logs(_filestr="config/record.logs/*.log"):
    """
    Read all log files and store them each records
    to a list of dictionary
    """
    files = sorted(glob.glob(_filestr))
    _modes = pd.DataFrame()
    for f in files:
        logger.info(f"Log file {f}")
        _modes = pd.concat([_modes, pd.read_csv(f, parse_dates=["date"])])
    _modes.ut_start = [
        r["date"] + dt.timedelta(hours=r["ut_start"]) for i, r in _modes.iterrows()
    ]
    _modes.ut_end = [
        r["date"] + dt.timedelta(hours=r["ut_end"]) for i, r in _modes.iterrows()
    ]
    logger.info(f"Total _mode log entry {len(_modes)}")
    return _modes


def get_gridded_parameters(
    q, xparam="beam", yparam="slist", zparam="v", r=0, rounding=True
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    else:
        plotParamDF[xparam] = plotParamDF[xparam].tolist()
        plotParamDF[yparam] = plotParamDF[yparam].tolist()
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(xparam, yparam)
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z

def datetimeToEpoch(my_date):
    """reads in a datetime and outputs the equivalent epoch time
    Parameters
    ----------
    my_date : datetime
        a datetime object
    Returns
    -------
    my_epoch : float or list of floats
        an epoch time equal to the datetime object
    Example
    -------
        import datetime as dt
        my_epoch = utils.timeUtils.datetimeToEpoch(dt.datetime(2012,7,10))
    Written by AJ 20120914
    Modified by Nathaniel Frissell 20130729 - Added list support.
    Modified by AJ 20130925 - fixed local time bug with conversion
    Modified by Nathaniel Frissell 20130930 - Fixed list support.
    Modified by ASR 20151120
    """
    from datetime import datetime
    import calendar
    import numpy as np

    assert(isinstance(my_date, (list, np.ndarray, datetime))), logging.error(
        'input must be of type datetime or list of datetimes')
    if isinstance(my_date, (list,np.ndarray)):
        for dt in my_date:
            assert(isinstance(dt, datetime)), logging.error(
                'each member of my_date list must be of type datetime')

    if isinstance(my_date, datetime):
        unx = calendar.timegm(my_date.timetuple()) + my_date.microsecond / 1e6
    else:
        unx = [calendar.timegm(dt.timetuple()) +
               dt.microsecond / 1e6 for dt in my_date]

    return unx


def getParamDict(param):
    """Get information about a parameter, including units, default ranges,
    and axis labels.
    Parameters
    ----------
    param : str
        name of parameter
    Returns
    -------
    paramDict : str
        dictionary containing information about the chosen parameter
    Example
    -------
        paramDict = getParamDict('w_l')
    written by Nathaniel Frissell, 2013-07
    """
    import numpy as np

    # Create empty dictionary.
    paramDict = {}

    if param == 'p_l' or param == 'power':
        paramDict['param'] = 'power'
        paramDict['label'] = r'$\lambda$ Power'
        paramDict['unit'] = 'dB'
        paramDict['range'] = (0, 30)
    elif param == 'p_s':
        paramDict['param'] = 'power'
        paramDict['label'] = r'$\sigma$ Power'
        paramDict['unit'] = 'dB'
        paramDict['range'] = (0, 30)
    elif param == 'v' or param == 'velocity':
        paramDict['param'] = 'velocity'
        paramDict['label'] = 'Velocity'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (-500, 500)
    elif param.find('vheight') >= 0:
        paramDict['param'] = 'height'
        paramDict['label'] = "h'"
        paramDict['unit'] = 'km'
        paramDict['range'] = (75.0, 900.0)
    elif param == 'w_l' or param == 'width':
        paramDict['param'] = 'width'
        paramDict['label'] = r'$\lambda$ Spectral Width'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (0, 100)
    elif param == 'w_s':
        paramDict['param'] = 'width'
        paramDict['label'] = r'$\sigma$ Spectral Width'
        paramDict['unit'] = 'm s^{-1}'
        paramDict['range'] = (0, 100)
    elif param.find('elv') >= 0:
        paramDict['param'] = 'elevation'
        paramDict['label'] = 'Elevation'
        paramDict['unit'] = 'degrees'
        paramDict['range'] = (10, 30)
    elif param == 'phi0':
        paramDict['param'] = 'phi'
        paramDict['label'] = r'$\phi_0$'
        paramDict['unit'] = 'radians'
        paramDict['range'] = (-np.pi, np.pi)
    return paramDict

def epem(date):
    """
    input: date - datetime object (assumed UTC)
    ouput: gha - Greenwich hour angle, the angle between the Greenwich
           meridian and the meridian containing the subsolar point.
           dec - solar declination.
    """
    dg2rad = np.pi/180.
    rad2dg = 1./dg2rad
    # compute julian day from UTC datetime object.
    # datetime objects use proleptic gregorian calendar.
    jday = JulianDayFromDate(date,calendar='proleptic_gregorian')
    jd = np.floor(jday) # truncate to integer.
    # utc hour.
    ut = date.hour + date.minute/60. + date.second/3600.
    # calculate number of centuries from J2000
    t = (jd + (ut/24.) - 2451545.0) / 36525.
    # mean longitude corrected for aberration
    l = (280.460 + 36000.770 * t) % 360
    # mean anomaly
    g = 357.528 + 35999.050 * t
    # ecliptic longitude
    lm = l + 1.915 * np.sin(g*dg2rad) + 0.020 * np.sin(2*g*dg2rad)
    # obliquity of the ecliptic
    ep = 23.4393 - 0.01300 * t
    # equation of time
    eqtime = -1.915*np.sin(g*dg2rad) - 0.020*np.sin(2*g*dg2rad) \
            + 2.466*np.sin(2*lm*dg2rad) - 0.053*np.sin(4*lm*dg2rad)
    # Greenwich hour angle
    gha = 15*ut - 180 + eqtime
    # declination of sun
    dec = np.arcsin(np.sin(ep*dg2rad) * np.sin(lm*dg2rad)) * rad2dg
    return gha, dec

def JulianDayFromDate(date,calendar='standard'):

    """
creates a Julian Day from a 'datetime-like' object.  Returns the fractional
Julian Day (resolution 1 second).
if calendar='standard' or 'gregorian' (default), Julian day follows Julian
Calendar on and before 1582-10-5, Gregorian calendar after 1582-10-15.
if calendar='proleptic_gregorian', Julian Day follows gregorian calendar.
if calendar='julian', Julian Day follows julian calendar.
Algorithm:
Meeus, Jean (1998) Astronomical Algorithms (2nd Edition). Willmann-Bell,
Virginia. p. 63
    """
    # based on redate.py by David Finlayson.
    year=date.year; month=date.month; day=date.day
    hour=date.hour; minute=date.minute; second=date.second
    # Convert time to fractions of a day
    day = day + hour/24.0 + minute/1440.0 + second/86400.0
    # Start Meeus algorithm (variables are in his notation)
    if (month < 3):
        month = month + 12
        year = year - 1
    A = int(year/100)
    jd = int(365.25 * (year + 4716)) + int(30.6001 * (month + 1)) + \
         day - 1524.5
    # optionally adjust the jd for the switch from
    # the Julian to Gregorian Calendar
    # here assumed to have occurred the day after 1582 October 4
    if calendar in ['standard','gregorian']:
        if jd >= 2299170.5:
            # 1582 October 15 (Gregorian Calendar)
            B = 2 - A + int(A/4)
        elif jd < 2299160.5:
            # 1582 October 5 (Julian Calendar)
            B = 0
        else:
            raise ValueError('impossible date (falls in gap between end of Julian calendar and beginning of Gregorian calendar')
    elif calendar == 'proleptic_gregorian':
        B = 2 - A + int(A/4)
    elif calendar == 'julian':
        B = 0
    else:
        raise ValueError('unknown calendar, must be one of julian,standard,gregorian,proleptic_gregorian, got %s' % calendar)
    # adjust for Julian calendar if necessary
    jd = jd + B
    return jd

def genCmap(param, scale, colors='lasse', lowGray=False):
    """Generates a colormap and returns the necessary components to use it
    Parameters
    ----------
    param : str
        the parameter being plotted ('velocity' and 'phi0' are special cases,
        anything else gets the same color scale)
    scale : list
        a list with the [min,max] values of the color scale
    colors : Optional[str]
        a string indicating which colorbar to use, valid inputs are 
        'lasse', 'aj'.  default = 'lasse'
    lowGray : Optional[boolean]
        a flag indicating whether to plot low velocities (|v| < 15 m/s) in
        gray.  default = False
    Returns
    -------
    cmap : matplotlib.colors.ListedColormap
        the colormap generated.  This then gets passed to the mpl plotting
        function (e.g. scatter, plot, LineCollection, etc.)
    norm : matplotlib.colors.BoundaryNorm
        the colormap index.  This then gets passed to the mpl plotting
        function (e.g. scatter, plot, LineCollection, etc.)
    bounds : list
        the boundaries of each of the colormap segments.  This can be used
        to manually label the colorbar, for example.
    Example
    -------
        cmap,norm,bounds = genCmap('velocity', [-200,200], colors='aj', lowGray=True)
    Written by AJ 20120820
    """
    import matplotlib,numpy
    import matplotlib.pyplot as plot

    #the MPL colormaps we will be using

    cmj = matplotlib.cm.jet
    cmpr = matplotlib.cm.prism

    #check for a velocity plot

    if(param == 'velocity'):
        #check for what color scale we want to use
        if(colors == 'aj'):
            if(not lowGray):
                #define our discrete colorbar
                cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                         cmpr(.11), cmpr(.1),
                                                         cmpr(.175), cmpr(.158),
                                                         cmj(.32), cmj(.37)])
            else:
                cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                         cmpr(.11), cmpr(.1),
                                                         '.6', cmpr(.175),
                                                         cmpr(.158), cmj(.32),
                                                         cmj(.37)])
        else:
            if(not lowGray):
                #define our discrete colorbar
                cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8),
                                                         cmj(.7), cmj(.65),
                                                         cmpr(.142), cmj(.45),
                                                         cmj(.3), cmj(.1)])
            else:
                cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8),
                                                         cmj(.7), cmj(.65),
                                                         '.6', cmpr(.142),
                                                         cmj(.45), cmj(.3),
                                                         cmj(.1)])

        #define the boundaries for color assignments
        bounds = numpy.round(numpy.linspace(scale[0],scale[1],7))
        if(lowGray):
            bounds[3] = -15.
            bounds = numpy.insert(bounds,4,15.)
        bounds = numpy.insert(bounds,0,-50000.)
        bounds = numpy.append(bounds,50000.)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    elif(param == 'phi0'):
        #check for what color scale we want to use
        if(colors == 'aj'):
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmpr(.142), cmpr(.125),
                                                     cmpr(.11), cmpr(.1),
                                                     cmpr(.18), cmpr(.16),
                                                     cmj(.32), cmj(.37)])
        else:
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmj(.9), cmj(.8), cmj(.7),
                                                     cmj(.65), cmpr(.142),
                                                     cmj(.45), cmj(.3),
                                                     cmj(.1)])

        #define the boundaries for color assignments
        bounds = numpy.linspace(scale[0],scale[1],9)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    elif(param == 'grid'):
        #check what color scale we want to use
        if(colors == 'aj'):
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmpr(.175), cmpr(.17),
                                                     cmj(.32), cmj(.37),
                                                     cmpr(.142), cmpr(.13),
                                                     cmpr(.11), cmpr(.10)])
        else:
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmj(.1), cmj(.3), cmj(.45),
                                                     cmpr(.142), cmj(.65),
                                                     cmj(.7), cmj(.8), cmj(.9)])

        #define the boundaries for color assignments
        bounds = numpy.round(numpy.linspace(scale[0],scale[1],8))
        bounds = numpy.append(bounds,50000.)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
    else:
        # If its a non-velocity plot, check what color scale we want to use
        if(colors == 'aj'):
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmpr(.175), cmpr(.158),
                                                     cmj(.32), cmj(.37),
                                                     cmpr(.142), cmpr(.13),
                                                     cmpr(.11), cmpr(.10)])
        else:
            #define our discrete colorbar
            cmap = matplotlib.colors.ListedColormap([cmj(.1), cmj(.3), cmj(.45),
                                                     cmpr(.142), cmj(.65),
                                                     cmj(.7), cmj(.8), cmj(.9)])

        #define the boundaries for color assignments
        bounds = numpy.round(numpy.linspace(scale[0],scale[1],8))
        bounds = numpy.append(bounds,50000.)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    cmap.set_bad('w',1.0)
    cmap.set_over('w',1.0)
    cmap.set_under('.6',1.0)

    return cmap,norm,bounds

def greatCircleAzm(lat1, lon1, lat2, lon2):
    """Calculates the azimuth from the coordinates of a start point to and end
    point along a great circle path.
    Parameters
    ----------
    lat1 : float
        latitude [deg]
    lon1 : float
        longitude [deg]
    lat2 : float
        latitude [deg]
    lon2 : float
        longitude [deg]
    Returns
    -------
    azm : float
        azimuth [deg]
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)
    
    dlon = lon2 - lon1
    y = np.sin(dlon) * np.cos(lat2)
    x = np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)
    azm = np.degrees(np.arctan2(y,x))

    return azm


def greatCircleDist(lat1, lon1, lat2, lon2):
    """Calculates the distance in radians along a great circle path between two
    points.
    Parameters
    ----------
    lat1 : float
        latitude [deg]
    lon1 : float
        longitude [deg]
    lat2 : float
        latitude [deg]
    lon2 : float
        longitude [deg]
    Returns
    -------
    radDist : float
        distance [radians]
    """
    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlat = (lat2 - lat1) / 2.0
    dlon = (lon2 - lon1) / 2.0
    a = np.sin(dlat)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon)**2
    radDist = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    return radDist