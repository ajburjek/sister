#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import astropy.units as u
import numpy as np
from datetime import datetime
"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def daysToMonDayHrMinSec(year, days):
    """DAYSTOMONDAYHRMINSEC Convert days to month day, hour, minute, second
    See Vallado for copy pasta, www.celestrak.com"""
    
    lmonth = [31,28,31,30,31,30,31,31,30,31,30,31] #last day of month
    dayofyear = int(days//1.0)
    #  ----------------- find month and day of month ----------------
    if(year % 4) == 0:
        lmonth = (31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
    i = 1
    inttemp = 0
    while ((dayofyear> (inttemp + lmonth[i-1])) and (i<12)):
        inttemp = inttemp + lmonth[i-1]
        i += 1
        
    mon = i
    day = dayofyear - inttemp
    
    #  ----------------- find hours minutes and seconds -------------
    temp = (days - dayofyear) * 24.0
    hr   = int(temp // 1.0)
    temp = (temp - hr) * 60.0
    minute  = int(temp // 1.0)
    sec  = (temp - minute) * 60.0

    return mon, day, hr, minute, sec

def jday(year, mon, day, hr, minute, sec):
    """JDAY Convert year, month, day, hour, minute, second to Julian day
    See Vallado for copy pasta, www.celestrak.com, 
    Julian date is elapsed date from noon on Jan 1, 4713 B.C."""
    return (367.0 * year -
          7.0 * (year + ((mon + 9.0) // 12.0)) * 0.25 // 1.0 +
          275.0 * mon // 9.0 +
          day + 1721013.5 +
          ((sec / 60.0 + minute) / 60.0 + hr) / 24.0  #  ut in days
          #  - 0.5*sgn(100.0*year + mon - 190002.5) + 0.5;
          )
    