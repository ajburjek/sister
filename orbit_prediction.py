#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
from astropy.time import Time
from poliastro.core.angles import D_to_nu
#from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import mean_motion
#custom
import readtle
from downloadfileurl import download_file
from time_helper import daysToMonDayHrMinSec, jday

"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson: christopher.r.simpson@simpsonaerospace.com
   University of Alabama
   Andrew Burjek:          ajburjek@crimson.ua.edu
   William Patton:         wjpatton1@crimson.ua.edu"""
#------------------------------------------------------------------------------
#Download TLEs and find classical orbital elements
rooturl = 'http://www.celestrak.com/NORAD/elements/supplemental/'
iss_tles = download_file(rooturl+'iss.txt')
coe, yd = readtle.readtle(iss_tles)

julianday = np.empty([1,yd.size])
for i in range(yd.shape[1]):
    mon, day, hr, minute, sec = daysToMonDayHrMinSec(yd[0,i],yd[1,i])
    julianday = np.insert(julianday, i, jday(yd[0,i], mon, day, hr, minute, sec))
        
epoch_tles = Time(julianday, format='jd', scale='ut1')

#Estimate mean state and state covariance from TLEs
nu = D_to_nu(coe[5,yd.shape[1]-1])
Xf_tle_n = Orbit.from_classical(Earth,
                                coe[0,yd.shape[1]-1]*u.m,
                                coe[1,yd.shape[1]-1]*u.one,
                                coe[2,yd.shape[1]-1]*u.deg,
                                coe[3,yd.shape[1]-1]*u.deg,
                                coe[4,yd.shape[1]-1]*u.deg, 
                                nu*u.deg,
                                epoch = epoch_tles[yd.shape[1]-1])

Xf_n = mean_motion(Earth.k, Xf_tle_n.r, Xf_tle_n.v, np.linspace(0,0,num=2)*u.s)
# fill orbit
i=0
X_res = np.empty([6,yd.shape[1]-1])
for i in range(yd.shape[1]-1):
    nu = D_to_nu(coe[5,i])
    X0_tle_i = Orbit.from_classical(Earth,
                                    coe[0,i]*u.m,
                                    coe[1,i]*u.one,
                                    coe[2,i]*u.deg,
                                    coe[3,i]*u.deg,
                                    coe[4,i]*u.deg, 
                                    nu*u.deg,
                                    epoch = epoch_tles[i])
    tof = np.linspace((epoch_tles[yd.shape[1]-1]-epoch_tles[i]).sec,(epoch_tles[yd.shape[1]-1]-epoch_tles[i]).sec, num=2)*u.s
    Xf_tle_i = mean_motion(Earth.k, X0_tle_i.r, X0_tle_i.v, tof)
    temp = np.array([[Xf_tle_i[0][0,:]-Xf_n[0][0,:]],[Xf_tle_i[1][0,:]-Xf_n[1][0,:]]]).reshape(6,0)
    X_res = np.insert(X_res, i, temp)

mean_X = X_res/(yd.shape[1]-2)
    