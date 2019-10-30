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
#custom
import readtle
from downloadfileurl import download_file
from time_helper import daysToMonDayHrMinSec, jday
from poliastro.core.propagation import mean_motion as mean_motion_fast

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

def mean_motion(k, r, v, tofs, **kwargs):
    
    k = Earth.k
    r0 = Xf_tle_n.r
    v0 = Xf_tle_n.v
    tofs = np.linspace(0.0,0.0,num=1)*u.s
    
    results = [mean_motion_fast(k, r0, v0, tof) for tof in tofs]
    return (
        [result[0] for result in results] * u.km,
        [result[1] for result in results] * u.km / u.s,
    )
    

julianday = np.empty([0])
i=0
for i in range(yd.shape[1]):
    mon, day, hr, minute, sec = daysToMonDayHrMinSec(yd[0,i],yd[1,i])
    julianday = np.hstack([julianday, jday(yd[0,i], mon, day, hr, minute, sec)])
        
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

k = Earth.k
r = Xf_tle_n.r
v = Xf_tle_n.v
tofs = np.linspace(0.0,0.0,num=1)*u.s

Xf_n = mean_motion(k,r,v,tofs)

# fill orbit
i=0
X_res = np.empty((6,1))
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
    tof = np.linspace((epoch_tles[yd.shape[1]-1]-epoch_tles[i]).sec,
                      (epoch_tles[yd.shape[1]-1]-epoch_tles[i]).sec, num=1)*u.s
    Xf_tle_i = mean_motion(Earth.k, X0_tle_i.r, X0_tle_i.v, tof)
    temp_r = np.array([Xf_tle_i[0][0,:]-Xf_n[0][0,:]]).reshape(3,1)
    temp_v = np.array([Xf_tle_i[1][0,:]-Xf_n[1][0,:]]).reshape(3,1)
    X_res = np.hstack([X_res,np.vstack([temp_r,temp_v])])
    del(temp_r)
    del(temp_v)

X_res = np.delete(X_res, 0, axis=1)
mean_X = X_res.sum(axis=1)/(yd.shape[1]-2) #km, km/sec
cov_x = (X_res-mean_X.reshape(6,1)).sum(axis=1).reshape(6,1)*((X_res-mean_X.reshape(6,1)).sum(axis=1))/(yd.shape[1]-3)

    