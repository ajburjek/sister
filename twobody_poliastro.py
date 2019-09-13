#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import mean_motion
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
from astropy.time import Time, TimeDelta
from Calc_Residuals import Time_Received, Alter_Time, Residuals


"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#Hardaway Hall
loc_hh = EarthLocation.from_geodetic(-87.545*u.deg,33.213*u.deg,height=19.25*u.m,ellipsoid='WGS84')
#time
t_ob = Time(Time.now(),scale='tai')
delta_t = TimeDelta(np.linspace(0,86400,num=86401)*u.second,scale='tai')
t_del_nounit = np.linspace(0,86400,num=86401)*u.second

t_interval = 10*u.second
#topocentric frame at Hardaway Hall
frame_now_hh = AltAz(obstime=t_ob+delta_t,location=loc_hh)
#ISS
a    = 6799.18 * u.km
ecc  = 0.0008424 * u.one
inc  = 51.6403 * u.deg
raan = 292.1260 * u.deg
aop  = 28.6368 * u.deg
nu   = 331.5244 * u.deg
ss = Orbit.from_classical(Earth, a, ecc, inc, raan, aop, nu, epoch=t_ob)
#coordinates are printed as a tuple (r,v)
X = mean_motion(Earth.k, ss.r, ss.v, t_del_nounit)
ss_SkyCoord = SkyCoord(X[0][:,0], X[0][:,1], X[0][:,2], frame = 'gcrs', representation_type='cartesian')
ssaltaz_hh = ss_SkyCoord.transform_to(frame_now_hh)
for i in range(0,ssaltaz_hh.alt.size):
    if ssaltaz_hh.alt[i]>=0*u.deg:
        time_rec = Time_Received(ssaltaz_hh.distance,t_interval)
        time_alt = Alter_Time(time_rec)
        residual = Residuals(time_rec,time_alt)
