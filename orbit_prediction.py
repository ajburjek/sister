#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Predict Orbit Based On Tx Beacon at Known/Fixed Interval"""
import numpy as np
import astropy.units as u
from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.twobody.propagation import mean_motion
#custom
import readtle
from downloadfileurl import download_file

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
time = yd[1,:]*u.h

#Estimate mean state and state covariance from TLEs
# fill orbit

#iss_tle = Orbit.from_classical(Earth,
#                               coe[0,:]*u.m,
#                               coe[1,:]*u.one,
#                               coe[2,:]*u.deg,
#                               coe[3,:]*u.deg,
#                               coe[4,:]*u.deg,
#                               coe[5,:]*u.deg)