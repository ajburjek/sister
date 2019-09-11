#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import astropy.units as u
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter2D

"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
a    = 6799.18 * u.km
ecc  = 0.0008424 * u.one
inc  = 51.6403 * u.deg
raan = 292.1260 * u.deg
aop  = 28.6368 * u.deg
nu   = 331.5244 * u.deg
ss = Orbit.from_classical(Earth, a, ecc, inc, raan, aop, nu)
ss.plot()
X = np.empty([6,86400])
ss_update = ss.propagate(30*u.second)
for t in range(0,86401):
    ss_update = ss.propagate(float(t)*u.second)
    X[0:6,t] = [[ss_update.r[0]], [ss_update.r[1]], [ss_update.r[2]], [ss_update.v[0]], [ss_update.v[1]], [ss_update.v[2]]]