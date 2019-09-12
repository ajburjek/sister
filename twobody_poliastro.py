#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import astropy.units as u
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit

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

X = np.empty([6,1])
t = np.linspace(0,86400,num=86401)
for i in range(0,t.size):
    ss_update = ss.propagate(t[i]*u.second)
#    if i==0:
#        X = np.insert(X, 0, [[ss_update.r[0].value], [ss_update.r[1].value], [ss_update.r[2].value], [ss_update.v[0].value], [ss_update.v[1].value], [ss_update.v[2].value]], axis=1)
#    else:
    X = np.append(X, [[ss_update.r[0].value], [ss_update.r[1].value], [ss_update.r[2].value], [ss_update.v[0].value], [ss_update.v[1].value], [ss_update.v[2].value]], axis=1)