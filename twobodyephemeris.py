#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import numpy as np
#custom
import rv2oe
import keplerseqn
"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def twobodyephem(X0, t):
    #TWOBODYEPHEM Produce two-body ephemeris table from time and state
    #Using two-body point force, propagate from initial state through given    
    #time to produce ephemeris.
    #INPUTS:
    #   X0 - Given state in inertial frame at time t. [meters, meters/seconds]
    #   t  - State will be predicted for each point. [seconds]
    
    mu = 3.9860044e+14;  #m^3/s^2 Earth gravitational parameter
    we = 2*np.pi()/86164;#rad/sec Earth average rotational rate
    re = 6378137;#meters, spherical Earth radius
    ag = 0;#Greenwich angle
    [oe, T, rp, ra, f0] = rv2oe(X0, mu)
    
    a    = oe[1]#semimajor axis, meters
    e    = oe[2]#eccentricity
    i    = oe[3]#inclination, radians
    Si = np.sin(i) 
    Ci = np.cos(i)
    raan = oe[4]#%ascending node, rad
    Sraan = np.sin(raan)
    Craan = np.cos(raan)
    w    = oe[5]#aop, radians
    Sw = np.sin(w)
    Cw = np.cos(w)
    M0    = oe[6]#mean anomaly, radians
    Sm0 = np.sin(M0)
    Cm0 = np.cos(M0)
    
    n = np.sqrt(mu/a**3)#mean motion
    for i in range(0,np.size(t)):
        Et = keplerseqn(e,n,t[i],M0)
        Mt = Et - (e*np.sin(Et))
        ft = np.arccos((np.cos(Et)-e)/(1-(e*np.cos(Et))))
        oe[6] = Mt
        Xt[:,i] = oe2rv(oe,ft)
        