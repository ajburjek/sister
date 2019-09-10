#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import numpy as np
"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
def rv2oe(rv,mu):
    #RV2OE Convert from position and velocity to orbital elements
    #Given a known gravitational parameter for a planet, an artificial satellit
    #e initial state is provided in a nonrotating coordinate system.
    
    rholder = rv[0:3]
    vholder = rv[3:6]
    del rv
    rv = np.vstack((rholder,vholder))
    
    #magnitudes
    rmag = np.linalg.norm(rholder)#m
    vmag = np.linalg.norm(vholder)#m/s
    
    #specific angular momentum
    h = np.cross(rholder,vholder)#m^2/s
    hmag = np.linalg.norm(h)
    #inclination
    i = np.arccos(h[2]/hmag)#rad
    #raan
    hxy = np.sqrt(h[0]**2 + h[1]**2)
    sinraan = h[0]/hxy
    cosraan = -h[1]/hxy
    raan = np.arctan2(sinraan,cosraan)
    
    #specific energy per unit mass
    energy = ((vmag**2)/2) - (mu/rmag)
    
    #semimajor axis
    a = -mu/(2*energy)
    
    #eccentricity
    e = np.sqrt(1 + ((2*energy*hmag**2)/mu**2))
    evec = np.cross(vholder,(h/mu)) - rholder/rmag
    
    #semiminor axis
    b = a*np.sqrt(1 - e**2)
    
    #semi-latus rectum
    p = a*(1 - e**2)
    
    #true anomaly
    if np.dot(rholder,vholder)<0:
        f = (2*np.pi() - np.arccos(np.dot(evec,rholder)/(e*rmag)))
    else:
        f = np.arccos(np.dot(evec,rholder)/(e*rmag))
    
    #perifocus
    sinwf = rholder[2]/(rmag*np.sin(i))
    coswf = ((rholder[0]/rmag)*np.cos(raan))+((rholder[1]/rmag)*sin(raan))
    wf = np.arctan2(sinwf,coswf)
    w = wf - f
    
    #eccentric anomaly
    cosE0 = ((rmag/a)*np.cos(f)) + e
    sinE0 = ((rmag/b)*np.sin(f))
    E0 = np.arctan2(sinE0,cosE0)
    
    #Mean Anomaly
    M0 = E0 - e*np.sin(E0)
    
    #Period
    T = 2*np.pi()*np.sqrt(a**3/mu)
    rp = a*(1-e)
    ra = a*(1+e)
    
    oe = np.array([a, e, i, raan, w, M0])
    return oe, T, rp, ra, f