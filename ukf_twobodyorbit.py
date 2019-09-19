#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import astropy.units as u
import numpy as np
#import matplotlib.pyplot as plt
from filterpy.kalman import UnscentedKalmanFilter, unscented_transform
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
#from poliastro.twobody.propagation import mean_motion
#from astropy.coordinates import EarthLocation, SkyCoord, AltAz, GCRS, ICRS
#from astropy.time import Time, TimeDelta
#from Calc_Residuals import Time_Received, Alter_Time, Residuals


"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#X[0][:,0], X[0][:,1], X[0][:,2]
def stm_update(attractbody, state):
    r_mag = np.linalg.norm(state.r)
    r_mag5 = r_mag*r_mag*r_mag*r_mag*r_mag
#    v_mag = np.linalg.norm(state.v)
    update_pos = -(attractbody.k/(r_mag*r_mag*r_mag))
    lam_ii = (3*attractbody.k*state.r[0]*state.r[0])/(r_mag5)
    lam_jj = (3*attractbody.k*state.r[1]*state.r[1])/(r_mag5)
    lam_kk = (3*attractbody.k*state.r[2]*state.r[2])/(r_mag5)
    lam_ij = (3*attractbody.k*state.r[0]*state.r[1])/(r_mag5)
    lam_ik = (3*attractbody.k*state.r[0]*state.r[2])/(r_mag5)
    lam_jk = (3*attractbody.k*state.r[1]*state.r[2])/(r_mag5)
    lam = np.array([[lam_ii+update_pos, lam_ij, lam_ik],[lam_ij, lam_jj+update_pos, lam_jk], [lam_ik, lam_jk, lam_kk+update_pos]])
    return np.array([[0, np.identity(6)],[lam, 0]])

def stm(attractbody, state):
    r_mag = np.linalg.norm(state.r)
    r_mag5 = r_mag*r_mag*r_mag*r_mag*r_mag
#    v_mag = np.linalg.norm(state.v)
    update_pos = -(attractbody.k/(r_mag*r_mag*r_mag))
    lam_ii = (3*attractbody.k*state.r[0]*state.r[0])/(r_mag5)
    lam_jj = (3*attractbody.k*state.r[1]*state.r[1])/(r_mag5)
    lam_kk = (3*attractbody.k*state.r[2]*state.r[2])/(r_mag5)
    lam_ij = (3*attractbody.k*state.r[0]*state.r[1])/(r_mag5)
    lam_ik = (3*attractbody.k*state.r[0]*state.r[2])/(r_mag5)
    lam_jk = (3*attractbody.k*state.r[1]*state.r[2])/(r_mag5)
    lam = np.array([[lam_ii+update_pos, lam_ij, lam_ik],[lam_ij, lam_jj+update_pos, lam_jk], [lam_ik, lam_jk, lam_kk+update_pos]])
    return lam
