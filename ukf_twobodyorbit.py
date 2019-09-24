#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Two-Body Orbit"""
import numpy as np
import astropy.units as u
from twobody_poliastro import hh_obs
from poliastro.bodies import Earth
#from filterpy.kalman import UnscentedKalmanFilter, unscented_transform

"""Python 3.7
   Simpson Aerospace (c) 2019
   Christopher R. Simpson
   christopher.r.simpson@simpsonaerospace.com"""
#------------------------------------------------------------------------------
#X[0][:,0], X[0][:,1], X[0][:,2]
def stm_update(state):
    r_mag = np.sqrt(state.r[0]**2 + state.r[1]**2 + state.r[2]**2) 
    r_mag5 = r_mag*r_mag*r_mag*r_mag*r_mag
#    v_mag = np.linalg.norm(state.v)
    update_pos = -(state.attractor.k/(r_mag*r_mag*r_mag))
    lam_ii = (3*state.attractor.k*state.r[0]*state.r[0])/(r_mag5)
    lam_jj = (3*state.attractor.k*state.r[1]*state.r[1])/(r_mag5)
    lam_kk = (3*state.attractor.k*state.r[2]*state.r[2])/(r_mag5)
    lam_ij = (3*state.attractor.k*state.r[0]*state.r[1])/(r_mag5)
    lam_ik = (3*state.attractor.k*state.r[0]*state.r[2])/(r_mag5)
    lam_jk = (3*state.attractor.k*state.r[1]*state.r[2])/(r_mag5)
    lam = np.array([[lam_ii+update_pos, lam_ij, lam_ik],[lam_ij, lam_jj+update_pos, lam_jk], [lam_ik, lam_jk, lam_kk+update_pos]])
    return np.array([[0, np.identity(6)],[lam, 0]])

def stm(state):
    r_mag = np.sqrt(state.r[0]**2 + state.r[1]**2 + state.r[2]**2) 
    r_mag = r_mag.to(u.m)
    r_mag5 = r_mag*r_mag*r_mag*r_mag*r_mag
#    v_mag = np.linalg.norm(state.v)
    update_pos = -(state.attractor.k/(r_mag*r_mag*r_mag))
    lam_ii = ((3*state.attractor.k*state.r[0]*state.r[0])/(r_mag5))
    lam_jj = ((3*state.attractor.k*state.r[1]*state.r[1])/(r_mag5))
    lam_kk = ((3*state.attractor.k*state.r[2]*state.r[2])/(r_mag5))
    lam_ij = ((3*state.attractor.k*state.r[0]*state.r[1])/(r_mag5))
    lam_ik = ((3*state.attractor.k*state.r[0]*state.r[2])/(r_mag5))
    lam_jk = ((3*state.attractor.k*state.r[1]*state.r[2])/(r_mag5))
    lam = np.array([[lam_ii+update_pos, lam_ij, lam_ik],[lam_ij, lam_jj+update_pos, lam_jk], [lam_ik, lam_jk, lam_kk+update_pos]])
    return lam

Y_comp, X_comp = hh_obs()
X_STM = stm(X_comp)