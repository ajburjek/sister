#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 12:57:08 2019

@author: William Patton

Attempt at Kalman Filter
"""

import pylab as pl
import numpy as np
from pykalman import AdditiveUnscentedKalmanFilter

'Initialize Parameters'

def transition_function(X):
    ...

def observation_function(X):
    ...

transition_covariance = np.eye(2)
    
observation_covariance = np.eye(2) + something

initial_state_mean = [0, 0]

initial_state_covariance = [[1, 0.1], [ 0.1, 1]]

akf = AdditiveUnscentedKalmanFilter(transition_function,observation_function, \
    transition_covariance,observation_covariance,initial_state_mean, \
    initial_state_covariance)


akf_state_estimates = akf.filter(timesteps,states)

pl.figure()
lines_true = pl.plot(states, color='b')
lines_akf = pl.plot(akf_state_estimates, color='g', ls='-.')
pl.show()