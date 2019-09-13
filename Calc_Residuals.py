#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:06:26 2019

@author: willjpatton
"""
import numpy as np
from astropy import constants as const

def Time_Received (rho,interval):
        # time received from range calculation
        time_received = interval + rho/const.c;
        return time_received

def Alter_Time (time_received):
    #change time recieved by some value for testing
    time_altered = time_received*np.random.random_sample()
    return time_altered

def Residuals (time_received,time_altered):
    #return time residuals by taking differences of time recieved and altered
    residual = np.abs(time_received-time_altered)
    return residual
