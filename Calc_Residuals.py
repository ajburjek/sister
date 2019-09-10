#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 10 14:06:26 2019

@author: willjpatton
"""

import numpy as np


'Ground Station Location
Xg = 
Yg = 
Zg = 

interval = 1 #sec

def Range ():
    #x y and z values from range and ground station
    range = np.sqrt(((X-Xg)**2)+((Y-Yg)**2)+((Z-Zg)**2))
    return range

#elevation

def Time_Recieved (range,interval):
        # time recieved from range calculation
        c = 299792458 #m/s
        time_received = interval + range/c;
        return time_recieved

def Alter_Time (time_received):
    #change time recieved by some value for testing
    time_altered = time_recieved*1.1
    return time_altered

def Residuals (time_received,time_altered):
    #return time residuals by taking differences of time recieved and altered
    residual = np.abs(time_received-time_altered)
    return residual
'''
need rt and zt to find elevation


if -90 < El < 90 deg
elevation = np.asin(z/r)
else
end
'''

range = Range()

time_recieved = Time_Received(range,interval)

time_altered = Alter_Time(time_received)

residual = Residuals(time_received,time_altered)



