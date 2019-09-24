#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Repeating Ground Track Orbits in High-Fidelity Geopotential"""
import numpy as np
import os
import string
from astropy import constants as const
from astropy.time import Time
"""Python 3.7
   Simpson Aerospace, Copyright 2019
   Christopher R. Simpson
   simpsonchristo@gmail.com """
#------------------------------------------------------------------------------
def ismember(A,B):
    """ISMEMBER Returns an array containing logical 1 (true) where the data in A is found in B."""
    return [np.sum(a==B) for a in A]

def chksum(linewords):
    """Checksum (Modulo 10): Letters, blanks, periods, plus signs = 0; minus signs = 1"""
    L = linewords
    c = 0
    result = 0
    
    for k in range(68):
        if(L[k]>'0' and L[k]<='9'):
            c += int(L[k]) - 48
        elif(L[k]=='-'):
            c+=1
    if(np.remainder(c,10)==int(L[68])-48):
        result = 1
    return result

def readtle(file, catalog=np.zeros(1)):
    """Read NORAD two-line element (TLE) file. Adapted from Mr. Brett Pantalone's MATLAB function."""
    #   Adapted from Mr. Brett Pantalone's MATLAB readtle function.
    #   INPUTS:
    #   file    - Path to any standard two-line element file.
    #   catalog - Optional array of NORAD catalog numbers for the satellites of
    #            interest. The default action is to display data from every
    #            satellite in the file.
  
    #   Brett Pantalone
    #   mailto:bapantal@ncsu.edu
    #   http://research.ece.ncsu.edu/osl/
    if(catalog==np.zeros(1)):
        catalog = np.array([])
    
    try:
        fd = open(file)
    except(IOError):
        file = file[0:-3]+".tle"
        fd = open(file)
    
    assert(os.path.isfile('./'+file), "File doesn''t exist in this directory.")
    
    kiter = 0
    A0 = fd.readline().rstrip()
    A1 = fd.readline().rstrip()
    A2 = fd.readline().rstrip()
    oe = np.array([0, 0, 0, 0, 0, 0])
    epoch = np.array([0, 0])
    
    try:
        while(isinstance(A2,str)==1):
            kiter+=1
            satnum = np.array([float(A1[2:6])])
            if(catalog.size==0 or ismember(satnum,catalog)==1):
                if(kiter==1):
                    print('-'*50)
                    print('Satellite: %s' % A0)
                    assert(chksum(A1), 'Checksum failure on line 1')
                    assert(chksum(A2), 'Checksum failure on line 2')
                    print("Catalog Number: %f" % satnum)
                    epochyear = np.array([float('20'+A1[18:20])])
                    epochday  = np.array([float(A1[20:31])])
                    epoch = np.array([epochyear,epochday])
                    print("Epoch time: %s" % A1[18:31]) #YYDDD.DDDDDDDD
                    inc = np.array([float(A2[8:15])])
                    print("Inclination: %f deg" % inc)
                    raan = np.array([float(A2[17:24])])
                    print("Right Ascension of the Ascending Node: %f deg" % raan)
                    ecc = np.array([float('.' + A2[26:32])])
                    print("Eccentricity: %f" % ecc)
                    aop = np.array([float(A2[34:41])])
                    print("Argument of perigee: %f deg" % aop)
                    M = np.array([float(A2[43:50])])
                    print("Mean Anomaly: %f deg" % M)
                    n = np.array([float(A2[52:62])])
                    print("Mean motion: %f rev/day" % n)
                    T = 86400/n;
                    print("Period of rev: %.0f s/rev" % T)
                    a = (((T/(2*np.pi))**2)*3.986004e+14)**(1/3);
                    print("Semimajor axis: %.0f meters" % a)
                    b = a*(1 - ecc**2)**(0.5)
                    print("Semiminor axis: %.0f meters" % b)
                    oe = np.array([a, ecc, inc, raan, aop, M])
                elif(kiter>1):
                    print('-'*50)
                    print('Satellite: %s' % A0)
                    assert(chksum(A1), 'Checksum failure on line 1')
                    assert(chksum(A2), 'Checksum failure on line 2')
                    print("Catalog Number: %f" % satnum)
                    epochyear = np.array([float('20'+A1[18:20])])
                    epochday  = np.array([float(A1[20:31])])
                    epoch_new = np.array([epochyear,epochday])
                    print("Epoch time: %s" % A1[18:31]) #YYDDD.DDDDDDDD
                    inc = np.array([float(A2[8:15])])
                    print("Inclination: %f deg" % inc)
                    raan = np.array([float(A2[17:24])])
                    print("Right Ascension of the Ascending Node: %f deg" % raan)
                    ecc = np.array([float('.' + A2[26:32])])
                    print("Eccentricity: %f" % ecc)
                    aop = np.array([float(A2[34:41])])
                    print("Argument of perigee: %f deg" % aop)
                    M = np.array([float(A2[43:50])])
                    print("Mean Anomaly: %f deg" % M)
                    n = np.array([float(A2[52:62])])
                    print("Mean motion: %f rev/day" % n)
                    T = 86400/n;
                    print("Period of rev: %.0f s/rev" % T)
                    a = (((T/(2*np.pi))**2)*3.986004e+14)**(1/3);
                    print("Semimajor axis: %.0f meters" % a)
                    b = a*(1 - ecc**2)**(0.5)
                    print("Semiminor axis: %.0f meters" % b)
                    oe_new = np.array([a, ecc, inc, raan, aop, M])
                    oe = np.concatenate((oe,oe_new), axis=1)
                    epoch = np.concatenate((epoch,epoch_new),axis=1)
            A0 = fd.readline().rstrip()
            A1 = fd.readline().rstrip()
            A2 = fd.readline().rstrip()
    except:
        fd.close()
    return oe, epoch