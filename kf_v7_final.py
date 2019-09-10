# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 09:42:02 2019

@author: klbryan

clean version of kf_v6.py
working linear Kalman Filter
p779 of Vallado
"""

import numpy as np

class KalmanFilter(object):
    
    def __init__(self, X, P, sigma, H, Q, R):
        self.X = X
        self.P = P
        self.sigma = sigma
        self.H = H
        self.Q = Q
        self.R = R
        self.I = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
	
    def predict(self, z, delta_t):
        self.z = z
        self.delta_t = delta_t
        self.phi =  np.array([[1., self.delta_t, ((self.delta_t) ** 2)/2], [0., 1., self.delta_t], [0., 0., 1.]])
        
        self.X_bar = np.dot(self.phi, self.X)
        self.P_bar = np.dot(np.dot(self.phi, self.P), self.phi.T) + self.Q
        
    def update(self):
        H_t = np.array([[0.],[0.],[0.]])
        for i in range(len(self.H)):   
            H_t[i] = self.H[i]
        
        K_int = 1/((np.dot(np.dot(self.H, self.P_bar), H_t)) + self.R)
        self.K = np.dot(self.P_bar, H_t)*K_int
        self.b = self.z - np.dot(self.H, self.X_bar)
        self.X_hat = self.X_bar + self.K*self.b
        self.H = self.H.reshape(3,1)
        P_hat_int = self.I - np.outer(self.K, self.H)
        self.P_hat = np.dot(P_hat_int, self.P_bar)
        self.std = np.sqrt(self.P_hat)
        
    def printresults(self):
        print "State estimate = ", self.X_hat
        print "Standard deviation [pos, v, a] = [", self.std[0][0], " ", self.std[1][1], " ", self.std[2][2], "]"
        print "Residual (m) = ", self.b

def main():
    t = 0
    
    X = np.array([[600.], [-45.], [1.5]])
    std = np.array([[15., 0., 0.], [0., 5., 0.], [0., 0., 1.]])
    P = std ** 2
    sigma = 5.
    H = np.array([1., 0., 0.])
    noise = np.array([[1., 0., 0.], [0., 0.1, 0.], [0., 0., 0.1]])
    Q = noise ** 2
    R = 25.
    
    z = input('Enter the range observed (m), or press Enter to exit: ')
    delta_t = input('Enter the time interval between the most recent range measurements (s): ')
    delta_t = float(delta_t)
    
    while z:
        t = t + delta_t
        z = float(z)
        KF = KalmanFilter(X, P, sigma, H, Q, R)
        KF.predict(z, delta_t)
        KF.update()
        print "\nt = ", t
        KF.printresults()
        
        #next measurement 
        X = KF.X_hat
        P = KF.P_hat
        z = input('Enter the range observed (m), or press Enter to exit: ')
        delta_t = input('Enter the time interval between the most recent range measurements (s): ')
        delta_t = float(delta_t)

if __name__ == "__main__":
    main()