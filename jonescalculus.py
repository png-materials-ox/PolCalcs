# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 19:09:48 2024

@author: mans4209
"""

import numpy as np
import scipy.constants as sc

from matplotlib import pyplot as plt

class JonesCalculus:
    
    def __init__(self, Evec=[1,1,1], wlen=600e-09, t=[1], z=[1],
                phix = 0, phiy=120):
        self.Ex = Evec[0]
        self.Ey = Evec[1]
        self.Ez = Evec[2]
    
        self.wlen = wlen
        self.w = sc.c/wlen
        self.k = (2*np.pi)/wlen
        self.t = t
        self.z = z
        
        self.phix = phix
        self.phiy = phiy
        
        self.E0x = self.Ex * np.cos(np.deg2rad(self.phix))
        self.E0y = self.Ey * np.cos(np.deg2rad(self.phiy))
        self.E0 = np.matrix([[self.E0x], [self.E0y]])
    
    def qwp(self):
        return np.exp(-1j*np.pi/4) * np.matrix('1,0;0,1j')
    
    def qwp_rotating(self, theta):
        Jqwp = []
        for a in theta:
            ang = np.deg2rad(a)
            Jqwp.append(np.exp(-1j*np.pi/4) * np.matrix([[np.cos(ang)**2+1j*np.sin(ang)**2, (1-1j)*np.sin(ang)*np.cos(ang)],
                                [(1-1j)*np.sin(ang)*np.cos(ang), np.sin(ang)**2+1j*np.cos(ang)**2]]))
        return Jqwp
    
    def hwp_rotating_general(self, theta):
        Jhwp = []
        for a in theta:
            ang = np.deg2rad(a)
            Jhwp.append(np.exp(-1j*np.pi/2) * np.matrix([np.cos(ang)**2-np.sin(ang)**2, 2*np.cos(ang)*np.sin(ang)],
                                            [2*np.cos(ang)*np.sin(ang), np.sin(ang)**2-np.cos(ang)**2]))
        return Jhwp
            
    
    def hwp_rotating(self, theta):
        Jhwp = []
        for ang in theta:
            Jhwp.append(np.matrix([[np.cos(2*np.deg2rad(ang)) , np.sin(2*np.deg2rad(ang))], 
                      [np.sin(2*np.deg2rad(ang)), -np.cos(2*np.deg2rad(ang))]]))
        return Jhwp
    
    def lp(self):
        return np.matrix('1,0;0,0')
    
    def birefringent_material(self, eta=1, theta=[], phi=0):
        J = []

        for ang in theta:
            ang = np.deg2rad(ang)
            term1 = np.cos(ang)**2 + np.exp(1j*eta) * np.sin(ang)**2
            term2 = (1-np.exp(1j*eta)) * np.exp(-1j*phi)*np.cos(ang)*np.sin(ang)
            term3 = (1-np.exp(1j*eta)) * np.exp(1j*phi)*np.cos(ang)*np.sin(ang)
            term4 = np.sin(ang)**2 + np.exp(1j*eta) * np.cos(ang)**2

            J.append(np.exp(-1j*eta/2) * np.matrix([[term1, term2],[term3, term4]]))
        
        return J
