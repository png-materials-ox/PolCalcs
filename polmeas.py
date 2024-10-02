# -*- coding: utf-8 -*-
"""
Created on Sat Jul 13 19:09:48 2024

@author: mans4209
"""

import numpy as np
import scipy.constants as sc

from matplotlib import pyplot as plt
from jonescalculus import JonesCalculus

    
class PolMeasurement:
    
    def __init__(self, Evec=[1,1,0], wlen=600e-09, t=[1], z=[1],
                phix = 0, phiy=120):
        
        self.jones = JonesCalculus(Evec=Evec, wlen=622e-09, t=[1], z=[1],
                phix = 20, phiy=20)
        
        self.E0 = self.jones.E0
        
    def Eqwp(self, theta=[]):
        if len(theta) == 0:
            return self.jones.qwp() * self.E0
        else:
            return self.jones.qwp_rotating(theta)[0] * self.E0
            
        
    def Ehwp(self, E, theta=[0]):
        Jhwp = self.jones.hwp_rotating(theta)
        Ehwp = []

        for i in range(len(theta)):
            Ehwp.append(Jhwp[i] * E)
        return Ehwp
    
    def Elp(self, E):
        Elp = []
        for e in E:
            Elp.append(self.jones.lp() * e)
        return Elp
    
    def E_birefringent_material(self, E, theta=[0], eta=1, phi=0):
        Jbr = self.jones.birefringent_material(eta=eta, theta=theta, phi=phi)
        Ebr = []

        for i in range(len(theta)):
            Ebr.append(Jbr[i] * self.E0)
        return Ebr


    def I(self, E):
        I = []
        for Ef in E:
            In = Ef[0]**2
            I.append(In.tolist()[0][0])
        return I

    def plot_polar(self, I):
        # Generating sample data
        theta = np.linspace(0, 2*np.pi, len(I))
        r = I

        # Creating the polar scatter plot
        fig = plt.figure(figsize=(8, 8))
        ax = plt.subplot(111, polar=True)
        ax.plot(theta, r)

        plt.show()
        
    
