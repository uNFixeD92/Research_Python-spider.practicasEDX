#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 14:58:21 2021

@author: fmarcosdev
"""

# random walks :o
import random;
import matplotlib.pyplot as plt;
import numpy as np;


recorrido = np.random.normal(0,1,(2,5)) #crea matriz de 2x5 con valores del 0 al 1 en cada puesto
X_origin = np.array([[0],[0]]) # posicion base para mi array

X_walks = np.cumsum(recorrido , axis = 1)
X = np.concatenate((X_origin,X_walks),axis=1)	# concatenar numpy arrays 
plt.plot(X[0],X[1],"ro-") # grafica x0 en horizontal y x1 en vertical
