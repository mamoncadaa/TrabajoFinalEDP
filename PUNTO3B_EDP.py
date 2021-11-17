# -*- coding: utf-8 -*-
"""Punto3b.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1ki58LUsup_lLNlr7Z00A-oXHb6R5bHnl
"""

import numpy as np
from sympy import *
import math as m
from scipy.sparse import spdiags
import matplotlib.pyplot as plt
#1.b
c = 2
u_x_0 = lambda x: exp(-x)
ut_x_0 = lambda x: -2*exp(-x)
u_0_t = lambda t:exp(-2*t)
u_1_t = lambda t: exp(-1-2*t)

h = 2**(-7)
k = h/(c*2)
sigma = c*k/h

inth = int(1/h)
intk = int(1/k)

n = np.ones((inth,inth))
diag = np.identity(inth)*(2-2*sigma**2)
diagSup = np.triu(n, k=1) - np.triu(n, k=2) 
diagInf = np.tril(n,k=-1) - np.tril(n,k=-2)
A = diag + (sigma**2)*diagSup + (sigma**2)*diagInf

#j=1 y x in (0,1)
xi = np.linspace(0,1,inth)
ti = np.linspace(0,1,intk)

fxi = np.array([u_x_0(fx) for fx in xi])
gxi = np.array([ut_x_0(gx) for gx in xi])

tj = np.zeros(inth)
tj[0] = u_0_t(0)
tj[-1] = u_1_t(0)

wij = np.zeros((intk,inth))
wi0 = fxi
wi1 = wi0
wi1 = 1/2*np.matmul(A,fxi) + k*gxi + 1/2*sigma**2*tj
wij[0] = wi0
wij[1] = wi1
#j!=1
for i in range(2,intk):
  tj[0] = u_0_t(ti[i-1])
  tj[-1] = u_1_t(ti[i-1])
  wij[i] = np.matmul(A,wij[i-1]) - wij[i-2] + sigma**2*tj

wij_tras = wij.T

X,Y = np.meshgrid(np.linspace(0,1,intk), np.linspace(0,1,inth))
Z = wij_tras

graf = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_wireframe(X,Y,Z)
ax.set_title('Aproximación')

tj

a = np.linspace (0 , 1, inth)
b = np.linspace (0 , 1, intk)
X, Y = np.meshgrid(a, b)
Zreal = f(X, Y)
Z = wij

def f(a, b):
    return np.exp(-a-2*b)

graf = plt.figure()
ax = plt.axes(projection ='3d')
ax.plot_wireframe(X, Y, Zreal)
ax.set_title('Función real');

a.size