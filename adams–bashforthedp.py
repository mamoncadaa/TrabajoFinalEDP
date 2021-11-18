# -*- coding: utf-8 -*-
"""Adams–BashforthEDP.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1VjcIHXUr93u9HuxjxuaJHGi5V1OC0sQ5
"""

from math import *
from tabulate import tabulate
import numpy as np
import matplotlib.pyplot as plt
import scipy as sip
import sympy as sp
from sympy import *

'''
Este código hace uso del método RKF que se encuentra en este mismo repositorio 
'''
##Funciones PVI
def func1(t,y1,y2):
  return (9*y1)+(24*y2)+(5*cos(t))-((1/3)*sin(t))
def func2(t,y1,y2):
  return -(24*y1)-(51*y2)-(9*cos(t))+((1/3)*sin(t))
  
##Solución PVI
def funcSol1(t):
  return (2*e**(-3*t))-(e**(-39*t))+(1/3*cos(t))
def funcSol2(t):
  return (-e**(-3*t))+(2*e**(-39*t))-(1/3*cos(t))

def AdBash3(func1,func2,funcSol1,funcSol2,a,b,n,y1,y2):
    h = b-a/n
    x, rk1, rk2 = RKF(func1,func2,funcSol1,funcSol2,a,b,n,y1,y2)
    re1 = np.zeros(n+1)
    re2 = np.zeros(n+1)
    re1[0:3] = rk1[0:3]
    re2[0:3] = rk2[0:3]
    Ey1 = np.zeros(n+1)
    Ey2 = np.zeros(n+1) 
    for i in range(3,n+1):
      re1[i] = re1[i-1] + h*((5/12)*func1(x[i-3],re1[i-3], re2[i-3])-(4/3)*func1(x[i-2],re1[i-2], re2[i-2])+(23/12)*func1(x[i-1],re1[i-1], re2[i-1]))
      re2[i] = re1[i-1] + h*((5/12)*func2(x[i-3],re1[i-3], re2[i-3])-(4/3)*func2(x[i-2],re1[i-2], re2[i-2])+(23/12)*func2(x[i-1],re1[i-1], re2[i-1])) 
      Ey1[i] = abs(re1[i] -funcSol1(x[i]))
      Ey2[i] = abs(re2[i] -funcSol2(x[i]))
    print(tabulate({'x': x, 'y1': re1,'Error y1':Ey1,'y2': re2,'Error y2':Ey2 }, headers="keys", tablefmt='fancy_grid'))
    return x, re1, re2
AdBash3(func1,func2,funcSol1,funcSol2,0,1,4,4/3,2/3)

x = np.linspace(0,1,1000)
real1 = np.zeros(len(x))
real2 = np.zeros(len(x))
for i in range(len(x)):
  real1[i] = funcSol1(x[i])
  real2[i] = funcSol2(x[i])

x1,y11,y21=AdBash3(func1,func2,funcSol1,funcSol2,0,1,2,4/3,2/3)
x2,y12,y22=AdBash3(func1,func2,funcSol1,funcSol2,0,1,4,4/3,2/3)
x3,y13,y23=AdBash3(func1,func2,funcSol1,funcSol2,0,1,8,4/3,2/3)
x4,y14,y24=AdBash3(func1,func2,funcSol1,funcSol2,0,1,16,4/3,2/3)
x5,y15,y25=AdBash3(func1,func2,funcSol1,funcSol2,0,1,32,4/3,2/3)
x6,y16,y26=AdBash3(func1,func2,funcSol1,funcSol2,0,1,64,4/3,2/3)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
plt.plot(x,real1,'b', label='Exacta')

plt.plot(x1,y11,'g', label='h=1/2')

plt.plot(x2,y12,'r', label='h=1/4')
plt.plot(x3,y13,'c', label='h=1/8')
plt.plot(x4,y14,'m', label='h=1/16')
plt.plot(x5,y15,'y', label='h=1/32')
plt.plot(x6,y16,'k', label='h=1/64')
plt.ylim(0,2)
plt.legend(loc='upper right')
plt.show()