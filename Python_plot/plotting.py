# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 15:37:43 2019

@author: Kristian
"""

import matplotlib.pyplot as plt
import numpy as np


def f(t):
    return 1-(1-np.exp(-10))*t-np.exp(-10.0*t)

f1=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output1.txt","r")
f2=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output2.txt","r")
f3=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output3.txt","r")

file=f2

flines=file.readlines()
xval = np.array([])
solution = np.array([])
exact = np.array([])
error = np.array([])

for x in flines:
    n=0
    for i in x.split(' '):
        if i=='':
            pass
        else:
            j=float(i)
            if n==0:
                xval=np.append(xval,j)
            elif n==1:
                solution=np.append(solution,j)
            elif n==2:
                exact=np.append(exact,j)
            else:
                error=np.append(error,j)
            n+=1

t1=np.arange(xval[0],xval[-1],0.01)
print(exact[0])


plt.plot(xval,solution,'r.',markersize=3)
plt.title("n=100")
plt.plot(t1,f(t1),'g')
plt.xlabel('Position x')
plt.ylabel('u(x)')
plt.savefig('n100plot.png')