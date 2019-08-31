# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 15:37:43 2019

@author: Kristian
"""

import matplotlib.pyplot as plt

f1=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output1.txt","r")
f2=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output2.txt","r")
f3=open("E:\\Dokumenter\\compphys\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output3.txt","r")

flines=f2.readlines()
xval = []
solution = []
exact = []
error = []


for x in flines:
    n=0
    for i in x.split(' '):
        if i=='':
            pass
        else:
            j=float(i)
            if n==0:
                xval.append(j)
            elif n==1:
                solution.append(j)
            elif n==2:
                exact.append(j)
            else:
                error.append(j)
            n+=1

plt.plot(xval,solution)
