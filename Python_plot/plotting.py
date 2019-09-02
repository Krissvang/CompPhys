# -*- coding: utf-8 -*-
"""
Created on Sat Aug 31 15:37:43 2019

@author: Kristian
"""

import matplotlib.pyplot as plt
import numpy as np


    

def f(t):
    return 1-(1-np.exp(-10))*t-np.exp(-10.0*t)
n=3

#Lists for names
filenames=[]
read_files=[]

#Create file paths
path="..\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output"
for i in range(n):
    filenames.append(path+"%d.txt"%(i+1))
    read_files.append(open(filenames[i],"r"))

    flines=read_files[i].readlines()
    xval, solution, exact, error = np.loadtxt(path+"%d.txt"%(i+1), unpack=True) 
    
    t1=np.arange(xval[0],xval[-1],0.01)
    plt.plot(xval,solution,'r.',markersize=3)
    
    plt.title("n=%d"%(pow(10,i)))
    plt.plot(t1,f(t1),'g')
    plt.xlabel('Position x')
    plt.ylabel('u(x)')
    plt.savefig('n%dplot.png'%(pow(10,i)))
    plt.show()