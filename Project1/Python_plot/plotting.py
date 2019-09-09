import matplotlib.pyplot as plt
import numpy as np
import os
    
#Creating exact solution
def f(t):
    return 1-(1-np.exp(-10))*t-np.exp(-10.0*t)

#Choose the number of gridpoints by gridpoints = 10^n
n=6
#Choose which algorithm to run, 0=general, 1 =specialized and 2)=LU.
alg=1
Title=""
#Lists for names
filenames=[]
read_files=[]

#Creates title depending on algorithm choosen
if alg==0:
    Title="Generelized algorith"
elif alg==1:
    Title="Spesialized algorithm"

elif alg==2:
    Title="LU decomposition"


#Create file paths
path="..\\..\\project1\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\output"
path2="..\\..\\project1\\build-Project1-Desktop_Qt_5_13_0_MinGW_64_bit-Debug\\time"
maxerrors = list()
#Create file names
plt.close('all')
x_val=[]

Errors=[-1.1797, -3.08804, -5.08005, -7.07927, -9.0791, -10.791, -9.54234]

def plot_solution():
    
    for i in range(n):
        filenames.append(path+"-alg-%d-n=%d.txt"%(alg,pow(10,i+1))) 
        #Open the files
        read_files.append(open(filenames[i],"r"))                   
        #Create vectors with values                     
        xval, solution, exact, error = np.loadtxt(
        path+"-alg-%d-n=%d.txt"%(alg,pow(10,i+1)), unpack=True)
        np.delete(error,0)
        error=error[~np.isnan(error)]
        maxerrors.append(np.max(error))
        #Interval for exact solution
        t1=np.arange(xval[0],xval[-1],0.001)
        #Plot the numerical solution in red dots
        plt.figure(i+1)
        plt.plot(xval,solution,'r.',markersize=3, label="Numerical solution")
        #Plot exact solution in green
        plt.plot(t1,f(t1),'g', label= "Exact solution")
        #Plot configs
        plt.grid(1)
        plt.title(Title + " n=%d"%(pow(10,i+1)))
        plt.xlabel('Position x')
        plt.ylabel('u(x)')
        plt.legend()
        #Save as png in Latex folder
        plt.savefig('../../Latex/img/alg-%d-n%dplot.png'%(alg,pow(10,i+1)),dpi=199)
        plt.show()
        #Close files
    for i in range(n):
        read_files[i].close()

def relative_error():
    #Create values for the x axis.
    for i in range(7):
        x_val.append(np.log10((pow(10,i+1)+2)))
    #Plots and configs of plot
    plt.plot(x_val,Errors,'ro-')
    plt.title("Relative error")
    plt.grid(1)
    plt.xlabel('$\log{(n)}$')
    plt.ylabel('$\epsilon$')
    plt.savefig('../../Latex/img/relative_error.png',dpi=199)


def timing():
    #Prints out the mean times and the standard deviations.
    for i in range(n):
        time = np.loadtxt(
        path2+"-alg-%d-n=%d.txt"%(alg,pow(10,i+1)), unpack=True)
        time=time*pow(10,-9)
        mean_time=np.mean(time)
        mean_time=round(mean_time,10)
        std=round(np.std(time),10)
        print("Mean is ", mean_time, " and std is ", std )

            
timing()
        
        
        
        
