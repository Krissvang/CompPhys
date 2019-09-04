import matplotlib.pyplot as plt
import numpy as np

    
#Creating exact solution
def f(t):
    return 1-(1-np.exp(-10))*t-np.exp(-10.0*t)
n=3
alg=2

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
for i in range(n):
    #Create file names
    filenames.append(path+"-alg-%d-n=%d.txt"%(alg,pow(10,i+1))) 
    #Open the files
    read_files.append(open(filenames[i],"r"))                   
    #Create vectors with values                     
    xval, solution, exact, error = np.loadtxt(
    path+"-alg-%d-n=%d.txt"%(alg,pow(10,i+1)), unpack=True)     
    #Interval for exact solution
    t1=np.arange(xval[0],xval[-1],0.01)
    #Plot the numerical solution in red dots
    plt.plot(xval,solution,'r.',markersize=3)
    #Plot exact solution in green
    plt.plot(t1,f(t1),'g')
    #Plot configs
    plt.title(Title + " n=%d"%(pow(10,i+1)))
    plt.xlabel('Position x')
    plt.ylabel('u(x)')
    #Save as png in Latex folder
    plt.savefig('../../Latex/img/alg-%d-n%dplot.png'%(alg,pow(10,i+1)))
    plt.show()
    
    #Close files
for i in range(n):
    read_files[i].close()