import numpy as np



filenames=[]
read_files=[]


n=5
alg=1

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
    print(max(error))
    