#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
// use namespace for output and input
using namespace std;

// object for output files
ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

//Triagonal matrix algorithm

double * t1(double* a, double* b, double* d, double* solution, double* g, int n){
    // Forward substitution
    for (int i = 0; i < n; ++i) {
        d[i+1]-=a[i]*b[i]/d[i];
        g[i+1]-=a[i]*g[i]/d[i];
    }
    // Backward substitution
    solution[n-1]=g[n-1]/d[n-1];
    for (int i = 2; i < n+1; ++i) {
        solution[n-i]=(g[n-i]-b[n-i]*solution[n+1-i])/d[n-i];}
    return solution;
}

//Specialized matrix algorithm
//Here the diagonal elements are identical, and the of diagonal elements are identical.
double * t2(double* a, double *d, double* solution, double* g, int n){
    double e=a[0];
    double ee=e*e;
    // Forward substitution
    for (int i = 0; i < n; ++i) {
        d[i+1]-=ee/d[i];
        g[i+1]-=e*g[i]/d[i];
    }
    //Backward substitution
    solution[n-1]=g[n-1]/d[n-1];
    for (int i = 2; i < n+1; ++i) {
        solution[n-i]=(g[n-i]-e*solution[n+1-i])/d[n-i];
    }
    return solution;

}



// Begin main program
int main(int argc, char *argv[]){
  int exponent;
  clock_t start1, finish1, start2, finish2;
    string filename;
    // We read also the basic name for the output file and the highest power of 10^n we want
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max power 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
    }
    // Loop over powers of 10
    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      // Declare new file name
      string fileout = filename;
      // Convert the power 10^i to a string
      string argument = to_string(i);
      // Final filename as filename-i-
      fileout.append(argument);
      fileout.append(".txt");
      double h = 1.0/(n+1);
      double hh = h*h;
      // Set up arrays for the simple case
      double *d = new double [n]; double *g = new double [n]; double *a = new double [n-1];
      double *b = new double [n-1];
      double *solution = new double [n+2]; double *x = new double[n+2];
      // Quick setup of updated diagonal elements and value of
      solution[0] = solution[n+1] = 0.0;

      for (int i = 0; i <= n+1; i++){
          x[i]= i*h;
      }

      for (int i = 0; i < n; i++) {
          d[i] = 2;             // We get that d[0]=d_1
          g[i] = hh*f(x[i+1]);  // We get that g[0]=g_1
      };

      for (int i=0; i < n-1; i++) {
          a[i]=-1;
          b[i]=-1;
      }

      if(atoi(argv[3])==0){
        start1=clock();
        solution=t1(a,b,d,solution,g,n);
        finish1=clock();
        double time1=(double) (finish1-start1)/(CLOCKS_PER_SEC);
        cout<<"When n = " << n << ", time used for general algorithm is = "<< time1 << endl;
      }
      else if (atoi(argv[3])==1) {
        start2=clock();
        solution=t2(a,d,solution,g,n);
        finish2=clock();
        double time2=(double) (finish2-start2)/(CLOCKS_PER_SEC);
        cout <<"When n = "<< n << ", time used for specialized algorithm is "<< time2 << endl;
      }
      else{
          cout << "Please choose a algorithm to run. Choose with the third command line argument."<<endl;
          cout << "0 = General tridiagonal matrix set of linear equations"<<endl;
          exit(0);
      }
      //exact
      double diff = 0.0;
      for (int i = 0; i <= n; i++) {
          diff += fabs(solution[i]-exact(x[i]));
      }
      //cout << "By averaging we get that the average error is "<< diff/n << endl;

      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 1; i < n;i++) {
    double xval = x[i];
     double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
         ofile << setw(15) << setprecision(8) << xval;
         ofile << setw(15) << setprecision(8) << solution[i];
         ofile << setw(15) << setprecision(8) << exact(xval);
         ofile << setw(15) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();
      delete [] x; delete [] d; delete [] g; delete [] a; delete [] solution;
    }
    return 0;
}
