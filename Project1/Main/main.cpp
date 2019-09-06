#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include "time.h"
#include <armadillo>
// use namespace for output and input
using namespace std;
using namespace arma;
// object for output files
ofstream ofile;
// Functions used
inline double f(double x){return 100.0*exp(-10.0*x);
}
inline double exact(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}

//General triagonal matrix algorithm
double * t1(double* a, double* b, double* d, double* solution, double* g, int n){
    // Forward substitution
    for (int i = 0; i < n; ++i) {
        d[i+1]-=a[i]*b[i]/d[i];
        g[i+1]-=a[i]*g[i]/d[i];
    }
    // Backward substitution
    solution[n]=g[n-1]/d[n-1];
    //Choose i=2 to i=n because that gives the least amount of operations.
    for (int i = 2; i < n+1; ++i) {
        solution[n-i+1]=(g[n-i]-b[n-i]*solution[n+2-i])/d[n-i];}
    return solution;
}

//Specialized matrix algorithm
//Here the diagonal elements are identical, and the of diagonal elements are identical.
double * t2(double* a, double *d, double *d2, double* solution, double* g, double* g2, int n){
    //Saves n FLOPS to calculate e*e before the substitutions.
    // Forward substitution
    int e = a[0];
    int ee=e*e;

    //We have a special algorithm when the off diagonal
    // elements are -1 and the diagonal elements are 2.
        if(e==-1 && d[0]==2){
            //Forward substitution
            for (int i = 0; i < n; ++i) {

                g[i+1]+=g[i]/d2[i];
            }
            //Backward substitution
            solution[n]=g[n-1]/d2[n-1];
            for (int i = 2; i < n+1; ++i) {
                solution[n-i+1]=(g[n-i]+solution[n+2-i])/d2[n-i];

            }

        }
        else {
            //Forward substitution
            for (int i = 0; i < n; ++i) {
                d[i+1]-=ee/d[i];
                g[i+1]-=e*g[i]/d[i];
            }
            //Backward substitution
            solution[n]=g[n-1]/d[n-1];
            for (int i = 2; i < n+1; ++i) {
                solution[n+1-i]=(g[n-i]+solution[n+2-i])/d[n-i];
            }

        }

    return solution;

}

//Armadillo LU decomposition
double * t3(double* a, double * b, double *d, double *solution, double* g, unsigned int n){
    Mat<double> A(n,n,fill::zeros);
    for (unsigned int i = 0; i < n; ++i) {
        A(i,i)=d[i];
    }
    for (unsigned int i = 0; i < n-1; ++i) {
        A(i,i+1)=a[i];
        A(i+1,i)=b[i];
    }
    Mat<double> L, U;
    lu(L, U, A);
    double * y = new double [n];

    //Solving Ly=g
    y[0]=g[0]/L(0,0);
    for (unsigned int i = 1; i < n; ++i) {
        y[i]= g[i];
        for (unsigned int j = 0; j < i; ++j) {
            y[i]-=L(i,j)*y[j];
        }
        y[i]/=L(i,i);
    }
    //Solving Ux=y
    solution[n]=y[n-1]/U(n-1,n-1);
    for (unsigned int i = 2; i < n+1; ++i) {
        solution[n-i+1]=y[n-i];
        for (unsigned int j = 1; j < i; ++j) {
            solution[n-i+1]-=U(n-i,n-j)*solution[n-j+1];
        }
        solution[n-i+1]/=U(n-i,n-i);
    }
    delete [] y;
    return solution;
}

// Begin main program
int main(int argc, char *argv[]){
  int exponent;
  clock_t start1, start2, start3, finish1, finish2, finish3;
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
      fileout.append("-");
      // Convert the power 10^i to a string
      string argument = to_string(int(pow(10,i)));
      // Final filename as filename-i-

      double h = 1.0/(n+1);
      double hh = h*h;

      // Set up arrays for the simple case
      double *d = new double [n]; double *g = new double [n]; double *a = new double [n-1];
      double *b = new double [n-1]; double *d2 = new double [n]; double *g2 = new double [n];
      double *solution = new double [n+2]; double *x = new double[n+2];
      // Quick setup of updated diagonal elements and value of
      solution[0] = solution[n+1] = 0.0;

      for (int i = 0; i <= n+1; i++){
          x[i]= i*h;
      }

      for (int i = 0; i < n; i++) {
          d[i] = 2;             // We get that d[0]=d_1
          g[i] = hh*f(x[i+1]);  // We get that g[0]=g_1
          d2[i]=(double)(i+2)/(i+1);
          g2[i]=g[i]/d2[i];
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
        fileout.append("alg-0-n=");
      }
      else if (atoi(argv[3])==1) {
        start2=clock();
        solution=t2(a,d,d2,solution,g,g2,n);
        finish2=clock();
        double time2=(double) (finish2-start2)/(CLOCKS_PER_SEC);
        cout <<"When n = "<< n << ", time used for specialized algorithm is "<< time2 << endl;
        fileout.append("alg-1-n=");
      }
      else if (atoi(argv[3])==2) {
          start3=clock();
          solution=t3(a,b,d,solution,g,n);
          finish3=clock();
          double time3 = (double) (finish3-start3)/(CLOCKS_PER_SEC);
          cout << "When n = "<< n << ", time used for LU algorighm is " << time3 << endl;
          fileout.append("alg-2-n=");
      }
      else{
          cout << "Please choose a algorithm to run. Choose with the third command line argument."<<endl;
          cout << "0 = General tridiagonal matrix set of linear equations."<<endl;
          cout << "1 = Specialized tridiagonal matrix set of linear equations."<<endl;
          cout << "2 = LU decomposition of linear equations."<<endl;
          exit(0);
      }
      //exact
      fileout.append(argument);
      fileout.append(".txt");

      ofile.open(fileout);
      ofile << setiosflags(ios::showpoint | ios::uppercase);
      //      ofile << "       x:             approx:          exact:       relative error" << endl;
      for (int i = 0; i < n+1;i++) {
    double xval = x[i];
     double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
         ofile << setw(20) << setprecision(8) << xval;
         ofile << setw(20) << setprecision(8) << solution[i];
         ofile << setw(20) << setprecision(8) << exact(xval);
         ofile << setw(20) << setprecision(8) << log10(RelativeError) << endl;
      }
      ofile.close();

      delete [] x; delete [] d; delete [] g; delete [] a; delete [] solution;
    }
    return 0;
}
