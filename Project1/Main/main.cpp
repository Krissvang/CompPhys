#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <chrono>
#include <armadillo>
// use namespace and armadillo for output and input
using namespace std;
using namespace arma;
// object for output files
ofstream ofile;

//Constants depending on Iather the program is going to performe
// certain actions. If 1 it will run, if 0 it will not.
int write_out = 0;
int relative_error_write = 0;
int take_time = 1;

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
double * t2(double* a, double *d, double *d2, double* solution, double* g,int n){
    //Saves n FLOPS to calculate e*e before the substitutions.
    // Forward substitution
    int e = a[0];
    int ee=e*e;

    //I have a special algorithm when the off diagonal
    // elements are -1 and the diagonal elements are 2.
        if(e==-1 && d[0]==2.0){
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
    //Create matrix
    Mat<double> A(n,n,fill::zeros);
    for (unsigned int i = 0; i < n; ++i) {
        A(i,i)=d[i];
    }
    for (unsigned int i = 0; i < n-1; ++i) {
        A(i,i+1)=a[i];
        A(i+1,i)=b[i];
    }
    //Create the L and U matrices
    Mat<double> L, U;
    lu(L, U, A);

    vec y(n);

    for (int i = 0; i < n; i++) {
        y(i)=g[i];

    }
    vec z = solve(L,y);
    vec x = solve(U,z);
    for (int i = 0; i < n; i++) {
        solution[i]=x[i];
    }
    /*
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
    */
    return solution;
}

// Begin main program
int main(int argc, char *argv[]){
  int exponent;
  string filename;

    // I read also the basic name for the output file and the highest poIr of 10^n I want
    if( argc <= 1 ){
          cout << "Bad Usage: " << argv[0] <<
              " read also file name on same line and max poIr 10^n" << endl;
          exit(1);
    }
        else{
        filename = argv[1]; // first command line argument after name of program
        exponent = atoi(argv[2]);
    }
    // Loop over poIrs of 10

    for (int i = 1; i <= exponent; i++){
      int  n = (int) pow(10.0,i);
      // Declare new file names
      string time_taken = "time";
      string fileout = filename;
      fileout.append("-");
      // Convert the poIr 10^i to a string
      string argument = to_string(int(pow(10,i)));


      double h = 1.0/(n+1);
      double hh = h*h;

      // Set up arrays for the program
      double *d = new double [n]; double *g = new double [n]; double *a = new double [n-1];
      double *b = new double [n-1]; double *d2 = new double [n];
      double *solution = new double [n+2]; double *x = new double[n+2];
      // Set boundary conditions
      solution[0] = solution[n+1] = 0.0;
      //Create value of position
      for (int i = 0; i <= n+1; i++){
          x[i]= i*h;
      }
      // Make "matrix" elements for our problem
      for (int i = 0; i < n; i++) {
          d[i] = 2;             // I get that d[0]=d_1 in the mathematical sense
          g[i] = hh*f(x[i+1]);  // I get that g[0]=g_1
          d2[i]=(double)(i+2)/(i+1);
      };
      for (int i=0; i < n-1; i++) {
          a[i]=-1;
          b[i]=-1;
      }
      //Runs General algorithm
      if(atoi(argv[3])==0){
          //Takes time if take_time = 1 and writes to time_taken file
          if(take_time==1){
              time_taken.append("-alg-0-n=");
              time_taken.append(argument);
              time_taken.append(".txt");
              ofile.open(time_taken);
              for (int i = 0; i < 10; ++i) {
                  auto start=chrono::steady_clock::now();
                  solution=t1(a,b,d,solution,g,n);
                  auto finish = chrono::steady_clock::now();
                  auto time=chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
                  ofile << setw(20) << setprecision(8) << time <<endl;
              }
            ofile.close();
            fileout.append("alg-0-n=");}
          else {
              solution=t1(a,b,d,solution,g,n);
              fileout.append("alg-0-n=");}

      }
        //Runs specialize algorithm
      else if (atoi(argv[3])==1) {
          //Takes time if take_time = 1 and writes to time_taken file
          if(take_time==1){
              time_taken.append("-alg-1-n=");
              time_taken.append(argument);
              time_taken.append(".txt");
              ofile.open(time_taken);
              for (int i = 0; i < 10; ++i) {
                  auto start=chrono::steady_clock::now();
                  solution=t2(a,d,d2,solution,g,n);
                  auto finish = chrono::steady_clock::now();
                  auto time=chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
                  ofile << setw(20) << setprecision(8) << time <<endl;

              }
            ofile.close();
            fileout.append("alg-1-n=");}
          else {
              solution=t2(a,d,d2,solution,g,n);
              fileout.append("alg-1-n=");}

      }

      //Runs LU decomposition
      else if (atoi(argv[3])==2) {
          //Takes time if take_time = 1 and writes to time_taken file
          if(take_time==1){
            time_taken.append("-alg-2-n=");
            time_taken.append(argument);
            time_taken.append(".txt");
            ofile.open(time_taken);
            for (int i = 0; i < 10; ++i) {
                auto start=chrono::steady_clock::now();
                solution=t3(a,b,d,solution,g,n);
                auto finish = chrono::steady_clock::now();
                auto time=chrono::duration_cast<chrono::nanoseconds>(finish - start).count();
                ofile << setw(20) << setprecision(8) << time <<endl;
            }
          ofile.close();
          fileout.append("alg-1-n=");}
        else {
            solution=t3(a,b,d,solution,g,n);
            fileout.append("alg-1-n=");}
    }

      else{
          //If no third argument this is written out in console
          cout << "Please choose a algorithm to run. Choose with the third command line argument."<<endl;
          cout << "0 = General tridiagonal matrix set of linear equations."<<endl;
          cout << "1 = Specialized tridiagonal matrix set of linear equations."<<endl;
          cout << "2 = LU decomposition of linear equations."<<endl;
          exit(0);
      }
      //exact
      //Sets filename
      fileout.append(argument);
      fileout.append(".txt");
      //Writes out the different values in a file.
      if(write_out==1){

          ofile.open(fileout);
          ofile << setiosflags(ios::showpoint | ios::uppercase);
          for (int i = 0; i < n+1;i++) {
        double xval = x[i];
         double RelativeError = fabs((exact(xval)-solution[i])/exact(xval));
             ofile << setw(20) << setprecision(8) << xval;
             ofile << setw(20) << setprecision(8) << solution[i];
             ofile << setw(20) << setprecision(8) << exact(xval);
             ofile << setw(20) << setprecision(8) << log10(RelativeError) << endl;
          }
          ofile.close();
       }
      //Writes out relative_error in a console.
      if(relative_error_write==1) {
            double MaxError=-100;
            double x_error=0;
            double* px=&x_error;
            double* p1=&MaxError;
            for (int i = 1; i < n; ++i) {
                double RelativeError = log10(fabs((exact(x[i])-solution[i])/exact(x[i])));

                double* p2=&RelativeError;
                if(MaxError < RelativeError){
                    *p1=*p2;
                    *px=x[i];


                }

            }
            cout << "Max error is: "<< MaxError << " and the x value is "<<x_error << endl;
            cout<< endl;



      }
      delete [] x; delete [] d; delete [] g; delete [] a; delete [] solution;
    }


    return 0;
}
