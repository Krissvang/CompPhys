#include <iostream>
#include <armadillo>

using namespace std;

int main()
{
    cout << "Hello World!" << endl;

    unsigned int n = 10;
    arma::Mat<double> A(n, n, arma::fill::zeros);

    for(unsigned int i=0; i<n; i++){
        A(i,i) = 10;
    }

    for(unsigned int i=0; i<n-3; i++){
        A(i+3,i) = -5;
        A(i,i+3) = -5;
    }

    arma::Mat<double> P, L, U;

    arma::lu(P, L, U, A);

    cout << P << endl;

    return 0;
}
