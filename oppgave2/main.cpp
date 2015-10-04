/*
Jacobi's method for finding eigenvalues
 eigenvectors of the symetric matrix A.

The eigenvalues of A will be on the diagonal
 of A, with eigenvalue i being A[i][i].
 The j-th component of the i-th eigenvector
 is stored in R[i][j].

A: input matrix (n x n)
R: empty matrix for eigenvectors (n x n)
n: dimention of matrices
*/
#include <iomanip>
#include <QCoreApplication>
#include <iostream>
#include <cmath>
#include <armadillo>
#include <fstream>


using namespace std;
using namespace arma;



ofstream ofile;

int main(){
    ofstream myfile;

        //define the tridiagonal matrix
        int n = 500;
        double p_max = 15.0;
        double p_min = 0.0;
        double ** A = new double*[n+1];
        double ** R = new double*[n+1];
        double h = (p_max-p_min)/(n);
        double hh = h*h;
        double wr = 0.1;//0.5, 1, 5

        for (int i = 0; i < n+1; i++){
            A[i] = new double[n+1];
            R[i] = new double[n+1];
        }

        //make sure every thing is zero
        for (int i = 0; i < n+1; i++){
            for(int j = 0; j < n+1; j++){
                R[i][j]=0;
                A[i][j]=0;
            }
        }

        //first point
        A[0][0]= wr*wr*hh+1/h+2/hh;;
        A[1][0]= -1/hh;
        for (int i = 1; i < n; i++){
            for(int j = 0; j < n; j++){
                if (i == j){
                    A[i][j]= wr*wr*(i+1)*(i+1)*hh+1/(h*(i+1))+2/hh;
                    A[i-1][j] = -1/hh;
                    A[i+1][j] = -1/hh;
                }
            }

        }
        //last point
        A[n][n]= wr*wr*(n+1)*(n+1)*hh+1/((n+1)*h)+2/hh;
        A[n-1][n] = -1/hh;

        //armidilo method ----------------
        mat A1 = zeros<mat>(n+1,n+1);
        for (int i = 0; i < n+1; i++){
            for(int j = 0; j < n+1; j++){
                A1(i,j) = A[i][j];
            }
            }
                vec eigenval;
                mat eigenvec;

                eig_sym(eigenval,eigenvec,A1);

                cout << eigenval(0) << endl;
                cout << eigenval(1) << endl;
                cout << eigenval(2) << endl;


     ofile.open("C:\\Users\\richard\\Documents\\Fys4150Project2\\oppgave2d\\text.txt");
     //--------------------------------------------------------------------
     // write results to the output file
    for(int i = 0; i < n; i++){

         ofile << setprecision(8) << i*h;
         ofile << setw(15) << setprecision(8) << eigenvec(i,0)*eigenvec(i,0);
         ofile << setw(15) << setprecision(8) << eigenvec(i,1)*eigenvec(i,1)<<endl;
     }
     ofile.close();

            // return a.exec();
        }


