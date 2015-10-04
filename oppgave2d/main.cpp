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




void jacobi_method ( double ** A, double ** R, int n, int *numberOfIterations, int z);
double maxoffdiag ( double ** A, int * k, int * l, int n );
void rotate ( double ** A, double ** R, int k, int l, int n );
void findMini ( double ** A, int n );

ofstream ofile;

int main(){
    ofstream myfile;
    int *numberOfIterations = new int [200];
    for (int z = 5; z < 200; z++){
        //define the tridiagonal matrix
        int n = z;
        double p_max = 7.0;
        double p_min = 0.0;
        double ** A = new double*[n+1];
        double ** R = new double*[n+1];
        double h = (p_max-p_min)/(n);
        double hh = h*h;

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
        A[0][0]= hh+2/hh;;
        A[1][0]= -1/hh;
        for (int i = 1; i < n; i++){
            for(int j = 0; j < n; j++){
                if (i == j){
                    A[i][j]= (i+1)*(i+1)*hh+2/hh;
                    A[i-1][j] = -1/hh;
                    A[i+1][j] = -1/hh;
                }
            }

        }
        //last point
        A[n][n]= (n+1)*(n+1)*hh+2/hh;
        A[n-1][n] = -1/hh;
        jacobi_method ( A, R, n, numberOfIterations, z );
    }

     ofile.open("C:\\Users\\richard\\Documents\\Fys4150Project2\\oppgave2\\text.txt");
     //--------------------------------------------------------------------
     // write results to the output file
    for(int i = 5; i < 200; i++){

         ofile << setprecision(8) << i;
         ofile << setw(15) << setprecision(8) << numberOfIterations[i]<<endl;
     }
     ofile.close();

            // return a.exec();
        }

        //------------------------------(jacobi_mehtod)----------------------
        void jacobi_method ( double ** A, double ** R, int n , int *numberOfIterations, int z){

            // Setting up the eigenvector matrix
            for ( int i = 0; i < n; i++ ) {
                for ( int j = 0; j < n; j++ ) {
                    if ( i == j ) {
                        R[i][j] = 1.0;
                    } else {
                        R[i][j] = 0.0;
                    }
                }
            }
            int k, l;
            double epsilon = 1.0e-8;
            double max_number_iterations = (double) n * (double) n * (double) n;
            int iterations = 0;
            double max_offdiag = maxoffdiag ( A, &k, &l, n );

            while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ) {
                max_offdiag = maxoffdiag ( A, &k, &l, n );
                rotate ( A, R, k, l, n );
                iterations++;
            }
            numberOfIterations[z] = iterations;
            //std::cout << "Number of iterations: " << iterations << "\n";
           // findMini ( A, n );
            return;
        }

        //Function to find the manimum eigenvalue element.
        void findMini ( double ** A, int n ){

        double mini1 = 100;
        double mini2 = 100;
        double mini3 = 100;

            for ( int i = 0; i < n; i++ ) {

                    if (A[i][i] < mini1){
                        mini3 = mini2 ;
                        mini2 = mini1;
                        mini1 = A[i][i];
                    } else if (A[i][i] < mini2){
                         mini3 = mini2 ;
                         mini2 = A[i][i];
                    } else if (A[i][i] < mini3){
                         mini3 = A[i][i];
                    }
            }
            cout<<mini1<<endl;
            cout<<mini2<<endl;
            cout<<mini3<<endl;

        }

        // Function to find the maximum matrix element. Can you figure out a more
        // elegant algorithm?
        double maxoffdiag ( double ** A, int * k, int * l, int n ) {
            double max = 0.0;
            for ( int i = 0; i < n; i++ ) {
                for ( int j = i + 1; j < n; j++ ) {
                    if ( fabs(A[i][j]) > max ) {
                        max = fabs(A[i][j]);
                        *l = i;
                        *k = j;
                    }
                }
            }
            return max;
        }
        // Function to find the values of cos and sin
        void rotate ( double ** A, double ** R, int k, int l, int n ) {
            double s, c;
            if ( A[k][l] != 0.0 ) {

                double t, tau;
                tau = (A[l][l] - A[k][k])/(2*A[k][l]);
                if ( tau > 0 ) {
                    t = 1.0/(tau + sqrt(1.0 + tau*tau));
                } else {
                    t = -1.0/( -tau + sqrt(1.0 + tau*tau));
                }
                c = 1/sqrt(1+t*t);
                s = c*t;
            }
            else {
                c = 1.0;
                s = 0.0;
            }
            double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
            a_kk = A[k][k];
            a_ll = A[l][l];

            // changing the matrix elements with indices k and l
            A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
            A[l][l] = s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
            A[k][l] = 0.0; // hard-coding of the zeros
            A[l][k] = 0.0;

            // and then we change the remaining elements
            for ( int i = 0; i < n; i++ ) {

                if ( i != k && i != l ) { a_ik = A[i][k];
                    a_il = A[i][l];
                    A[i][k] = c*a_ik - s*a_il;
                    A[k][i] = A[i][k];
                    A[i][l] = c*a_il + s*a_ik;
                    A[l][i] = A[i][l];
                }

                // Finally, we compute the new eigenvectors
                r_ik = R[i][k];
                r_il = R[i][l];
                R[i][k] = c*r_ik - s*r_il;
                R[i][l] = c*r_il + s*r_ik;
            } return;
        }


