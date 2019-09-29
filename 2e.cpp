#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <string>
#include <armadillo>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
using namespace std;
using namespace arma;

int find_max(int n, mat& A, int& p, int& q, double& max) {
    // Finds the maximum non diagonal element in the matrix
    for (int i = 0; i < n; i++) {
        for (int j = i+1; j < n; j++) {
            double aij = fabs(A(i,j));
            if (aij > max) {
                 max = aij; 
                 p = i; 
                 q = j;
            }
        }
    }
    return 0;
}

int jacobi(int n, mat& A, mat& R, double e, int maxiter) {
    // Finds the approximated eigenvalues and eigenvectors through jacobi's method

    int iter = 0, p=0, q=1;
    double max = 0;
    double off = 1;

    while (off > e || iter<=maxiter){
        find_max(n, A, p, q, max);
        
        off = max;
        max = 0;

        // Find Jacobi rotation
        double s, c;
        double t;

        if (A(p,q) != 0.0) {

            double tau = (A(q,q) - A(p,p))/(2*A(p,q));
        
            if (tau >= 0) {
                t = 1.0/(tau + sqrt(1.0 + tau*tau));
            }
            else {
                t = -1.0/(-tau + sqrt(1.0 + tau*tau));
            }

            c = 1/sqrt(1+t*t);
            s = c*t;
        }
        else {
            c = 1.0;
            s = 0.0;
        }

        double a_pp, a_qq, a_ip, a_iq, r_ip, r_iq;
        a_pp = A(p,p);  // p=0, q=1
        a_qq = A(q,q);
        A(p,p) = c*c*a_pp - 2.0*c*s*A(p,q) + s*s*a_qq;
        A(q,q) = s*s*a_pp + 2.0*c*s*A(p,q) + c*c*a_qq;
        A(p,q) = 0.0; // hard-coding non-diagonal elements by hand
        A(q,p) = 0.0; 
        for ( int i = 0; i < n; i++ ) {
            if ( i != p && i != q ) {
                a_ip = A(i,p);
                a_iq = A(i,q);
                A(i,p) = c*a_ip - s*a_iq;
                A(p,i) = A(i,p);
                A(i,q) = c*a_iq + s*a_ip;
                A(q,i) = A(i,q);
            }
            // And finally the new eigenvectors
            r_ip = R(i,p);
            r_iq = R(i,q);
            R(i,p) = c*r_ip - s*r_iq;
            R(i,q) = c*r_iq + s*r_ip;
        }
        iter++;
    }
    return 0;
}

int set_matrix(mat& A, mat& R, int n, int r_max, double omega) {
    // Sets all the values of the matrix

    double h = (double) r_max/(n+1);
    double hh = h*h;
    
    // Intialize variables and set values for the first elements
    
    double dia0 = 2.0/hh;
    double nodia = -1.0/hh;
    double dia = dia0;
    double rho;

    double omega2 = omega*omega;
    A(0,0) = dia + (hh*omega2)+(1/h);
    R(0,0) = 1.0;
    A(0,1) = nodia;

    // Set values along the entire matrix
    for (int i = 1; i < n-1; i++) {
        rho = (i+1)*h;
        dia = dia0 + (rho*rho*omega2) + (1/rho);
        A(i,i-1)  = nodia;
        A(i,i)    = dia;
        R(i,i)    = 1;
        A(i,i+1)  = nodia;
    }

    rho = rho + (n-1)*h;
    A(n-1,n-1) = dia + (rho*rho*omega2) + (1/rho); 
    R(n-1,n-1) = 1;
    A(n-2,n-1) = nodia; 
    A(n-1,n-2) = nodia;

    return 0;
}

int write_file(int n, mat A, mat R, string fileout, vec eigact) {
    // Writes out the eigenvalues, eigenvectors and exact values found from the calculations
    
    ofstream outfile;
    string fileout1 = fileout;
    fileout1.append("A");
    outfile.open(fileout1);
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < n; i++) {
        
        for (int j = 0; j < n; j++) {
            outfile << setw(15) << setprecision(8) << A(i,j);
        }
        outfile << endl;
    }

    outfile.close();

    string fileout2 = fileout;
    fileout2.append("R");
    outfile.open(fileout2);
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            outfile << setw(15) << setprecision(8) << R(i,j);
        }
        outfile << endl;
    }

    outfile.close();
    

    string fileoutX = fileout;
    fileoutX.append("Exact");
    outfile.open(fileoutX);
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    for (int i = 0; i < n; i++) {
        outfile << setw(15) << setprecision(8) << eigact(i) << endl;
    }

    outfile.close();
}

int main(int argc, char *argv[]){
	string filename;
	int n;
	if(argc <= 1){
    	printf("Error: Please include filename and n dimension of matrix\n");
    	exit(1);
    }
    else{
        filename = argv[1]; 
        n = atoi(argv[2]);
    }

    string fileout = filename;
    fileout.append("_e_");

    int r_max = 10;

    // Set the different values for the different frequencies, for the different values, uncomment them
    double omega = 0.01; fileout.append("w1_");
    // double omega = 0.5; fileout.append("w2_");
    // double omega = 1.0; fileout.append("w3_");
    // double omega = 5.0; fileout.append("w4_");

    mat A = zeros<mat>(n,n);
    mat R = zeros<mat>(n,n);

    set_matrix(A, R, n, r_max, omega);

    // Unit test to check if the algorithm has found the largest non-diagonal element in the matrix
    double max = 0;
    int p = 0, q = 1;
    find_max(n, A, p, q, max);
    unsigned test1 = true;     // Test value for the unit test to check if our algorithm found the largest non-diagonal element
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i!=j) {
                if (max<fabs(A(i,j))) {
                    test1 = false;
                }
            }
        }
    }
    if (!test1) {
        printf("Unit test 1 failed, algorithm did not find the largest non-diagonal element in the matrix \n");
        printf("Exiting program...\n");
        exit(1);
    }

    clock_t start1, start2, finish1, finish2;  //  declare start and final time for each exponent to test the time of the algorithm

    start1 = clock();
    vec eigact = eig_sym(A);        // Find the exact values of eigenvalues
    finish1 = clock();
    printf("For Armadillo's exact: %f\n", ((finish1 - start1)/(double) CLOCKS_PER_SEC ));

    double e = 1.0e-15;         // Threshold that is required by the largest non-diagonal element to be below
    int maxiter = 10;

    start2 = clock();
    jacobi(n, A, R, e, maxiter);    // Calculates the 
    finish2 = clock();
    printf("For the jacobi %f\n", ((finish2 - start2)/(double) CLOCKS_PER_SEC ));

    // Write the results in a text file
    write_file(n, A, R, fileout, eigact);

    // Unit test to check if the approximated eigenvalues are close enough to the exact values
    unsigned test2 = false;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (fabs(A(i,i) - eigact(j)) < 1.0e-5) {
                test2 = true;
            }
        }
        if (!test2) {
            printf("Unit test 2 failed, approximated eigenvalues are not close enough to the exact eigenvalues \n");
            printf("Exiting program...\n");
            exit(1);
        }
    }

    // Unit test to check if all eigenvectors are orthogonal
    unsigned test3 = true;
    mat B= zeros<mat>(n,n); B.randu();
    for (int i = 0; i < n; i++) {
        if ((dot(R.col(i), R.col(i))-1) > 1.0e-5) {
            test3 = false;
        }
    }
    if (!test3) {
        printf("Unit test 3 failed, approximated eigenvectors are not orthogonal \n");
        printf("Exiting program...\n");
        exit(1);
    }
	return 0;
}