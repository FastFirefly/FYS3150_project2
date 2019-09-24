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

int main(int argc, char *argv[]){
	ofstream outfile;
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

	clock_t start, finish;  //  declare start and final time for each exponent to test the time of the algorithm
    start = clock();

    string fileout = filename;
    fileout.append("_d_");

    int r_max = 8;         // r_max equal to 8 give the closest approximation for the one-electron energies
    int r_min = 0;
    double h = (double) r_max/(n+1);
    double hh = h*h;
   	
   	// Intialize variables and set values
    mat A = zeros<mat>(n,n);
    mat R = zeros<mat>(n,n);
   	double dia0 = 2.0/hh;
   	double nodia = -1.0/hh;
    double dia = dia0;
    double rho;
   	A(0,0) = dia + h;
   	R(0,0) = 1.0;
   	A(0,1) = nodia;

    // Set values along the entire matrix
   	for (int i = 1; i < n-1; i++) {
        rho = (i+1)*h;
        dia = dia0 + (rho*rho);
        A(i,i-1)  = nodia;
        A(i,i)    = dia;
        R(i,i)    = 1;
        A(i,i+1)  = nodia;
    }

    rho = rho + (n-1)*h;
    A(n-1,n-1) = dia + rho; 
    R(n-1,n-1) = 1;
    A(n-2,n-1) = nodia; 
    A(n-1,n-2) = nodia;

    vec eigact = eig_sym(A);

    double e = 1.0e-5;

    int iter = 0, p=0, q=1;
    double max = 0;
    double off = 1;
    while (off > e){
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
		a_pp = A(p,p);	// p=0, q=1
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
    // Write the results in a text file
    string fileout1 = fileout;
    fileout1.append("A");
    outfile.open(fileout1);
    outfile << setiosflags(ios::showpoint | ios::uppercase);
	for (int i = 0; i < n; i++) {
	    
	    for (int j = 0; j < n; j++) {
            /*
	    	if (A(i,j) < e) {
	    		A(i,j) = 0;
	    	}
	    */	
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

    // Unit test for eigenvalues
    
    unsigned test = false;
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            if (fabs(A(i,i) - eigact(j)) < 1.0e-5) {
                test = true;
            }
        }
        if (!test) {
            printf("Unit test failed, approximated eigenvalues are not close enough to the exact eigenvalues \n");
            printf("Exiting program...\n");
            exit(1);
        }
    }
    printf("Number of iterations = %d\n", iter); 
    /* For dimension 200, it takes 50014 itegration points to reproduce the analytical 
       results with 4 leading digits after the decimal point. Result might vary depending on
       dimension.
       For 400, it takes 206084 integration points
    */
    // Print out time used for that exponent
    finish = clock();
    printf("%f\n", ((finish - start)/(double) CLOCKS_PER_SEC ));

	return 0;
}