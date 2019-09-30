# FYS3150_project2

## To compile the programs use the following: 

\>c++ -O3 [Program] -o [Desired name for program] -larmadillo -std=c++11

## Then run the program with this: 

\>./[Desired name for program] [Name of output text file] [Desired dimension]

### Example:
Compile with: 

\>c++ -O3 2b.cpp -o 2b.x -larmadillo -std=c++11

Run with:

\>./2b.x output 200

The output you will get in the terminal is the number of iterations it took the program and how long it took to finish the code

Now you should be able to run the code :)

In the folder named results, you will find the output results for all the programs and plots of their |wavefunction|^2 in their own respective folder. 

## How to plot

First compile and run the respective c++ programs:

plot_d.py requires that 2b.cpp and 2d.cpp both have been compiled and executed

plot_e.py requires that all the varients of omega/(w1=0.01, w2=0.5, w3=1, w4=5) have been computed. 
In 2e.cpp at row 190, there will be 4 lines where 3 are commented out. To run for the different values of omega, simply uncomment the wanted omega and comment out the previous one.

After that simply run the python3 program with either:

\>plot_d.py or \>plot_e.py
