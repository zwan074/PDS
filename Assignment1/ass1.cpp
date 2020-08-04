#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"

#define a 1664525
#define m 4294967296
#define c 1013904223
#define n0 12345
#define N 10000000
#define sidelen 65536

typedef unsigned long ULONG;

ULONG modlin(ULONG x)
{
    return (a * x + c) % m;
}
// Put integer n in range x1 , x2 with the maximum integer value
double rescale(ULONG n, double x1, double x2)
{
    double f = static_cast<double>(n) / static_cast<double>(sidelen);
    return x1 + f * (x2 - x1);
}

bool is_in_circle (double x,double y) {
    
    return true ? x * x + y * y <= 0.25 : false;

}

ULONG compute_A (int k) {
    ULONG A = 1;
    for(int i = 0; i < k; i++){
        A = (A * i) % m ;
    }
    return A;
}

//always use argc and argv, as mpirun will pass the appropriate parms.
int main(int argc,char* argv[])
{
    MPI::Init(argc,argv);

    // What is my ID and how many processes are in this pool?
    int myid = MPI::COMM_WORLD.Get_rank();
    int numproc = MPI::COMM_WORLD.Get_size();
    std::cout << "This is id " << myid << " out of " << numproc << std::endl;
    ULONG number_of_random_numbers = N/numproc;

    if (myid == 0) {

        
        ULONG n_prev = n0;
        ULONG number_in_circle0 = 0;

        for(int i = 1; i < numproc; i++){
            //std::vector<ULONG> vbuffer;
            ULONG n_next = (a * n_prev + c) % m;
            n_prev = n_next;
            //vbuffer.push_back(n_next);
            //vbuffer.push_back(i);
            MPI::COMM_WORLD.Send(&n_next, 1, MPI::long, i,0);
        }

        int k = numproc;
        ULONG A = compute_A (k);
        ULONG temp = 0 ;
        
        for (int i = k - 1  ; i > 0 ; i--){
            temp = temp + compute_A (i);
        }

        ULONG C = (c * temp) % m ;
        n_prev = n0;

        for(int i = 0 ; i < number_of_random_numbers; i++){
           
            ULONG n_next = A * n_prev + C;
            n_prev = n_next;
            // Scale the random number to a random 2−d position
            ULONG ix = n_next % sidelen;
            ULONG iy = n_next / sidelen;
            // Scale current random integer to value from 0−1
            double x = rescale( ix, -1, 1);
            double y = rescale( iy, -1, 1);
            if ( is_in_circle ( x, y) ) 
                number_in_circle0++;

        }

        for (int i=1;i<numproc;i++) {//receive from all nodes
            ULONG number_in_circle1 = 0;
            MPI::COMM_WORLD.Recv(&number_in_circle1, 1, MPI::long, i, 0);
            number_in_circle0 += number_in_circle1;

        }
        double result = number_in_circle0 / N ;
        std::cout << "result " << result << std::endl;


    } 

    else if (myid == 1) {
        //std::vector<ULONG> vbuffer;

        ULONG n_next;
        MPI::COMM_WORLD.Recv(&n_next, 1, MPI::ULONG, 0, 0);
        ULONG n_prev = n_next;

        int k = numproc;
        ULONG number_in_circle1 = 0;
        ULONG A = compute_A (k);
        ULONG temp = 0 ;
        
        for (int i = k - 1  ; i > 0 ; i--){
            temp = temp + compute_A (i);
        }

        ULONG C = (c * temp) % m ;

        for(int i = 0 ; i < number_of_random_numbers; i++){
            ULONG n_next = A * n_prev + C;
            n_prev = n_next;
            // Scale the random number to a random 2−d position
            ULONG ix = n_next % sidelen;
            ULONG iy = n_next / sidelen;
            // Scale current random integer to value from 0−1
            double x = rescale( ix, -1, 1);
            double y = rescale( iy, -1, 1);
            if ( is_in_circle ( x, y) ) 
                number_in_circle1++;

        }


        // Slave sends 'sum1' to master
        MPI::COMM_WORLD.Send(&number_in_circle1, 1, MPI::long, 0,0);
    }
    MPI::Finalize();
}
