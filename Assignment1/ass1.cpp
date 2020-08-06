#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "mpi.h"

#define a 1664525
#define m 4294967296
#define c 1013904223
#define n0 12345
#define N 10000//10000000
#define sidelen 65536

typedef unsigned long ULONG;

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
        A = (A * a) % m ;
    }
    return A;
}

//always use argc and argv, as mpirun will pass the appropriate parms.
int main(int argc,char* argv[])
{
    int numproc, myid, namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    ULONG n_next,n_prev,number_in_circle0,number_in_circle1,ix,iy;
    double T0,T1;
    MPI_Status Stat;//status variable, so operations can be checked

    MPI_Init(&argc,&argv);//INITIALIZE
    MPI_Comm_size(MPI_COMM_WORLD, &numproc); //how many processors??
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);    //what is THIS processor-ID?

    //what is THIS processor name (hostname)?
    MPI_Get_processor_name(processor_name,&namelen);
    fprintf(stdout, "Processor ID = %d: %s %d\n", myid, processor_name, numproc);

    ULONG number_of_random_numbers = N/numproc;
    fprintf(stdout, "number_of_random_numbers = %d\n", number_of_random_numbers);  
    /*
    int k = numproc;
    ULONG A = compute_A (k);
    ULONG temp = 0 ;
    
    for (int i = k - 1  ; i > 0 ; i--){
        temp = temp + compute_A (i);
    }

    ULONG C = (c * temp) % m ;
    */
    if (myid == 0) {
        n_prev = n0;
        number_in_circle0 = 0;
        fprintf(stdout,"check pt 1 \n");
        T0 = MPI_Wtime();
        for(int i = 1; i < numproc; i++){
            //std::vector<ULONG> vbuffer;
            n_next = (a * n_prev + c) % m;
            n_prev = n_next;
            //vbuffer.push_back(n_next);
            //vbuffer.push_back(i);
            //MPI::COMM_WORLD.Send(&n_next, 1, MPI::UNSIGNED_LONG, i,0);
            fprintf(stdout, "send %ld to slave %d\n",n_next,i);
            MPI_Send(&n_next, 1, MPI_UNSIGNED_LONG, i,0, MPI_COMM_WORLD);
        }
        fprintf(stdout,"check pt 2 \n");
        n_prev = n0;
        for(int i = 0 ; i < number_of_random_numbers; i++){
            //fprintf(stdout, "master i = %d\n", i); 
            n_next = (a * n_prev + c) % m;
            n_prev = n_next;
            // Scale the random number to a random 2−d position
            ix = n_next % sidelen;
            iy = n_next / sidelen;
            // Scale current random integer to value from 0−1
            double x = rescale( ix, -1, 1);
            double y = rescale( iy, -1, 1);
            if ( is_in_circle ( x, y) ) 
                number_in_circle0++;

        }
        fprintf(stdout,"number_in_circle0 = %d \n", number_in_circle0);
        fprintf(stdout,"check pt 3 \n");
        T0 = MPI_Wtime() - T0 ;
        fprintf(stdout,"Master total time:  %f s\n", T0);
        for (int i=1;i<numproc;i++) {//receive from all nodes
            MPI_Recv(&number_in_circle1, 1, MPI_UNSIGNED_LONG, i,0, MPI_COMM_WORLD, &Stat);   
            fprintf(stdout,"number_in_circle1 = %d \n", number_in_circle1);
            number_in_circle0 += number_in_circle1;
        }
        fprintf(stdout,"final points = %d \n", number_in_circle0); 
        double result = (4 * number_in_circle0) / N ;
        fprintf(stdout,"The final result is %f \n",result);

    } 

    else {
        //std::vector<ULONG> vbuffer;
        fprintf(stdout,"slave check pt 1 \n");
        T1 = MPI_Wtime();
        MPI_Recv(&n_next, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD, &Stat);
        fprintf(stdout,"n_next = %ld \n", n_next);
        n_prev = n_next;
        number_in_circle1 = 0;

        fprintf(stdout,"slave check pt 2 \n");
        for(int i = 0 ; i < number_of_random_numbers; i++){
            //fprintf(stdout, "slave i = %d\n", i); 
            n_next = (a * n_prev + c) % m;
            n_prev = n_next;
            // Scale the random number to a random 2−d position
            ix = n_next % sidelen;
            iy = n_next / sidelen;
            // Scale current random integer to value from 0−1
            double x = rescale( ix, -1, 1);
            double y = rescale( iy, -1, 1);
            if ( is_in_circle ( x, y) ) 
                number_in_circle1++;

        }
        T1 = MPI_Wtime() - T1 ;
        fprintf(stdout,"slave check pt 3 \n");
        fprintf(stdout,"Slave total time:  %f s\n", T1);
        fprintf(stdout,"number_in_circle1 = %d \n", number_in_circle1);
        // Slave sends 'sum1' to master
        //MPI::COMM_WORLD.Send(&number_in_circle1, 1, MPI::UNSIGNED_LONG, 0,0);
        MPI_Send(&number_in_circle1, 1, MPI_UNSIGNED_LONG, 0, 0, MPI_COMM_WORLD);
    }
    MPI::Finalize();
}
