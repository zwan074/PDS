#include <iostream>
#include "mpi.h"
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

using namespace std;

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

int check(vector<float> data,int nitems) {
  double sum=0;
  int sorted=1;
  int i;

  for(i=0;i<nitems;i++) {
     sum+=data[i];
     if(i && data[i]<data[i-1]) sorted=0;
  }
  fprintf(stdout,"sum=%f, sorted=%d\n",sum,sorted);
}

int main(int argc, char *argv[])
{
    const int ndata = 500000;
    const float xmin = 1.0;
    const float xmax = 100.0;
    double T0,T1,T2,T3;
    MPI::Init(argc, argv);
    int numproc = MPI::COMM_WORLD.Get_size(); // number of buckets
    int myid    = MPI::COMM_WORLD.Get_rank();

    // Fill up the array with data to send to the destination node. Note
    // that the contents of the array will
    int * sendbuf_rand_nums = new int[ndata * numproc];
    int * recvbuf_rand_nums = new int[ndata];
    
    float stepsize = (xmax - xmin) / numproc;

    int root = 0;

    //generate random numbers on master proc.
    if (myid == root) {
        T0 = MPI_Wtime();
        for (int i = 0; i < ndata * numproc; ++i) {
            sendbuf_rand_nums[i] = drand48()*(xmax-xmin-1)+xmin;
        }   
        cout << "Processor " << myid << " generating random number Time: " << MPI_Wtime() - T0 << endl;
    }
    
    //divide them equally on slaves and master procs.
    T1 = MPI_Wtime();
    MPI::COMM_WORLD.Scatter(sendbuf_rand_nums, ndata, MPI_FLOAT, recvbuf_rand_nums, ndata, MPI_FLOAT, root);


    //For each proc. allocate numbers inside small buckets
    vector<vector<float> > small_bucket_2d;
    vector<float>  small_bucket_1d;
    vector<int> nitems;

    for (int i = 0; i < numproc; ++i) {
        vector<float> sub_small_bucket;
        small_bucket_2d.push_back(sub_small_bucket);
    }

    for (int i = 0; i < ndata; ++i) {
        int bktno = (int)floor((recvbuf_rand_nums[i] - xmin) / stepsize);
        small_bucket_2d[bktno].push_back(recvbuf_rand_nums[i]);
    }

    for (int i = 0; i < small_bucket_2d.size(); ++i) {
        nitems.push_back(small_bucket_2d[i].size());
        for (int j = 0; j < small_bucket_2d[i].size() ; ++j){
            small_bucket_1d.push_back(small_bucket_2d[i][j]);
        }
    }


    MPI_Barrier(MPI_COMM_WORLD);
    
    //put numbers into the big buckets


    vector<int> recvcnt(numproc);
    MPI::COMM_WORLD.Alltoall(&nitems[0], 1, MPI_INT, &recvcnt[0], 1, MPI_INT);
    vector<int> recvoff(numproc);
    recvoff[0] = 0;

    for (int n = 1; n < numproc; ++n) {
        recvoff[n] = recvoff[n-1] + recvcnt[n-1];
    }

    int big_bucket_size = 0 ;
    for (int i = 0; i < recvcnt.size() ; i++) {
        big_bucket_size += recvcnt[i];
    }

    vector<float> big_bucket( big_bucket_size );
    vector<int> sendoff(numproc);
    sendoff[0] = 0;

    for (int n = 1; n < numproc; ++n) {
        sendoff[n] = sendoff[n-1] + nitems[n-1];
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MPI::COMM_WORLD.Alltoallv(
        &small_bucket_1d[0], &nitems[0], &sendoff[0], MPI_FLOAT,
        &big_bucket[0], &recvcnt[0], &recvoff[0], MPI_FLOAT);
    
    //quick sort numbers in each big bucket
    qsort (&big_bucket[0], big_bucket.size(), sizeof(float), compare);
    
    T1 = MPI_Wtime() - T1;
    cout << "Processor " << myid << " : " << T1 << endl;

    MPI_Barrier(MPI_COMM_WORLD);

    //send them back to the master proc.
    vector<float> final_sorted_vector( ndata * numproc );
    vector<int> sendcnt_final(1);
    vector<int> recvcnt_final(numproc);
    sendcnt_final[0] = big_bucket.size();

    MPI::COMM_WORLD.Gather( &sendcnt_final[0], 1, MPI_INT, &recvcnt_final[0], 1, MPI_INT, root);

    vector<int> recvcnt_final_off(numproc);
    recvcnt_final_off[0] = 0;
    for (int n = 1; n < numproc; ++n) {
        recvcnt_final_off[n] = recvcnt_final_off[n-1] + recvcnt_final[n-1];
    }

    MPI::COMM_WORLD.Gatherv( &big_bucket[0], big_bucket.size(), MPI_FLOAT, 
                            &final_sorted_vector[0], &recvcnt_final[0], &recvcnt_final_off[0] , MPI_FLOAT, 0);

    

    //lastly compute results and check if numbers have been sorted.
    if (myid == root) {
        cout << "final vector size: " << final_sorted_vector.size() << endl;
        check(final_sorted_vector,final_sorted_vector.size() ) ;

        cout << "Processor " << myid << " Total Time : " <<  MPI_Wtime() - T0 << endl;
        //for (int i = 0; i < final_sorted_vector.size() ; i++) {
            //cout << final_sorted_vector[i] << "," ;
        //}
        //cout << endl;
    }
    

    MPI::Finalize();

}
