#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"


double f(double x);

int main(int argc, char** argv){
  //MPI Variables 
  int 	my_rank, p, i, j;
  
  //Integral Variables
  int n_trapez, n, counter,a0,b0,n0;
  double x, h,loc_sum, integral;
  double a,b; 
  double local_a, local_b;
  double start, startl, finishl;
  double N=1000; //Averaging over N iterations. 

  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  //  printf("Time used for the funciton get data");
  for(counter=1; counter<=N; counter++){
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  Get_data(&a,&b,&n, my_rank, p);  
  MPI_Barrier(MPI_COMM_WORLD);
  finishl = MPI_Wtime();
 
  if(my_rank==0){
    printf("%d %d %e \n", p, n, finishl-start);
}
 }
 
  MPI_Finalize();
}


double f( double x){
double eval;
 eval =1/(1+x*x);
 return eval;
}
