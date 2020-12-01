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
  int n_trapez, n, counter;
  double x, h,loc_sum, integral=0;
  double a,b; 
  double local_a, local_b;
  double start, finish;
  double N=1000;
 
  MPI_Status  status;    /* return status for  receive  */                           MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p); 
  for(counter=1; counter<=N; counter++){

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  Get_data(&a,&b,&n, my_rank, p);  
  MPI_Barrier(MPI_COMM_WORLD);
  finish = MPI_Wtime();

 if(my_rank==0){	
 //printf("%d\t %.16lf\n", n, fabs(integral*4-M_PI));
 printf("%d %d %e %e \n", p, n, finish-start, fabs(integral*4-M_PI));
 }
}
  MPI_Finalize();
}


double f( double x){
double eval;
eval =1/(1+x*x);
return eval;
}
