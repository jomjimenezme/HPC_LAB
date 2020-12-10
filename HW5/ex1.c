#include<stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
int main(int argc, char** argv){
    
  int my_rank, p, i,j;    
  int local_n, local_m, m,n, dot;
  double norm;
  double* x; 
  double* local_y;
  double* local_A; 
  double* aux;
  double start, finish;   
  int N=1000;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
for (j=0; j<N; j++){ //So that we measure the time several times
   
//-----------------Allocation of local matrixes-------------------
  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
  local_n=n;

  local_A= malloc( local_m*local_n *sizeof( double ) );
  
  Read_matrix(local_A, m, n, my_rank, p, "matrix.d");

//------------------Allocation of local vectors----------------------
  
  x= malloc( n *sizeof( double ) );
  local_y= malloc(local_m*sizeof(double));
  Read_vector(x, n, my_rank, p);
MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();


  for(i=0; i<local_m; i++){
    local_y[i]= Serial_dot(&local_A[i*n], x , n);
  }

/*  for(i=0; i<local_m; i++){
 printf("rank= %d, localm= %d,  local_y %lf \n", my_rank, local_m, local_y[i]);
}*/

 MPI_Barrier(MPI_COMM_WORLD);
 finish = MPI_Wtime();


//------------------Comparing and printing------------------------------
    Read_norm(&norm, my_rank,p);
  dot=Parallel_dot(local_y, local_y, local_m);
  if(my_rank==0){
    printf("%d %lf %.16lf \n", p, finish-start,fabs(sqrt(dot)-norm));
  }



  free(local_A);
  free(local_y);
}
  MPI_Finalize();
      }
  
