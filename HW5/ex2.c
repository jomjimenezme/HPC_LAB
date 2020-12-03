#include<stdio.h>
#include <string.h>
#include "mpi.h"
#include "get_data.c"
int main(int argc, char** argv){
    
  int my_rank, p, source, dest;    
  int local_n, local_m, m,n;
  double* local_x;    
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
  local_n=n/p;
  local_x= malloc( local_n *sizeof( double ) );
  Read_vector(local_x, local_n, my_rank, p);
  

  
 // printf("Hello World, I am proces %d out of %d processes     %d  %d\n", my_rank, p,m, n); 
   printf("Hello World, I am proces %d out of %d processes     %lf  %lf\n", my_rank, p,local_x[0], local_x[1]); 
  free(local_x);

  
  MPI_Finalize();
      }
  
