#include<stdio.h>
#include <string.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
int main(int argc, char** argv){
    
  int my_rank, p, i;    
  int local_n, local_m, m,n;
  double* x; 
  double* local_y;   
  double* local_A; 
  double start, finish;   
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();
 
//-----------------Allocation of local matrixes-------------------
  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
/*  if(m%p!=0){
    if(my_rank<m%p){
	local_m = local_m +1; }}
*/
  local_n=n;

  local_A= malloc( local_m*local_n *sizeof( double ) );
  
  Read_matrix(local_A, m, n, my_rank, p, "matrix.d");

//------------------Allocation of local vectors----------------------
  
  x= malloc( n *sizeof( double ) );
  local_y= malloc(local_m*sizeof(double));
  Read_vector(x, n, my_rank, p);


  for(i=0; i<local_m; i++){
    local_y[i]= Serial_dot(&local_A[i*n], x , n);
  }

 /* for(i=0; i<local_m; i++){
 printf("rank= %d, localm= %d,  dot %lf \n", my_rank, local_m, local_y[i]);
}*/

 MPI_Barrier(MPI_COMM_WORLD);
 finish = MPI_Wtime();

  if(my_rank==0){
    printf("%d %lf \n", p, finish-start);
  }
 
  free(x);
  free(local_A);
  free(local_y);
  MPI_Finalize();
      }
  
