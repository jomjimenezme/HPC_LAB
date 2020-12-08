#include<stdio.h>
#include <string.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
int main(int argc, char** argv){
    
  int my_rank, p, source, dest;    
  int local_n, local_m, m,n;
  double* local_x;    
  double* local_A; 
  double result;   
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
//-----------------Allocation of local matrixes-------------------
  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
  if(m%p!=0){
    if(my_rank<m%p){
	local_m = local_m +1; }}
  local_n=n;

  local_x= malloc( local_m *sizeof( double ) );
  local_A= malloc( local_m*local_n *sizeof( double ) );
  
  Read_matrix(local_A, m, n, my_rank, p, "matrix.d");

//------------------Allocation of local vectors----------------------
  Read_matrix(local_x, m, 1, my_rank, p, "vector.d");




  result= Parallel_dot(local_x, local_x,m);
  printf("%lf\n",result);
 











if(my_rank==0){

  printf("%d %lf %lf \n", my_rank, local_A[0], local_A[15]);
} 
  free(local_x);
  free(local_A);
  MPI_Finalize();
      }
  
