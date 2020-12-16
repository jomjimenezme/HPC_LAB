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
  double* local_C;
  double* local_A;
  double* local_B; 
  double start, finish;   
  int N=1;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status; 
for (j=0; j<N; j++){ //So that we measure the time several times
 MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  
//-----------------Allocation of local matrixes-------------------
  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
  local_n=n;

  local_A= malloc( local_m*local_n *sizeof( double ) );
  local_B= malloc( local_m*local_n *sizeof( double ) );
  //printf("size %d x %d \n", m,n);

  Read_matrix(local_A, m, n, my_rank, p, "matrix.d");
  Read_matrix(local_B, m, n, my_rank, p, "matrix.d");
  
 
//------------------Computation of C----------------------
  
  local_C[0]= serial_dot( &local_A[0] ,     &local_B[0],         );
  printf("%d, %lf, %lf, c[0]= %lf  \n", my_rank, local_B[0], local_B[5], local_C[0]);
//  for(i=0; i<local_m; i++){

 // printf("I am %d and I'll send to %d, and receive from: %d \n", my_rank, (my_rank+1)%p, ( (my_rank-1)%p +p )%(p)  );


   MPI_Sendrecv_replace(local_B, n*local_m, MPI_DOUBLE, (my_rank+1)%p, 0, ( (my_rank-1)%p +p )%p , 0, MPI_COMM_WORLD, &status );

  printf("%d, %lf, %lf, c[0]= %lf  \n", my_rank, local_B[0], local_B[5], local_C[0]);
//}
 

/*  
  x= malloc( n *sizeof( double ) );
  local_y= malloc(local_m*sizeof(double));
  Read_vector(x, n, my_rank, p);

  for(i=0; i<local_m; i++){
    local_y[i]= Serial_dot(&local_A[i*n], x , n);
  }


 MPI_Barrier(MPI_COMM_WORLD);
 finish = MPI_Wtime();


//------------------Comparing and printing------------------------------
    Read_norm(&norm, my_rank,p);
  dot=Parallel_dot(local_y, local_y, local_m);
  if(my_rank==0){
    printf("%d %lf %.16lf \n", p, finish-start,fabs(sqrt(dot)-norm));
  }

*/

  free(local_A);
  free(local_B);
}
  MPI_Finalize();
      }
  
