#include<stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
int main(int argc, char** argv){
    
  int my_rank, p, i,j,k, c;    
  int local_n, local_m, m,n;
  double norm;
  int shift=0;
  double* local_C;
  double* local_A;
  double* local_B; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status; 

  
//-----------------Allocation of local matrixes-------------------
  Read_size(&m, &n, my_rank, p);
  local_m=m/p;
  local_n=n;

  local_A= malloc( local_m*local_n *sizeof( double ) );
  local_B= malloc( local_m*local_n *sizeof( double ) );
  local_C= malloc( local_m*local_n *sizeof( double ) );

  Read_matrix(local_A, m, n, my_rank, p, "matrix.d");
  Read_matrix(local_B, m, n, my_rank, p, "matrix.d");
  
 
//------------------Computation of C----------------------
  

  for(c=0; c<=local_m; c++){ //Loop for circular shift.
    shift= local_m * (c +  my_rank)%p  ; //start position for A after every stage of rotation

    for (i=0; i<local_m; i++){ // The three loops of Matrix-Matrix product
      for (j=0; j<local_n; j++){
	for(k=0; k< local_m; k++){
          local_C[  i*local_n +j ] += local_A[ ( shift ) +  i*local_n+k  ]  *  local_B[ k * local_n +j  ];
	}
      }
    }
    MPI_Sendrecv_replace(local_B, n*local_m, MPI_DOUBLE,   ( (my_rank-1)%p +p )%p  , 0, (my_rank+1)%p , 0, MPI_COMM_WORLD, &status );
  }
 


//------------------Silly printing------------------------------
  for (c=0; c<p; c++){
    if(my_rank==c){
      for(i=0; i<local_m; i++){
        for(j=0; j<local_n; j++){
  	  printf(" %lf ", local_C[i*local_n+j]);
        }
  	  printf("\n");
      }
    }
  }


  free(local_A);
  free(local_B);
  free(local_C);
  MPI_Finalize();
      
}
  
