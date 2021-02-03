#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "grid.c"
#include "get_data.c"
#include "boundary_and_source.c"
#include "sor.c"

void Matrix_print(double A[], int m, int n);
const double EPS= 1E-5;  //TOLERANCE

int main(int argc, char** argv){
  int i, j;
  int n; //n will be read and we will set to N=n+2;
  int N;
  double max, gmax;
  double h;
  double* local_A;
  double start, finish;
  GRID_INFO_T pgrid;

  //-------MPI initialization
  int my_rank,p;
  int ncols, nrows; // of process grid
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Setup_grid(&pgrid);
  MPI_Comm_size(pgrid.row_comm, &ncols);
  MPI_Comm_size(pgrid.col_comm, &nrows);

  // Initialize grid-related  values
  Read_N(&n, my_rank, p);  //n^2 Global Internal Points 
  N=(2+n)/nrows;           // Size of local matrix NxN  
  h = 1.0/(n+1.0);
  local_A= malloc( N*N *sizeof( double ) );
//printf("we read internals n=%d so local N=%d \n",n, N);

// for( j=1; j<10; j++){//loop for average in scaling tudy
    memset( local_A, 0.0, N*N * sizeof(double) );  
    boundary_conditions(local_A, N, h, pgrid);
  //-------Jacobi relaxation loop------
    gmax=20;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    
    //i=1;
    while (gmax>EPS){  
      max=sor(local_A,N,h,pgrid);
      MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
     // if(my_rank==0){printf("%d %e \n", i, gmax);}
      //i++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if(my_rank==0){    printf("%d %lf %.16lf \n", p, finish-start, gmax);  }
  //}

  parallel_print("matrix.d", N, N, local_A, ncols, pgrid);


///------ Freeing Memory-------------
  free(local_A);
  MPI_Finalize();
}



void Matrix_print(double A[], int m, int n){
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf("%lf ", A[i*n+j]);
    }
    printf("\n");
  }
}

/*
 * check edges 
 *
 if (pgrid.my_rank==3){
 Matrix_print(my_up, N,1);
 printf("....................\n");
 Matrix_print(my_down, N,1);
 printf("....................\n");
 Matrix_print(my_left, N,1);
 printf("....................\n");
 Matrix_print(my_right, N,1);
 printf("....................\n");
}
 
 */
