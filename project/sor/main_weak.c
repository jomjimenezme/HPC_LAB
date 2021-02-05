#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "grid.c"
#include "get_data.c"
#include "boundary_and_source.c"
#include "sor.c"

const double EPS= 1E-5;  //TOLERANCE

int main(int argc, char** argv){
  int i, j;
  int n; //n will be read and we will set to N=n+2;
  int N;
  double max, gmax;
  double h;
  double w;
  double* local_A;
  double start, finish;
  GRID_INFO_T pgrid;
  double* aux;
  double* my_up; double* my_down;
  double* my_left; double* my_right;
  double* buff_down; double* buff_up;
  double* buff_left; double* buff_right;

  //-------MPI initialization
  int my_rank,p;
  int ncols, nrows; // of process grid
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  Setup_grid(&pgrid);
  MPI_Comm_size(pgrid.row_comm, &ncols);
  MPI_Comm_size(pgrid.col_comm, &nrows);
n=64;
 for( j=1; j<=5; j++){//loop for average in scaling tudy
  //--------------------------- Initialize grid-related  values
  N=(n)/nrows;           // Size of local matrix NxN  
  h = 1.0/(n-1.0);
  local_A= malloc( N*N *sizeof( double ) );
  w= 1.0- 4.0*M_PI/(n-1.0);
    memset( local_A, 0.0, N*N * sizeof(double) );  
    boundary_conditions(local_A, N, h, pgrid);
  
//------------------------------------Buffers initialization-----------------
    aux= malloc( N*N*sizeof(double) ); //copy of the local grid
  
    my_up= malloc( N*sizeof(double) );   my_down= malloc( N*sizeof(double) );
    my_left= malloc( N*sizeof(double) ); my_right= malloc( N*sizeof(double) );
    //buffer for message
    buff_down  = malloc( N*sizeof(double) ); buff_left  = malloc( N*sizeof(double) ); 
    buff_up    = malloc( N*sizeof(double) ); buff_right = malloc( N*sizeof(double) ); 

//------------------------------SOR relaxation loop------
    gmax=20;
    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();
    i=1;
   while (gmax>EPS){  
      initialize_buffers(local_A, my_up, my_down, my_left, my_right, N);
      max=sor(local_A,N,h,pgrid,w, aux, my_up, my_down, my_left, my_right, buff_up, buff_down, buff_left, buff_right);
      MPI_Allreduce(&max, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    //  if(my_rank==0){printf("%d %e \n", i, gmax);}
   //   i++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if(my_rank==0){    printf("%d %lf %.16lf \n", n, finish-start, gmax);  }
  n=n*2;
  }

 //parallel_print("matrix.d", N, N, local_A, ncols, pgrid);


///------ Freeing Memory-------------
  free(aux); free(my_up); free(my_down); free(my_left); free(my_right);
  free(buff_up); free(buff_down); free(buff_left); free(buff_right);
  free(local_A);
  MPI_Finalize();
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
