#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"


double jacobi(double grid[], int N, double h, GRID_INFO_T pgrid)
{
  MPI_Status  status;
  double max=-1.0;
  double delta=100.0;
  double* aux;
  double* my_up; double* my_down;
  double* my_left; double* my_right;
  double* buff_down; double* buff_up;
  double* buff_left; double* buff_right;
  double x,y;
  int ii, jj;
  int i_start=0, j_start=0, i_end=N-1, j_end=N-1;


  aux= malloc( N*N*sizeof(double) ); //copy of the local grid
  memcpy(aux, grid, N*N*sizeof(double));

  buff_down  = malloc( N*sizeof(double) ); // Buffer for msg
  buff_left  = malloc( N*sizeof(double) ); // Buffer for msg
  buff_up    = malloc( N*sizeof(double) ); // Buffer for msg
  buff_right = malloc( N*sizeof(double) ); // Buffer for msg
  //------------------Initializing boundaries------- 

  my_up= malloc( N*sizeof(double) );   my_down= malloc( N*sizeof(double) );
  my_left= malloc( N*sizeof(double) ); my_right= malloc( N*sizeof(double) );
  memcpy(my_up, grid, N*sizeof(double));
  memcpy(my_down, grid+N*(N-1), N*sizeof(double));
  for(ii=0; ii<N; ii++){
    my_left[ii] = grid[ii*N ];
  }
  for(ii=0; ii<N; ii++){
    my_right[ii] = grid[ii*N + N-1 ];
  }

 //------------------Exchange UP and DOWN  edges-------------------------------
  if(pgrid.my_row!=pgrid.nrows-1){// Not in the last grid-row 
  MPI_Sendrecv(my_down, N, MPI_DOUBLE,    pgrid.my_row+1 , 0,     buff_down, N, MPI_DOUBLE,  pgrid.my_row+1  , 1, pgrid.col_comm, &status);
  }
  if(pgrid.my_row!=0){//Not in the first grid-row
  MPI_Sendrecv(my_up, N, MPI_DOUBLE,    pgrid.my_row-1 , 1,     buff_up, N, MPI_DOUBLE,    pgrid.my_row-1, 0, pgrid.col_comm, &status);
  }

 // printf (" COL COM: my rank is %d and I received this from DOWN: %lf, and this from UP:%lf\n", pgrid.my_rank, buff_down[0], buff_up[0]);

//------------------Exchange LEFT and RIGHT edges-------------------------------
  if(pgrid.my_col!=pgrid.ncols-1){// Not in the last grid-col 
  MPI_Sendrecv(my_right, N, MPI_DOUBLE,    pgrid.my_col+1 , 2,     buff_right, N, MPI_DOUBLE,  pgrid.my_col+1  , 3, pgrid.row_comm, &status);
  }
  if(pgrid.my_col!=0){//Not in the first grid-col
  MPI_Sendrecv(my_left, N, MPI_DOUBLE,    pgrid.my_col-1 , 3,     buff_left, N, MPI_DOUBLE,    pgrid.my_col-1, 2, pgrid.row_comm, &status);
  }
  //printf (" ROW COM: my rank is %d and I received this from LEFT: %lf, and this from RIGHT:%lf\n", pgrid.my_rank, buff_left[0], buff_right[0]);
//------------Defining the grid start point for p with B.Conditions----
  if(pgrid.my_row==0) {i_start=1;}
  if(pgrid.my_row==pgrid.nrows-1) {i_end=N-2;}
  if(pgrid.my_col==0) {j_start=1;}
  if(pgrid.my_col==pgrid.ncols-1) {j_end=N-2;}

//---------------------Jacobi step------------------
//if(pgrid.my_rank==0){
  for( ii = i_start; ii <= i_end; ++ii){
    y= 1 -N*pgrid.my_row*h   -ii*h;
    for( jj = j_start; jj <= j_end; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N+jj]=h*h*f(x,y);
      if(jj+1==N){grid[ii*N + jj] += buff_right[ii];}  else{grid[ii*N+jj]+= aux[ii*N+jj+1];} //Right
      if(jj-1<0){ grid[ii*N + jj] += buff_left[ii]; }  else{grid[ii*N+jj]+= aux[ii*N+jj-1];} //Left
      if(ii+1==N){grid[ii*N + jj] += buff_down[jj]; }   else{grid[ii*N+jj]+= aux[(ii+1)*N+jj];} //Down
      if(ii-1<0){ grid[ii*N + jj] += buff_up[jj];   }     else{grid[ii*N+jj]+= aux[(ii-1)*N+jj];} //UP
      grid[ii*N+jj]/=4.0;
      delta=  fabs( aux[ii*N + jj] -  grid[ii*N +jj] ) ;

      if( delta  > max  ){  max=  delta ; /*printf("p==%d i=%d, j=%d ___ delta=%e \n", pgrid.my_rank, ii,jj,max);*/}
//      printf("i=%d, j=%d ___ x=%lf y=%lf \n", ii,jj,x,y);
    }
  }
//printf("I am process %d and i_start=%d   i_end=%d   j_start=%d  j_end=%d\n",pgrid.my_rank, i_start, i_end, j_start, j_end);
//}
//------------------FREEE MEMORY!-------------------------
  free(aux); free(my_up); free(my_down); free(my_left); free(my_right);
  free(buff_up); free(buff_down); free(buff_left); free(buff_right);
  return max;
}

