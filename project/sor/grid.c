#include<stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

typedef struct{
  int p; // # of processes
  MPI_Comm comm; // grid communicator
  MPI_Comm row_comm; // communicator for my row
  MPI_Comm col_comm; // communicator for my column
  int q; // order of grid
  int my_row; // my row’s coordinate
  int my_col; // my column’s coordinate
  int my_rank; // my rank in the grid
  int ncols, nrows;
} GRID_INFO_T;


void Setup_grid(GRID_INFO_T* grid){
  int wrap_around[2];
  int coordinates[2];
  int free_coords[2];
  int dimensions[2];    /* assuming p = 6  */
  int old_rank;
  /* set up global grid information */

  MPI_Comm_size(MPI_COMM_WORLD, &grid->p);
  MPI_Comm_rank(MPI_COMM_WORLD, &old_rank);

  grid->q = sqrt(grid->p);
  dimensions[0]= grid->q;
  dimensions[1]= grid->q;
 /* circular shift in second dimension, also in first */
  wrap_around[0] = 0;
  wrap_around[1] = 0;

  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, wrap_around, 1, &grid->comm);
  MPI_Comm_rank(grid->comm, &grid->my_rank);

 


/* get process coordinates in grid communicator */
  MPI_Cart_coords(grid->comm, grid->my_rank, 2, coordinates);
  grid->my_row= coordinates[0];
  grid->my_col= coordinates[1];

  /* set up row communicator */
  free_coords[0] = 0;
  free_coords[1] = 1;
  MPI_Cart_sub(grid->comm, free_coords, &grid->row_comm);


  /* set up column communicator */
  free_coords[0] = 1;
  free_coords[1] = 0;
  MPI_Cart_sub(grid->comm, free_coords, &grid->col_comm);

  MPI_Comm_size(grid->row_comm, &grid->ncols);
  MPI_Comm_size(grid->col_comm, &grid->nrows);

  //printf("my_rank = %d,out of %d, grid size= %d\n Coordinates (x,y)= (%d,%d) \n\n", old_rank, grid->p ,grid->q, grid->my_row, grid->my_col);

}




void exchange_boundaries( double my_up[], double buff_up[], double my_down[], double buff_down[], double my_left[], double buff_left[], double my_right[], double buff_right[], int N, GRID_INFO_T pgrid){
  MPI_Status  status;

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
  }



//------------------Initializing boundaries------- 
void initialize_buffers(double grid[], double my_up[], double my_down[], 
	double my_left[], double my_right[] , int N)
{
  int ii;
  memcpy(my_up, grid, N*sizeof(double));
  memcpy(my_down, grid+N*(N-1), N*sizeof(double));
  for(ii=0; ii<N; ii++){
    my_left[ii] = grid[ii*N ];
  }
  for(ii=0; ii<N; ii++){
    my_right[ii] = grid[ii*N + N-1 ];
  }

}
