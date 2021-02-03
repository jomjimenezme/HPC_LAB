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



