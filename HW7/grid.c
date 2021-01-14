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
  wrap_around[0] = 1;
  wrap_around[1] = 1;

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

 
  //printf("my_rank = %d,out of %d, grid size= %d\n Coordinates (x,y)= (%d,%d) \n\n", old_rank, grid->p ,grid->q, grid->my_row, grid->my_col);

}




/*    
  int m_bar, n_bar;
  int m, n, l;
  int source, dest;
  int shift=0;
  double* local_C;
  double* local_A;
  double* local_B; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status; 

//-----------------Allocation of local matrixes-------------------
  
  Read_size(&m, &n, my_rank, p, "A.txt");
  Read_size(&n, &l, my_rank, p, "B.txt");
  m_bar=m/p;
  n_bar=n/p;

  local_A= malloc( m_bar*n *sizeof( double ) );
  local_B= malloc( n_bar*l *sizeof( double ) );
  local_C= malloc( m_bar*l *sizeof( double ) );

  Read_matrix(local_A, m, n, my_rank, p, "A.txt");
  Read_matrix(local_B, n, l, my_rank, p, "B.txt");
  
  memset( local_C, 0, m_bar *l * sizeof(double) ); 
//------------------Computation of C----------------------
  source = (my_rank+1)%p ;
  dest = ( my_rank-1+p )%p;

  for(c=0; c<p; c++){ //Loop for circular shift of B.
    shift= (c +  my_rank)%p  ; //start position for A after every stage of rotation

   if(c)  {MPI_Sendrecv_replace(local_B, n_bar*l, MPI_DOUBLE,  dest  , 0, source , 0, MPI_COMM_WORLD, &status );}
    for (i=0; i<m_bar; i++){ // The three loops of Matrix-Matrix product
      for (j=0; j<l; j++){
	for(k=0; k< n_bar; k++){
          local_C[  i*l +j ] += local_A[  i*n + k  + shift*n_bar ]  *  local_B[ k*l +j  ];
	}
      }
    }
   MPI_Barrier(MPI_COMM_WORLD);
   //  MPI_Sendrecv_replace(local_B, n_bar*l, MPI_DOUBLE,  dest  , 0, source , 0, MPI_COMM_WORLD, &status );
  }

//----------------- Printing------------------------------
Parallel_blockrow_print(local_C, m, m_bar, l, my_rank, p);

  free(local_A);
  free(local_B);
  free(local_C);
  MPI_Finalize();
      
}
 */ 
