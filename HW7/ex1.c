#include<stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
#include "grid.c"


int main(int argc, char** argv){

  int my_rank, p, i,j,k, stage; 
  int m, n, l;   
  int source, dest;
  int nrows, ncols;
  int k_bar;
  GRID_INFO_T grid;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status;
 
  Setup_grid(&grid);
  
  Read_size(&m, &n, my_rank, p, "A.txt");
  Read_size(&n, &l, my_rank, p, "B.txt");
 
  MPI_Comm_size(grid.row_comm, &ncols);
  MPI_Comm_size(grid.col_comm, &nrows);

  dest =  ( (grid.my_row + grid.q-1)%(grid.q) ) *ncols+ grid.my_col ;
  source = ( (grid.my_row+1)%grid.q )*ncols +grid.my_col;

//  printf("my_rank = %d, Coordinates (x,y)= (%d,%d) \n sending to (%d,%d), receiving from (%d,%d) \n\n", my_rank, grid.my_row, grid.my_col, dest, dest, source,source);

  for (stage=0; stage<grid.q; stage++){
    k_bar=(grid.my_row + stage)%grid.q;
    //MPI_Broadcast(data,count,datatyoe,root,comm);
    MPI_Bcast(&my_rank,1, MPI_INT,0, grid.col_comm);
  }
// printf("I got %d\n", my_rank);


  MPI_Finalize(); 
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
