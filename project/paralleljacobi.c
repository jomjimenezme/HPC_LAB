#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "grid.c"
//#include "get_data.c"
//#include "algebra.c"
double f(double x, double y);
double g(double x, double y);
void boundary_conditions(double grid[], int N, double h, GRID_INFO_T pgrid);
double jacobi(double grid[], int N, double h, GRID_INFO_T pgrid);
void Matrix_print(double A[], int m, int n);
void sillyp(double local_A[], int my_rank, int N, GRID_INFO_T pgrid);
const double EPS= 1E-12;  //TOLERANCE

int main(int argc, char** argv){
  int i,j;
  int n; //n will be read and we will set to N=n+2;
  int N;
  double h;
  double* local_A;
  GRID_INFO_T pgrid;

  //-------MPI initialization
  int my_rank,p;
  int ncols, nrows; // of process grid
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status; 

  // Initialize grid-related  values 
  N = (32)/p;  //Size of local grid NxN
  n = ( (N-1)*(N-1)*p ); //Global Internal Points
  h = 1.0/(N*p/2-1);
  local_A= malloc( N*N *sizeof( double ) );
  memset( local_A, 0, N*N * sizeof(double) ); 
  Setup_grid(&pgrid);
  MPI_Comm_size(pgrid.row_comm, &ncols);
  MPI_Comm_size(pgrid.col_comm, &nrows);

  //rank, grid order, row, col, grid rank
//  printf("%d %d %d %d %d \n", my_rank, pgrid.q, pgrid.my_row, pgrid.my_col, pgrid.my_rank);
  boundary_conditions(local_A, N, h, pgrid);
  sillyp(local_A, my_rank, N, pgrid);

  jacobi(local_A,N,h,pgrid);
/*
//-------Jacobi relaxation loop------
  Matrix_print(local_A,N,N);
  printf("\n"); 
  while (jacobi(local_A,N,h)> EPS){
   //printf("i ");
  i++;
  }
*/
//------ Freeing Memory-------------
  free(local_A);
  MPI_Finalize();
}

void boundary_conditions(double grid[], int N, double h, GRID_INFO_T pgrid)
{
  int ii, jj;
  double x, y;
  if (pgrid.my_row==0){
    ii = 0; //upper border
    y=1.0;//(N-1)*h; 
    for(jj = 0; jj < N; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if (pgrid.my_row==pgrid.nrows-1){
    ii = N-1; //lower border
    y = 0.0;
    for(jj = 0; jj < N; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if (pgrid.my_col==0){
    jj = 0;  //left border
    x=0.0;
    for(ii = 0; ii < N; ++ii){
      y= 1 -N*pgrid.my_row*h   -ii*h;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if(pgrid.my_col==pgrid.ncols-1){
    jj = N-1;//right border
    x = 1.0;
    for(ii = 1; ii < N; ++ii){
      y= 1 -N*pgrid.my_row*h   -ii*h;
      grid[ii*N + jj] = g(x,y);
    }
  }  
}

double f(double x, double y){
  return 2*( (1+x) * sin(x+y)*cos(x+y) );
}

double g(double x, double y){
  return (1+x)*sin(x+y);
}



double jacobi(double grid[], int N, double h, GRID_INFO_T pgrid)
{
  double max=-1.0;
  double delta=100.0;
  double* aux;
  double x,y;
  int ii, jj;
  aux= malloc( N*N*sizeof(double) ); 
  memcpy(aux, grid, N*N*sizeof(double));
  //------------Jacobi step----------
  for( ii = 1; ii <= N-2; ++ii){
    y= 1 -N*pgrid.my_row*h   -ii*h;
    for( jj = 1; jj <= N-2; ++jj){ 
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N + jj] = ( grid[(ii+1)*N + jj] + grid[(ii-1)*N + jj] + grid[ii*N + jj + 1] + grid[ii*N + jj - 1] + h*h*f(ii,jj)  )/4.0; 
    }
  }
 //----------Computing delta-------SILLY VERSION!!!! (NOT a problem at the moment)
 for( ii = 1; ii <= N-2; ++ii){
   for( jj = 1; jj <= N-2; ++jj){
      delta=  ( aux[ii*N + jj] -  grid[ii*N +jj] ) / grid[ii*N+jj] ;
     if( fabs (delta ) > max  ) //We get the maximum percentual difference
       max=  fabs(delta) ;
   }
 }

/* Matrix_print(aux,N,N);
 printf("\n"); 
 printf("\n"); 
 Matrix_print(grid,N,N);*/
 free(aux);
 return max;
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


void sillyp(double local_A[], int my_rank, int N, GRID_INFO_T pgrid){
int i,j, ncols,nrows;
MPI_Comm_size(pgrid.row_comm, &ncols);
  MPI_Comm_size(pgrid.col_comm, &nrows);

 for(i=0; i<nrows; i++){
    for(j=0; j<ncols; j++){
      sleep(1);
     MPI_Barrier( MPI_COMM_WORLD );

      if( pgrid.my_col==j && pgrid.my_col==i ){
        printf("%d, (%d,%d), %d  \n", my_rank, pgrid.my_row, pgrid.my_col, pgrid.my_rank);
        Matrix_print(local_A, N,N);
        printf("\n");
      }
     sleep(1);
     MPI_Barrier( MPI_COMM_WORLD );
         }
           }
    }
  
/*    
  int my_rank, p, i,j,k, c;    
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
