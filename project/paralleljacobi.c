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
  double* my_up; double* my_down; 
  double* my_left; double* my_right;
  double x,y;
  int ii, jj;
  int i_start=0, j_start=0, i_end=0, j_end=0;
  aux= malloc( N*N*sizeof(double) ); 
  memcpy(aux, grid, N*N*sizeof(double));

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
//-------------------------------------------------

  /*MPI_Sendrecv(left, N, MPI_DDOUBLE,
                 int dest, int sendtag,
                 void right, N, MPI_DOUBLE,
                 int source, int recvtag, MPI_Comm comm, MPI_Status * status);
 */
//Defining the iteration start point for p with B.Conditions
  if(pgrid.my_row==0) {j_start=1;}
  if(pgrid.my_row==pgrid.nrows-1) {j_end=N-2;}
  if(pgrid.my_col==0) {i_start=1;}
  if(pgrid.my_col==pgrid.ncols-1) {i_end=N-2;}

  //------------Jacobi step----------
  for( ii = i_start; ii <= i_end; ++ii){
    y= 1 -N*pgrid.my_row*h   -ii*h;
    for( jj = j_start; jj <= j_end; ++jj){ 
      x = N*pgrid.my_col*h +h*jj;
      //grid[ii*N + jj] = ( aux[(ii+1)*N + jj] + aux[(ii-1)*N + jj] + aux[ii*N + jj + 1] + aux[ii*N + jj - 1] + h*h*f(x,y)  )/4.0;       
    }
  }

/*
 //----------Computing delta-------SILLY VERSION!!!! (NOT a problem at the moment)
 for( ii = 1; ii <= N-2; ++ii){
   for( jj = 1; jj <= N-2; ++jj){
      delta=  ( aux[ii*N + jj] -  grid[ii*N +jj] ) / grid[ii*N+jj] ;
     if( fabs (delta ) > max  ) //We get the maximum percentual difference
       max=  fabs(delta) ;
   }
 }
*/
/* Matrix_print(aux,N,N);
 printf("\n"); 
 printf("\n"); 
 Matrix_print(grid,N,N);*/
 free(aux); free(my_up); free(my_down); free(my_left); free(my_right);
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
