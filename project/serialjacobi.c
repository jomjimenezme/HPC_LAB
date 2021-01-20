#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
//#include "get_data.c"
//#include "algebra.c"
const double h = 1/(13);
double f(double x, double y);
double g(double x, double y);
void boundary_conditions(double grid[], int N);
double jacobi(double grid[], int N);
void Matrix_print(double A[], int m, int n);

int main(int argc, char** argv){
  int i;
  int n; //n will be read and set to N=n+2;
  int N;
  double* A;

  // Initialize grid values to zero
  n=10;
  N=n+2;
  A= malloc( N*N *sizeof( double ) );
  memset( A, 0, N*N * sizeof(double) ); 
  boundary_conditions(A, N);
  
  Matrix_print(A,N,N);
 printf("\n"); 
/*while (jacobi(A,N)> 0.001){
 //printf("i ");
i++;
}*/
  //Matrix_print(A,N,N);


 jacobi(A,N);/*
  Matrix_print(A,N,N);
 printf("\n");
 printf("\n");
 jacobi(A,N);
  Matrix_print(A,N,N);
 printf("\n");
 printf("\n");
 jacobi(A,N);
  Matrix_print(A,N,N);
 printf("\n");
 printf("\n");
 printf("\n");
*/
  free(A);
}

void boundary_conditions(double grid[], int N)
{
  int ii, jj;
  
  ii = 0;
  for(jj = 0; jj < N; ++jj){
    grid[ii*N + jj] = g(ii,jj);
  }
  ii = N-1;
  for(jj = 0; jj < N; ++jj){
    grid[ii*N + jj] = g(ii,jj);
  }
  jj = 0;
  for(ii = 1; ii < N; ++ii){
    grid[ii*N + jj] = g(ii,jj);
  }
  jj = N-1;
  for(ii = 1; ii < N; ++ii){
    grid[ii*N + jj] = g(ii,jj);
  }
  
  
}

double f(double x, double y){
  return 2*( (1+x) * sin(x+y)*cos(x+y) );
}

double g(double x, double y){
  return (1+x)*sin(x+y);
}



double jacobi(double grid[], int N)
{
  double max=-1.0;
  double delta=100.0;
  double* aux;
  int ii, jj;
  aux= malloc( sizeof grid ); 
  memcpy(aux, grid, sizeof grid);
  //------------Jacobi step----------SILLY VERSION
  for( ii = 1; ii <= N-2; ++ii){
    for( jj = 1; jj <= N-2; ++jj){      
      grid[ii*N + jj] = ( grid[(ii+1)*N + jj] + grid[(ii-1)*N + jj] + grid[ii*N + jj + 1] + grid[ii*N + jj - 1] + h*h*f(ii,jj)  )/4.0; 
    }
  }
 //----------Computing delta-------
 /*for( ii = 1; ii <= N-2; ++ii){
   for( jj = 1; jj <= N-2; ++jj){
      delta=  ( aux[ii*N + jj] -  grid[ii*N +jj] ) / grid[ii*N+jj] ;
     if( fabs (delta ) > max  ) //busco el máximo de la diferencia porcentual entre la gridriz anterior y la nueva
       max=  fabs(delta) ;
   }
 }*/
 Matrix_print(aux,N,N);
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
