#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//#include "mpi.h"
//#include "get_data.c"
//#include "algebra.c"
double f(double x, double y);
double g(double x, double y);
void boundary_conditions(double grid[], int N);
double jacobi(double grid[], int N, double h);
void Matrix_print(double A[], int m, int n);
const double EPS= 1E-12;  //TOLERANCE
int main(int argc, char** argv){
  int i;
  int n; //n will be read and we will set to N=n+2;
  int N;
  double h;
  double* A;

    // Initialize grid-related  values 
  n=14;
  N=n+2;
  h=1/(n+1);
  A= malloc( N*N *sizeof( double ) );
  memset( A, 0, N*N * sizeof(double) ); 
  boundary_conditions(A, N);

//-------Jacobi relaxation loop------
  printf("\n");
  jacobi(A,N,h);

 Matrix_print(A,N,N);
  while (jacobi(A,N,h)> EPS){
  i++;
  }
  printf("%d ", i);
//------ Freeing Memory-------------
  free(A);
}

void boundary_conditions(double grid[], int N)
{
  int ii, jj;
  double x, y; 
  double h=1.0/(N-1); 
  ii = 0; //Upper  border
  y=(N-1)*h; 
  for(jj = 0; jj < N; ++jj){
    x=h*jj;
    grid[ii*N + jj] = g(x,y);
  }
  ii = N-1; //lower border
  y = 0.0;
  for(jj = 0; jj < N; ++jj){
    x=h*jj;
    grid[ii*N + jj] = g(x,y);
  }
  jj = 0; //left border
  x=0.0;
  for(ii = 0; ii < N; ++ii){
    y=(N-1)*h-ii*h;
    grid[ii*N + jj] = g(x,y);
  }
  jj = N-1; //right border
  x = h*jj;
  for(ii = 1; ii < N; ++ii){
    y=(N-1)*h-h*ii;
    grid[ii*N + jj] = g(x,y);
  }
  
}

double f(double x, double y){
  return 2*( (1+x) * sin(x+y)*cos(x+y) );
}

double g(double x, double y){
  return (1+x)*sin(x+y);
  //return x;
 // return y;
}



double jacobi(double grid[], int N, double h)
{
  double max=-1.0;
  double delta=100.0;
  double x=0.0; double y=1.0;
  double* aux;
  int ii, jj;
  aux= malloc( N*N*sizeof(double) ); 
  memcpy(aux, grid, N*N*sizeof(double));
  //------------Jacobi step----------
  for( ii = 1; ii <= N-2; ++ii){
    for( jj = 1; jj <= N-2; ++jj){      
      grid[ii*N + jj] = ( aux[(ii+1)*N + jj] + aux[(ii-1)*N + jj] + aux[ii*N + jj + 1] + aux[ii*N + jj - 1] + h*h*f(x,y)  )/4.0; 
      x+= h;
    }
    y-=h;
  }
 //----------Computing delta-------SILLY VERSION!!!! (NOT a problem at the moment)
 for( ii = 1; ii <= N-2; ++ii){
   for( jj = 1; jj <= N-2; ++jj){
      delta=  ( aux[ii*N + jj] -  grid[ii*N +jj] ) / grid[ii*N+jj] ;
     if( fabs (delta ) > max  ) //We get the maximum percentual difference
       max=  fabs(delta) ;
   }
 }
// Matrix_print(aux,N,N);
 //printf("\n"); 
// printf("\n"); 
 //Matrix_print(grid,N,N);
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

