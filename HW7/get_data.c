#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

void Read_norm(double* norm_ptr , int my_rank, int p){
  if(my_rank==0){
    FILE *norm;
    norm =fopen("norm.d","r");
    fscanf(norm, "%lf", norm_ptr);
    fclose(norm);
  }
  MPI_Bcast(norm_ptr,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

void Read_size(int* m_ptr, int* n_ptr, int my_rank, int p, char* name){
  if(my_rank==0){
    FILE *size;
    size =fopen(name, "r");
    fscanf(size, "%d", m_ptr);
    fscanf(size, "%d", n_ptr);
    fclose(size);
  }
  MPI_Bcast(m_ptr,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_ptr,1, MPI_INT, 0, MPI_COMM_WORLD);
}

void Read_vector(double x[],int n, int my_rank, int p){
  int i;


  if(my_rank==0){  
    FILE *vector;
    vector =fopen("vector.d","r");
    for(i=0; i<n; i++) 
      fscanf(vector, "%lf", &x[i]);
    fclose(vector);
   }

  MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

void Read_rowblock_matrix(double local_A[], int m, int n, int my_rank, int p, char* name){
 int i,j;
  int local_m=m/p, local_n=n;
  MPI_Status status;
  double* temp; 

  if(my_rank==0){
    temp= malloc( n*m *sizeof( double)); 
    FILE *matrix;
    matrix =fopen(name,"r");
    fscanf(matrix, "%lf", &temp[0]);
    fscanf(matrix, "%lf", &temp[0]);

      for(i=0; i<m; i++){ 
        for(j=0; j<n; j++){
          fscanf(matrix, "%lf", &temp[i*local_n+j]);
        }
      }

    fclose(matrix);
  }
  MPI_Scatter(temp, local_m*local_n, MPI_DOUBLE, local_A, local_m*local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
if(my_rank==0){free(temp);};
 }

void Matrix_print(double A[], int m, int n){
  int i,j;
  for(i=0; i<m; i++){
    for(j=0; j<n; j++){
      printf("%lf\n", A[i*n+j]);
    }
  }  

}

void Read_block_matrix(double local_A[], int m, int n, int my_rank, int p, char* name){
 int i,j,c,x,y;
  int m_bar=m/sqrt(p), n_bar=n/sqrt(p);
  MPI_Status status;
  double* temp; 

  if(my_rank==0){
    temp= malloc( n*m *sizeof( double)); 
    FILE *matrix;
    matrix =fopen(name,"r");
    fscanf(matrix, "%lf", &temp[0]);
    fscanf(matrix, "%lf", &temp[0]);
      for(i=0; i<m; i++){ 
        for(j=0; j<n; j++){
          fscanf(matrix, "%lf", &temp[i*n+j]);
        }
      }

    fclose(matrix);
  

  for(c=1; c<p; c++){
    
    for(i=x*m_bar; i<x*m_bar+m_bar; i++){
      for(j=y*n_bar; j<y*n_bar+n_bar;j++){
        temp[i*n+j]=local_A[i*n+j];  
      }
    }
    Matrix_print(local_A,m_bar,n_bar);
  }
}

//  MPI_Scatter(temp, local_m*local_n, MPI_DOUBLE, local_A, local_m*local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
if(my_rank==0){free(temp);};
 }



void Parallel_blockrow_print(double local_A[], int m, int m_bar, int l, int  my_rank, int p ){
  int c;
  MPI_Status  status; 
  if(my_rank!=0){MPI_Send(local_A, m_bar*l, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD); }
  if(my_rank==0){ 
    printf("%d\n%d\n", m, l);
    Matrix_print(local_A, m_bar,l);
    for (c=1; c<p; c++){
      MPI_Recv(local_A, m_bar*l, MPI_DOUBLE, c, c, MPI_COMM_WORLD, &status);
    Matrix_print(local_A, m_bar,l);
    }
  }

}


