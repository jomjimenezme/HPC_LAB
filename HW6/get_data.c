#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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


void Read_matrix(double local_A[], int m, int n, int my_rank, int p, char* name){
 int i,j;
  int local_m=m/p, local_n=n;
  MPI_Status status;
  double* temp; 

  if(my_rank==0){
    temp= malloc( n*m *sizeof( double)); 
    FILE *matrix;
    matrix =fopen(name,"r");
    fscanf(matrix, "%d", &temp[0]);
    fscanf(matrix, "%d", &temp[0]);

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






