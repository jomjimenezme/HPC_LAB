#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"


void Read_size(int* m_ptr, int* n_ptr, int my_rank, int p){
  if(my_rank==0){
    FILE *size;
    size =fopen("size.d","r");
    fscanf(size, "%d", m_ptr);
    fscanf(size, "%d", n_ptr);
    fclose(size);
  }
  MPI_Bcast(m_ptr,1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_ptr,1, MPI_INT, 0, MPI_COMM_WORLD);

}


void Read_vector(double local_x[],int local_n, int my_rank, int p){
  int i;
  int n=local_n*p;
  //Dynamically allocation of vector
  double* temp;
  temp = malloc( n * sizeof(double) ) ;

  if(my_rank==0){  
    FILE *vector;
    vector =fopen("vector.d","r");
    for(i=0; i<p*local_n; i++) 
      fscanf(vector, "%lf", &temp[i]);
    fclose(vector);
   }

  MPI_Scatter(temp, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE,0, MPI_COMM_WORLD);

  free(temp);
}

void Read_matrix(double local_A[], int local_m, int local_n, int my_rank, int p){

  int i,j;
  //int n=local_n*p;
  //int m= local_m*p;
  int dest=0;
  //MPI_Request request, request2;
  MPI_Status status;
  //Dynamically allocation of vector
  //temp = malloc( local_n*local_m * sizeof(double) ) ;

  if(my_rank==0){  
    FILE *matrix;
    matrix =fopen("matrix.d","r");
    while(dest<p){
      for(i=0; i<local_m; i++){ 
        for(j=0; j<local_n; j++){
          fscanf(matrix, "%lf", &local_A[i+j]);
        }
        //if(local_m%j==0){
		printf("%lf, %lf, %lf, %lf \n", local_A[0], local_A[1], local_A[2], local_A[3]);
          MPI_Send(local_A,local_m*local_n, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD);
          //MPI_Isend(&local_A, local_m*local_n, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD, &request); 
	  dest++; 
        //}
      }
        //MPI_Irecv(&local_A, local_m*local_n, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD, &request2); 
    }
    fclose(matrix);
  }  

    MPI_Recv(local_A, local_m*local_n, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);

}
