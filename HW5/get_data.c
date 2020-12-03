#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"


void Read_size(int* m_ptr, int* n_ptr, int my_rank, int p){
  int i;
  int ncols =0, nrows=0;
  char c;

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


