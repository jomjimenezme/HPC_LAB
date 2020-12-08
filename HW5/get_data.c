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


void Read_matrix(double local_A[], int m, int n, int my_rank, int p, char* name){
 int i,j;
  int local_m, local_n=n;
  int dest=0;
  MPI_Status status;
    
  if(my_rank==0){  
    FILE *matrix;
    matrix =fopen(name,"r");
    while(dest<p){
    local_m=m/p;
      if(m%p!=0){
	if(dest<m%p){
	  local_m = local_m+1;
	}
      }
      for(i=0; i<local_m; i++){ 
        for(j=0; j<local_n; j++){
          fscanf(matrix, "%lf", &local_A[i*local_n+j]);
        }
		//printf("%lf, %lf, %lf, %lf \n", local_A[0], local_A[1], local_A[2], local_A[3]);

      }

          MPI_Send(local_A,local_m*local_n, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD);

	  dest++; 
//printf("Heeey rank=%d  m=%d  localm=%d \n", my_rank,m,local_m);
    }
    fclose(matrix);
  }
 local_m=m/p;
      if(m%p!=0){
	if(my_rank<m%p){
	  local_m = local_m+1;
	}
      }
 
//printf("Heeey rank=%d  m=%d  localm=%d \n", my_rank,m,local_m);
    MPI_Recv(local_A, local_m*local_n, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);


 }


/*void Read_vector(double local_A[], int m, int n, int my_rank, int p){

  int i,j;
  int local_m, local_n=n;
  int dest=0;
  MPI_Status status;  
  local_m=m/p;
      if(m%p!=0){
	if(dest<m%p){
	  local_m = local_m+1;
	}
      }

  if(my_rank==0){  
    FILE *matrix;
    matrix =fopen("vector.d","r");
    while(dest<p){
    local_m=m/p;
      if(m%p!=0){
	if(dest<m%p){
	  local_m = local_m+1;
	}
      }
      for(i=0; i<local_m; i++){ 
        for(j=0; j<local_n; j++){
          fscanf(matrix, "%lf", &local_A[i*local_m+j]);
        }
          MPI_Send(local_A,local_m*local_n, MPI_DOUBLE, dest, dest, MPI_COMM_WORLD);
	  dest++; 
      }
    }
    fclose(matrix);
  }  

    MPI_Recv(local_A, local_m*local_n, MPI_DOUBLE, 0, my_rank, MPI_COMM_WORLD, &status);


}
*/



/*
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
}*/
