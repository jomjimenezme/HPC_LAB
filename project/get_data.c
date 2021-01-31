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



void parallel_print(char *fname, int ml, int nl, double *A, int q,  GRID_INFO_T pgrid){
  FILE *fp;
  MPI_Status status;
  int i, j, bi, bj;
  
  if(pgrid.my_rank==0){
    fp=fopen(fname, "w");
    printf("writing %ix%i matrix %s\n", ml*q, nl*q, fname);
    fprintf(fp, "%i\n", ml*q);
    fprintf(fp, "%i\n", nl*q);
  }


  for(bi=0; bi<q; bi++)
  {
    for(i=0;i<ml;i++)
    {
      for(bj=0;bj<q;bj++)
      {
	if(pgrid.my_rank==0)
	{
	  if((bi!=pgrid.my_row || (bj!=pgrid.my_col)))
	    MPI_Recv(A+i*nl, nl, MPI_DOUBLE, bj+bi*q, i, MPI_COMM_WORLD, &status);
	  for(j=0;j<nl;j++)
	  {
	    fprintf(fp, "%lf ", A[i*nl+j]);
	  }
	    //fprintf(fp,"\n");
	}else
	{
	  if((pgrid.my_row==bi)&& (pgrid.my_col==bj))
	    MPI_Send(A+i*nl, nl, MPI_DOUBLE, 0, i, MPI_COMM_WORLD);
	}	
      }
    }
  }

  if(pgrid.my_rank==0)
      fclose(fp);
}
