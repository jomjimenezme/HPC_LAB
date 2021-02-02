#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mpi.h"

void Read_N(int* N , int my_rank, int p){
  if(my_rank==0){
    FILE *fp;
    fp =fopen("N.d","r");
    fscanf(fp, "%d", N);
    fclose(fp);
  }
  MPI_Bcast(N,1,MPI_INT, 0, MPI_COMM_WORLD);
  }




void parallel_print(char *fname, int ml, int nl, double *A, int q,  GRID_INFO_T pgrid){
  FILE *fp;
  MPI_Status status;
  int i, j, bi, bj;
  
  if(pgrid.my_rank==0){
    fp=fopen(fname, "w");
//    printf("writing %ix%i matrix %s\n", ml*q, nl*q, fname);
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
