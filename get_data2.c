#include <stdio.h>
#include <string.h>
#include "mpi.h"


void Get_data( float* a_ptr, float* b_ptr, int* n_ptr, int my_rank, int p){

	int source =0;
	int dest;
	int tag;

	MPI_Status status;

	if(my_rank ==0){
		printf("Enter a, b, and n\n");
		scanf("%f %f %d", a_ptr, b_ptr, n_ptr);
	}

	MPI_Bcast(a_ptr, MPI_FLOAT, 0	, MPI_COMM_WORLD);
	MPI_Bcast(b_ptr, MPI_FLOAT, 0	, MPI_COMM_WORLD);
	MPI_Bcast(n_ptr, MPI_INT, 0	, MPI_COMM_WORLD);
}

