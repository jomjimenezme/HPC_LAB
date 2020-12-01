#include <stdio.h>
#include <string.h>
#include "mpi.h"


void Get_data( double* a_ptr, double* b_ptr, int* n_ptr, int my_rank, int p){

	int source =0;
	int dest;
	int tag;

	MPI_Status status;

	if(my_rank ==0){

    	FILE *input;
        input = fopen("input.d", "r");
        fscanf(input, "%lf", a_ptr);
        fscanf(input, "%lf", b_ptr);	
        fscanf(input, "%d", n_ptr);	
        fclose(input);

	}

	MPI_Bcast(a_ptr, 1,MPI_DOUBLE, 0	, MPI_COMM_WORLD);
	MPI_Bcast(b_ptr, 1, MPI_DOUBLE, 0	, MPI_COMM_WORLD);
	MPI_Bcast(n_ptr, 1, MPI_INT, 0	, MPI_COMM_WORLD);
}

