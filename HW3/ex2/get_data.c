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


		for (dest=1; dest<p; dest++){
			tag=0;
			MPI_Send(a_ptr, 1 , MPI_DOUBLE , dest , tag , MPI_COMM_WORLD);
	
			tag =1;
			MPI_Send(b_ptr, 1 , MPI_DOUBLE , dest , tag , MPI_COMM_WORLD);
		
			tag =2;
			MPI_Send(n_ptr, 1 , MPI_INT , dest, tag, MPI_COMM_WORLD);
		}
	}else{
		tag = 0;
		MPI_Recv(a_ptr, 1 , MPI_DOUBLE , source , tag , MPI_COMM_WORLD, &status);
                tag =1;
		MPI_Recv(b_ptr, 1 , MPI_DOUBLE , source , tag , MPI_COMM_WORLD, &status);
                tag =2;
		MPI_Recv(n_ptr, 1 , MPI_INT , source , tag , MPI_COMM_WORLD, &status);

	}

}

