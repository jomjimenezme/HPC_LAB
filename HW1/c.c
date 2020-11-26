#include<stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){
    
int 	my_rank, p, source, dest;          
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    printf("Hello World, I am proces %d out of %d processes\n", my_rank, p); //All procesors will print this

    MPI_Finalize();

}
