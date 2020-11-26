#include<stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){

    MPI_Init(&argc, &argv);
    
    printf("Hello World\n"); //All procesors will print this

    MPI_Finalize();

}
