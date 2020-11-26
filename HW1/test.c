#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){

    int 	my_rank, p, source, dest,i;          
    int         tag = 0;       
    char        message[100];  
    MPI_Status  status;    /* return status for  receive  */                               
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);


    
   for (i = p-1; i > 0; --i) {
	//MPI_Barrier(MPI_COMM_WORLD);
        if (i != my_rank) {
sleep(1/10);
}else{
            printf("this is proc. %d on a total of %d processes\n", my_rank, p);
      }  
	//}
      //  MPI_Barrier(MPI_COMM_WORLD);
    }
    MPI_Finalize();
}
