#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){

    int 	my_rank, p, source, dest;          
    int         tag = 0;       
    char        message[100];  
    MPI_Status  status;    /* return status for  receive  */                               
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if ( my_rank!= 0) {
        /*  Write into a character string buffer "MESSAGE*/
 	sprintf(message, "This is process %d, out of %d processes", my_rank,p);        
	dest = 0;
        /* Use strlen+1 so that '\0' gets transmitted */
        MPI_Send(message, strlen(message)+1, MPI_CHAR, 
            dest, tag, MPI_COMM_WORLD);
    } else { /* my_rank == 0 */
        for (source = p-1; source > 0; source--) {
            MPI_Recv(message,//buffer
		 100,//count of ellements in buffer
		 MPI_CHAR, //Type of elements in buffer
 		 source, //rank of receiver
		 tag, //Tag of message, 
                MPI_COMM_WORLD, &status); //communicator;
            printf("%s\n", message);
        }
    }

    MPI_Finalize();
}
