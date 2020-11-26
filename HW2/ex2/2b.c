#include<stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){
    
int 	my_rank, p, source, dest,tag, left, right;
char	message[100];          
    
char	ml[100]; 
char    mr[100];          

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

	MPI_Status  status; 
	left=my_rank-1;
	right= my_rank+1;
	if (left<0) left=p-1;
	if (right==p) right=0; 
	
	if(my_rank%2==0){
	sprintf(message, "This is process %d, Hallo!\n", my_rank);
	
	MPI_Send(message,strlen(message)+1,MPI_CHAR, right, 0, MPI_COMM_WORLD);
	MPI_Send(message,strlen(message)+1,MPI_CHAR, left, 0, MPI_COMM_WORLD);
	MPI_Recv(message, 100, MPI_CHAR, right, 1, MPI_COMM_WORLD, &status);
	memcpy(mr, message, sizeof message);	
	MPI_Recv(message, 100, MPI_CHAR, left, 1, MPI_COMM_WORLD, &status);
	memcpy(ml, message, sizeof message);	
	

	printf("I am p=%d, I got this messages: \n %s  %s\n", my_rank, ml,mr);
	}	

	else{

        MPI_Recv(message, 100, MPI_CHAR, left, 0, MPI_COMM_WORLD, &status);
	memcpy(ml, message, sizeof message);	
        MPI_Recv(message, 100, MPI_CHAR, right, 0, MPI_COMM_WORLD, &status);
	memcpy(mr, message, sizeof message);	
	sprintf(message, "This is process %d, Hallo!\n", my_rank);
	MPI_Send(message,strlen(message)+1,MPI_CHAR, left, 1, MPI_COMM_WORLD);
        MPI_Send(message,strlen(message)+1,MPI_CHAR, right, 1, MPI_COMM_WORLD);

	printf("I am p=%d, I got this messages: \n %s  %s\n", my_rank, ml,mr);
	}




    MPI_Finalize();

}
