#include <stdio.h>
#include <string.h>
#include "mpi.h"


void Read_matrix(int my_rank, int p){
	int i;
	int ncols =0, nrows=0;
	char c;
	if(my_rank ==0){
		FILE *input;
	        input = fopen("matrix.d", "r");
//--------------------------------Get Number of Columns---------------------
		for (c=getc(input); c!= EOF && c!='\n'; c= getc(input)){
		    if(c==' ')//count spaces instead of characters
			ncols+=1;
		}ncols+=1;
//--------------------------------Get Number of ROWS---------------------
		rewind(input);
		for (c = getc(input); c != EOF; c = getc(input)){
		  if(c=='\n'){nrows=nrows+1; }
		} 
        	fclose(input);
	printf("%d %d \n", nrows, ncols);
	}


    //MPI_scatter(temp, local_n, MPI_DOUBLE, local_x, local_n, MPI_DOUBLE,o, MPI_COMM_WORLD);
}


void Read_vector(int, my_rank, intp){
	int i;
	


}


