#include <stdio.h>
#include <string.h>
#include "mpi.h"

int main(int argc, char** argv){
  
  int 	my_rank, p, source, dest, i;          
  int         tag = 0;       
  char        message[100];
  
  char ch;
  int number_names = 0;
  char aux[100]="";
  char nam[100]="";
  
  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  if ( my_rank== 0) {//Every process send to process 0
   dest=1;
    FILE *fpr= fopen("input.d", "r");
    while ((ch = fgetc(fpr)) != EOF)
      {
	strncat(aux, &ch,1);
 	/* Getting number of names */
	if (ch == ' ' || ch == '\t'|| ch == '\0' || ch== EOF){
	 sprintf(nam, "%s",aux);
	 MPI_Send(nam, strlen(nam)+1, MPI_CHAR, dest, 0, MPI_COMM_WORLD);
	 dest++;
	 number_names++;
	 strcpy(nam, "");
	 strcpy(aux, "");
	}
      } 
    fclose(fpr);

  for(dest=1; dest<p; dest++){
      MPI_Send(&number_names, 1, MPI_INT, dest, 2000, MPI_COMM_WORLD);
  }
 
    
  }else{
    source=0;
    dest=0;
      MPI_Recv(&number_names,1, MPI_INT, source,2000, MPI_COMM_WORLD, &status);
      if(my_rank<=number_names){
          MPI_Recv(nam,100, MPI_CHAR, source, 0, MPI_COMM_WORLD, &status);
          sprintf(message, "Hallo %s! I am process ---> %d\n",nam, my_rank);
          MPI_Send(message, strlen(message)+1, MPI_CHAR, dest, 1, MPI_COMM_WORLD);
      }
    }
    
  
  
  if(my_rank==0){
	printf("I got this NAMES: %d", number_names);
    FILE *fp= fopen("output.txt", "w");
    for (source = 1; source <= number_names; source++) {
	MPI_Recv(message,100, MPI_CHAR, source, 1, MPI_COMM_WORLD, &status);
	fputs(message, fp );
    }
    fputs("Hello Everyone, I am the master process (0), I printed everything\n",fp);
    fclose(fp); 
  }
  
  
  MPI_Finalize();
}
