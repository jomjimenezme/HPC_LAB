#include<stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"
#include "algebra.c"
int main(int argc, char** argv){
    
  int my_rank, p, i,j,k, c, count;    
  int local_nA, local_mA, mA,nA;
  int local_nB, local_mB, mB,nB;
  int local_nC, local_mC, mC,nC;
  double start, finish, startc;
  int shift=0;
  int N=100;
  double* local_C;
  double* local_A;
  double* local_B; 
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Status  status; 

  for (count=0; count<N; count++){ //So that we measure the time several times

//-----------------Allocation of local matrixes------------------
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  Read_size(&mA, &nA, my_rank, p, "Ascaling.txt");
  Read_size(&mB, &nB, my_rank, p, "Bscaling.txt");
  local_mA=mA/p;
  local_nA=nA;
  local_mB=mB/p;
  local_nB=nB;
  
  mC=mA;
  nC=nB;
  local_mC= local_mA;
  local_nC= local_nB;

  local_A= malloc( local_mA*local_nA *sizeof( double ) );
  local_B= malloc( local_mB*local_nB *sizeof( double ) );
  local_C= malloc( local_mC*local_nC *sizeof( double ) );

  Read_matrix(local_A, mA, nA, my_rank, p, "A.txt");
  Read_matrix(local_B, mB, nB, my_rank, p, "B.txt");
  
  memset( local_C, 0, local_mC * local_nC * sizeof(double) ); 
//------------------Computation of C----------------------

  MPI_Barrier(MPI_COMM_WORLD);
  startc = MPI_Wtime();

  for(c=0; c<=local_mB; c++){ //Loop for circular shift of B.
    shift= (c +  my_rank)%p  ; //start position for A after every stage of rotation
    for (i=0; i<local_mC; i++){ // The three loops of Matrix-Matrix product
      for (j=0; j<local_nC; j++){
	for(k=0; k< local_mB; k++){
          local_C[  i*local_nC +j ] += local_A[ ( shift )*local_mB +  i*local_nA+k  ]  *  local_B[ k * local_nB +j  ];
	}
      }
    }
    MPI_Sendrecv_replace(local_B, nB*local_mB, MPI_DOUBLE,   ( (my_rank-1)%p +p )%p  , 0, (my_rank+1)%p , 0, MPI_COMM_WORLD, &status );
  }
 


  MPI_Barrier(MPI_COMM_WORLD);
  finish = MPI_Wtime();

  if(my_rank==0) printf("%d %d %d", p, finish-start, finish-startc);
 
  free(local_A);
  free(local_B);
  free(local_C);

}
  MPI_Finalize();
      
}


/*
 *
//------------------Silly printing------------------------------
 if(my_rank==0) printf("%d\n%d\n", mC, nC);  
MPI_Barrier(MPI_COMM_WORLD);
 for (c=0; c<p; c++){
    if(my_rank==c){
      for(i=0; i<local_mC; i++){
        for(j=0; j<local_nC; j++){
  	  printf("%lf\n", local_C[i*local_nC+j]);
        }
  	  //printf("\n");
      }
    }
  }

*/
