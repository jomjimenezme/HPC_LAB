#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"


double f(double x);

int main(int argc, char** argv){
  //MPI Variables 
  int 	my_rank, p, i, j;
  
  
  //Integral Variables
  int n_trapez, n, counter,a0,b0,n0;
  double x, h,loc_sum, integral;
  double a,b; 
  double local_a, local_b;
  double start, startl, finishl;
  double N=100; //Averaging over N iterations. 
  

  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  //a0=a;b0=b;n0=n;
  for(counter=1; counter<=N; counter++){
  //a=a0;b=b0;n=n0;
  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  Get_data(&a,&b,&n, my_rank, p);  
  
  n=n*p; 
  h=(b-a)/n;
  n_trapez=n/p;
  local_a=a+(my_rank)*n_trapez*h;
  local_b=a+(my_rank+1)*n_trapez*h;
  loc_sum= (f(local_a)+f(local_b))/2.0;
  x=local_a;
  integral=0;
 
  for (i=1; i< n_trapez; i++){
    x=x+h;
    loc_sum+=f(x);
  }
  loc_sum = loc_sum*h;
  //MPI_Barrier(MPI_COMM_WORLD);
  //startl = MPI_Wtime();

  if(my_rank!=0){
    MPI_Send(&loc_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }

  else{
    
    integral=loc_sum;
	for(i=1; i<p; i++){
	  MPI_Recv(&loc_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
	  integral+=loc_sum;
	}
  }
  MPI_Barrier(MPI_COMM_WORLD);
  finishl = MPI_Wtime();
  
  if(my_rank==0){
    //printf("%d %d %e %e %e \n", p, n, finishl-start,finishl-startl, fabs(integral*4-M_PI));
    printf("%d %d %e %e \n", p, n, finishl-start, fabs(integral*4-M_PI));
}
 }
 
  MPI_Finalize();
}


double f( double x){
double eval;
 eval =1/(1+x*x);
 return eval;
}
