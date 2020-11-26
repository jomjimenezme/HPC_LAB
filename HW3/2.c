#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#include "get_data.c"


double f(double x);

int main(int argc, char** argv){
 //MPI Variables 
  int 	my_rank, p, source, dest, i, j;         
 

//Integral Variables
  int n_trapez, n;
  double x, h,loc_sum, integral;
  double a,b; 
  double local_a, local_b;
  double start, finish;
  
  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p); 
MPI_Barrier(MPI_COMM_WORLD);
start = MPI_Wtime();
  Get_data(&a,&b,&n, my_rank, p);  
  
  
  h=(b-a)/n;
 // printf("localsum=  %d,    x=%d,    h=%d, n_trapez= %d\n",loc_sum,x, n_trapez);

  n_trapez=n/p;
  local_a=a+(my_rank)*n_trapez*h;
  local_b=a+(my_rank+1)*n_trapez*h;
  loc_sum= (f(local_a)+f(local_b))/2.0;
  x=local_a;
 
  for (i=1; i< n_trapez; i++){
	x=x+h;
	loc_sum+=f(x);
  }
  loc_sum = loc_sum*h;

 if(my_rank!=0){
	MPI_Send(&loc_sum, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
 }else{
	integral=loc_sum;
	for(i=1; i<p; i++){
		MPI_Recv(&loc_sum, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
		integral+=loc_sum;
	}
	
 printf("%d\t %.16lf\n", n, fabs(integral*4-M_PI));
 //printf("For n=%d, and p=%d processors, we found \n pi = %.16lf\n", n, p, integral*4-M_PI);
 }
MPI_Barrier(MPI_COMM_WORLD);
finish = MPI_Wtime();
if(my_rank==0)
printf("Elapsed time = %e seconds \n", finish-start);

  MPI_Finalize();
}


double f( double x){
double eval;
eval =1/(1+x*x);
return eval;
}
