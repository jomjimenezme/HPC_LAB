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
  int n_trapez, n;
  double x, h,loc_sum, integral;
  double aux_sum, TN, T2N, S2N=10;
  double a,b; 
  double local_a, local_b;
  double eps=1.0e-11;
  
  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p); 
  Get_data(&a,&b,&n, my_rank, p);  
  n=2*p; 
  
  local_a=a+(my_rank)*n_trapez*h;
  local_b=a+(my_rank+1)*n_trapez*h;
  loc_sum= (f(local_a)+f(local_b))/2.0;
  aux_sum=loc_sum+f(local_a+h);
  TN_loc=(aux_sum)*h; 
  MPI_Reduce(&loc_sum, &integral, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  TN=integral*h;
  while(fabs(S2N-M_PI)>0.00001){
  loc_sum=0;
  n=n*2; 
  h=(b-a)/n;
  n_trapez=n/p;
  x=local_a;
 
  for (i=1; i< n_trapez; i+=2){
	x=x+h;
	loc_sum+=f(x);
  }
 // loc_sum = loc_sum*h;
  MPI_Reduce(&loc_sum, &integral, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  T2N=(aux_sum+integral)*h; 
  S2N=(4.0/3.0)*T2N-(1.0/3.0)*TN;
  MPI_Broadcast(&S2N,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
 
   if(my_rank==0){	
 printf("%d\t %.16lf\n", n, fabs(S2N*4-M_PI));
 }
 TN=T2N;
 aux_sum=integral;

}
 //}
  MPI_Finalize();
}


double f( double x){
double eval;
eval =1/(1+x*x);
return eval;
}
