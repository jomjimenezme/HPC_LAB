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
  double x, h,loc_sum,integral;
  double aux_sum, sum, TN, T2N, S2N=10;
  double a,b; 
  double local_a, local_b;
  double eps=1.0e-15;
  
  MPI_Status  status;    /* return status for  receive  */                               
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p); 
  Get_data(&a,&b,&n, my_rank, p);  
  n=2*p; 
  
  h=(b-a)/n;
  n_trapez=n/p;
  local_a=a+(my_rank)*n_trapez*h;
  local_b=a+(my_rank+1)*n_trapez*h;
  loc_sum= (f(local_a)+f(local_b))/2.0;
  loc_sum=loc_sum+f(local_a+h);
  MPI_Reduce(&loc_sum, &sum, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&sum, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  integral=sum*h;
  TN=integral;
  S2N=TN;
  while(fabs(4*S2N-M_PI)>eps){
  loc_sum=0;
  n=n*2; 
  h=(b-a)/n;
  n_trapez=n/p;
  x=local_a;
 
  for (i=1; i< n_trapez; i+=2){
	x=local_a+i*h;
	loc_sum+=f(x);
  }
  MPI_Reduce(&loc_sum, &aux_sum, 1,  MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&aux_sum,1,MPI_DOUBLE, 0, MPI_COMM_WORLD);
  T2N=(aux_sum+sum)*h; 
  S2N=(4.0/3.0)*T2N-(1.0/3.0)*TN;
  sum= sum+aux_sum;
   if(my_rank==0){	
 printf("%d\t %.20e\n", n, fabs(S2N*4-M_PI));
 }
 TN=T2N;

}
 //}
  MPI_Finalize();
}


double f( double x){
double eval;
eval =1/(1+x*x);
return eval;
}
