#include <stdio.h>
#include <string.h>
#include "mpi.h"


double f(double x);
int main(int argc, char** argv){
 //MPI Variables 
double integral,a, b, h, x;
int n,i;
 	
a=0.0;
b=1.0;
n=32000;

h=(b-a)/n;
integral= (f(a)+f(b))/2.0;
x=a;
for (i=1; i<n; i++){
	x=x+h;
	integral=integral+f(x);
}
integral = integral*h; 
printf("%lf  %lf,  %d, \n, %lf", a,b,n,integral);

 
  MPI_Finalize();
}

double f( double x){
double eval;
eval =1;
//eval =1/(1+x*x);
return eval;
}
//
