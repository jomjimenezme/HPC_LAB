#include <stdio.h>
#include "mpi.h"


double Serial_dot(double  x[], double  y[], int  n) {
    int    i; 
    double  sum = 0.0;
    for (i = 0; i < n; i++)
        sum = sum + x[i]*y[i];
    return sum;
}


double Parallel_dot( double  local_x[],  double  local_y[], int  n_bar) {
    double  local_dot;
    double  dot = 0.0;

    local_dot = Serial_dot(local_x, local_y, n_bar);
    MPI_Reduce(&local_dot, &dot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    return dot;
}


