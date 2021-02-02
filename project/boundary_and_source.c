#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


double f(double x, double y){
  return 2*( (1+x) * sin(x+y)*cos(x+y) );
}

double g(double x, double y){
  return (1+x)*sin(x+y);
}


void boundary_conditions(double grid[], int N, double h, GRID_INFO_T pgrid)
{
  int ii, jj;
  double x, y;
  if (pgrid.my_row==0){
    ii = 0; //upper border
    y=1.0;//(N-1)*h; 
    for(jj = 0; jj < N; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if (pgrid.my_row==pgrid.nrows-1){
    ii = N-1; //lower border
    y = 0.0;
    for(jj = 0; jj < N; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if (pgrid.my_col==0){
    jj = 0;  //left border
    x=0.0;
    for(ii = 0; ii < N; ++ii){
      y= 1 -N*pgrid.my_row*h   -ii*h;
      grid[ii*N + jj] = g(x,y);
    }
  }
  if(pgrid.my_col==pgrid.ncols-1){
    jj = N-1;//right border
    x = 1.0;
    for(ii = 0; ii < N; ++ii){
      y= 1 -N*pgrid.my_row*h   -ii*h;
      grid[ii*N + jj] = g(x,y);
    }
  }
}


