#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

double sor(double grid[], int N, double h, GRID_INFO_T pgrid, double w, double aux[], 
	double my_up[],   double my_down[],   double my_left[],   double my_right[],
	double buff_up[], double buff_down[], double buff_left[], double buff_right[]
	)
{
  double max=-1.0;
  double delta=100.0;
  double x,y;
  int ii, jj;
  int i_start=0, j_start=0, i_end=N-1, j_end=N-1;

  memcpy(aux, grid, N*N*sizeof(double));  //for the iteration error
 
//------------Defining the grid start point for p with B.Conditions----
  if(pgrid.my_row==0) {i_start=1;}
  if(pgrid.my_row==pgrid.nrows-1) {i_end=N-2;}
  if(pgrid.my_col==0) {j_start=1;}
  if(pgrid.my_col==pgrid.ncols-1) {j_end=N-2;}

//---------------------SOR step------------------
//
  exchange_boundaries(my_up, buff_up, my_down, buff_down, my_left, buff_left, my_right, buff_right, N,  pgrid);
  int color; //0== REDi/even;  1==BLACK
//if(pgrid.my_rank==0){
  for( ii = i_start; ii <= i_end; ++ii){
    y= 1 -N*pgrid.my_row*h   -ii*h;
    for( jj = j_start; jj <= j_end; ++jj){
      x = N*pgrid.my_col*h +h*jj;
        color=1;
      if( (ii+jj)%2 ==0 ) {color=0; 
        grid[ii*N+jj] = 4*grid[ii*N+jj]* (1-w)/w + h*h*f(x,y) ;
        if(jj+1==N){grid[ii*N + jj] += buff_right[ii];}  else{grid[ii*N+jj]+= grid[ii*N+jj+1];} //Right
        if(jj-1<0){ grid[ii*N + jj] += buff_left[ii]; }  else{grid[ii*N+jj]+= grid[ii*N+jj-1];} //Left
        if(ii+1==N){grid[ii*N + jj] += buff_down[jj]; }  else{grid[ii*N+jj]+= grid[(ii+1)*N+jj];} //Down
        if(ii-1<0){ grid[ii*N + jj] += buff_up[jj];   }  else{grid[ii*N+jj]+= grid[(ii-1)*N+jj];} //UP
        grid[ii*N+jj]*=0.25*w;
        delta=  fabs( aux[ii*N + jj] -  grid[ii*N +jj] ) ;
        if( delta  > max  ){  max=  delta ;}
      } 
     }
   }

//  MPI_Barrier(MPI_COMM_WORLD);  
  exchange_boundaries(my_up, buff_up, my_down, buff_down, my_left, buff_left, my_right, buff_right, N,  pgrid);
  for( ii = i_start; ii <= i_end; ++ii){
    y= 1 -N*pgrid.my_row*h   -ii*h;
    for( jj = j_start; jj <= j_end; ++jj){
      x = N*pgrid.my_col*h +h*jj;
      color=1;
      if( (ii+jj)%2 !=0 ) {color=1; 

        grid[ii*N+jj] = 4*grid[ii*N+jj]* (1-w)/w + h*h*f(x,y) ;
        if(jj+1==N){grid[ii*N + jj] += buff_right[ii];}  else{grid[ii*N+jj]+= grid[ii*N+jj+1];} //Right
        if(jj-1<0){ grid[ii*N + jj] += buff_left[ii]; }  else{grid[ii*N+jj]+= grid[ii*N+jj-1];} //Left
        if(ii+1==N){grid[ii*N + jj] += buff_down[jj]; }  else{grid[ii*N+jj]+= grid[(ii+1)*N+jj];} //Down
        if(ii-1<0){ grid[ii*N + jj] += buff_up[jj];   }  else{grid[ii*N+jj]+= grid[(ii-1)*N+jj];} //UP
        grid[ii*N+jj]*=0.25*w;
        delta=  fabs( aux[ii*N + jj] -  grid[ii*N +jj] ) ;
        if( delta  > max  ){  max=  delta ;}

      } 
     }
   }


//}
  return max;
}

