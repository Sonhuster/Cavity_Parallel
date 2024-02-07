#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include "initialize.h"
#include "Decompose.h"
#include "Update.h"


#define NROWS 		10
#define NCOLS		10
#define Min_Procs 	4
#define MASTER 		0
#define trans_u 	1
#define trans_v 	2
#define trans_p 	3
#define trans_nei 	4
#define NONE 		0

int hx =1, hy =1;
float delta = 4.5, dt = 0.001, Re = 100;


//void initialize_const(int nx, int ny, float *uu, float value);
//int isPerfectSquare(int number);
//int Cal_dim_length(int num_task);
//void print_array(float *uu, int nx, int ny);

int main (int argc, char *argv[])
{	
	
	int num_tasks, my_rank, dim_mid[2];
	
	int *coords = Decompose(num_tasks, my_rank, argc, argv, &dim_mid[0]);

	//printf("\n position: %d ,%d \n", coords[0], coords[1]);
	
	MPI_Finalize();
	
}

