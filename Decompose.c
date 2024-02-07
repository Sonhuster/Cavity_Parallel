#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "Decompose.h" 
#include "initialize.h"
#include "Update.h"


#define NROWS 		11		//Grid Size X
#define NCOLS		11		//Grid Size Y
#define Iteration   2000		//Define Iteration
#define Tolerance   0.0001	//Define required error
#define trans_u 	1
#define trans_v 	2
#define trans_p 	3
#define trans_nei_u 4
#define trans_nei_v 5
#define trans_nei_p 6
#define Min_Procs 	4
#define MASTER 		0
#define hx 			0.01
#define hy 			0.01



void Check();
int isPerfectSquare(int number);
int Cal_dim_length(int num_tasks);
int *Sender (int dimx, int dimy, int TAG, int *coords, int my_rank, MPI_Comm comm);
void print_array(float *uu, int nx, int ny);
void initialize_test(int nx, int ny, float *uu, float value);
void Send_neighbor (float *uu, int nx, int ny, int left, int right, int top, int bottom, int *coords, int dimx , int dimy, int TAG, MPI_Comm comm);
void Recv_neighbor (float *uu, int nx, int ny, int left, int right, int top, int bottom, int *coords, int dimx , int dimy, int TAG, MPI_Comm comm);


int *Decompose(int num_tasks, int my_rank, int argc, char *argv[], int *dim_mid)
{
	int nx = NROWS, ny = NCOLS;
	

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
	MPI_Status status;
	// Create a 2D cartesian grid (2x2)
	int dimx, dimy;
	dimx = Cal_dim_length(num_tasks); dimy = num_tasks/dimx; 
	//printf ("Full array will be splited by %d blocks in X direction, %d blocks in Y direction\n",dimx, dimy);
	
	*dim_mid = dimx; *(dim_mid + 1) = dimy;
	int dim[2] = {dimx, dimy};
	int periods[2] = {0, 0};
	int reorder = 0;

	MPI_Comm comm;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim, periods, reorder, &comm);
	MPI_Comm_rank(comm, &my_rank);
	// Get the ranks of the left, right, bottom, and top neighbors
    int left, right, bottom, top;
    MPI_Cart_shift(comm, 0, 1, &top, &bottom);  
    MPI_Cart_shift(comm, 1, 1, &left, &right);
	
    if ((dim[0] > sqrt(NROWS)) ||(dim[1] > sqrt(NCOLS)) || (num_tasks < Min_Procs)) {
        printf("ERROR: the number of tasks must be between %d and %d.\n", Min_Procs, (int)(sqrt(NROWS)*sqrt(NCOLS)) );
        printf("Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
	
	
	// Get the coordinates of the current process in the cartesian grid
	int *coords = (int *)malloc(2*sizeof(int*));
	MPI_Cart_coords(comm, my_rank, 2, coords);
	int centerx = coords[0],
		centery = coords[1];
	
	/*-----------------------------------------------------*/
	
	if (my_rank < num_tasks)
	{

		//Send initialize
		int *usize = Sender (dimx, dimy, trans_u, coords, my_rank, comm);
		int *vsize = Sender (dimx, dimy, trans_v, coords, my_rank, comm);
		int *psize = Sender (dimx, dimy, trans_p, coords, my_rank, comm);
		//printf("Rank: %d, coordinates: %d, %d. Size(%d,%d)\n",my_rank, coords[0], coords[1],usize[0],usize[1]); 
		/*Recv procedure*/
		float sub_u[usize[0]][usize[1]], sub_v[vsize[0]][vsize[1]], sub_p[psize[0]][psize[1]];
		
		int source = my_rank; 
		MPI_Recv(&sub_u[0][0], usize[0]*usize[1], MPI_FLOAT, source, trans_u, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&sub_v[0][0], vsize[0]*vsize[1], MPI_FLOAT, source, trans_v, comm, MPI_STATUS_IGNORE);
		MPI_Recv(&sub_p[0][0], psize[0]*psize[1], MPI_FLOAT, source, trans_p, comm, MPI_STATUS_IGNORE);
		
		//Iteration Calculation
		for (int it = 0; it <= Iteration; it ++)
		{
		//Update N-s Equations
		update(&sub_u[0][0], &sub_v[0][0], &sub_p[0][0], psize[0], psize[1], top, bottom, left, right);
		//Send to Neighborhood
		Send_neighbor (&sub_u[0][0], usize[0], usize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_u, comm);
		Send_neighbor (&sub_v[0][0], vsize[0], vsize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_v, comm);
		Send_neighbor (&sub_p[0][0], psize[0], psize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_p, comm);
	
		//Recv from Neighborhood
		Recv_neighbor (&sub_u[0][0], usize[0], usize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_u, comm);
		Recv_neighbor (&sub_v[0][0], vsize[0], vsize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_v, comm);
		Recv_neighbor (&sub_p[0][0], psize[0], psize[1], left, right, top, bottom, coords, dimx, dimy, trans_nei_p, comm);
		//Sendback procudure
	
		int xstart, ystart, ux_end, uy_end, ux_size, uy_size,
							vx_end, vy_end, vx_size, vy_size,
							px_end, py_end, px_size, py_size;
		if(centerx==0)
		{
			xstart = 0; ux_end = usize[0] - 2; ux_size = ux_end - xstart + 1;
						vx_end = vsize[0] - 2; vx_size = vx_end - xstart + 1;
						px_end = psize[0] - 2; px_size = px_end - xstart + 1;
		}
		if(centery==0)
		{
			ystart = 0; uy_end = usize[1] - 2; uy_size = uy_end - ystart + 1;
						vy_end = vsize[1] - 2; vy_size = vy_end - ystart + 1;
						py_end = psize[1] - 2; py_size = py_end - ystart + 1;
		}
		if(centerx>0 && centerx<dim[0]-1)
		{
			xstart = 1; ux_end = usize[0] - 2; ux_size = ux_end - xstart + 1;
						vx_end = vsize[0] - 2; vx_size = vx_end - xstart + 1;
						px_end = psize[0] - 2; px_size = px_end - xstart + 1;
		}
		if(centery>0 && centery<dim[0]-1)
		{
			ystart = 1; uy_end = usize[1] - 2; uy_size = uy_end - ystart + 1;
						vy_end = vsize[1] - 2; vy_size = vy_end - ystart + 1;
						py_end = psize[1] - 2; py_size = py_end - ystart + 1;
		}
		if(centerx==dim[0]-1)
		{
			xstart = 2; ux_end = usize[0] - 1; ux_size = ux_end - xstart + 1;
						vx_end = vsize[0] - 1; vx_size = vx_end - xstart + 1;
						px_end = psize[0] - 1; px_size = px_end - xstart + 1;
		}
		if(centery==dim[1]-1)
		{
			ystart = 2; uy_end = usize[1] - 1; uy_size = uy_end - ystart + 1;
						vy_end = vsize[1] - 1; vy_size = vy_end - ystart + 1;
						py_end = psize[1] - 1; py_size = py_end - ystart + 1;
		}
		float u_sendback[ux_size][uy_size], v_sendback[vx_size][vy_size], p_sendback[px_size][py_size];
		int m = 0, n = 0;
	


		for (int i = xstart; i <= ux_end; i++)
		{
			for (int j = ystart; j <= uy_end; j++)
			{
				u_sendback[m][n] = sub_u[i][j];
				n++;
			}
			m++;
			n = 0;
		}

	
		m = 0;
		for (int i = xstart; i <= vx_end; i++)
		{
			for (int j = ystart; j <= vy_end; j++)
			{
				v_sendback[m][n] = sub_v[i][j];
				n++;
			}
			m++;
			n = 0;
		}

	
		m = 0;
		for (int i = xstart; i <= px_end; i++)
		{
			for (int j = ystart; j <= py_end; j++)
			{
				p_sendback[m][n] = sub_p[i][j];
				n++;
			}
			m++;
			n = 0;
		}
		if (my_rank != MASTER)
		{
		MPI_Send(&u_sendback[0][0], ux_size*uy_size, MPI_FLOAT, MASTER, trans_u, comm);
		MPI_Send(&v_sendback[0][0], vx_size*vy_size, MPI_FLOAT, MASTER, trans_v, comm);
		MPI_Send(&p_sendback[0][0], px_size*py_size, MPI_FLOAT, MASTER, trans_p, comm);
		}
		
		if(my_rank == MASTER)
		{
			int 	local_nrows = (NROWS%dim[0] == 0) ?NROWS/dim[0] : NROWS/dim[0] + 1; //const
			int 	local_ncols = (NCOLS%dim[1] == 0) ?NCOLS/dim[1] : NCOLS/dim[1] + 1; //const
			float 	Collected_P[(dim[0]-1)*(dim[1]-1)][local_nrows][local_ncols],
					Collected_U[(dim[0]-1)*(dim[1]-1)][local_nrows][local_ncols],
					Collected_V[(dim[0]-1)*(dim[1]-1)][local_nrows][local_ncols];
			float 	Full_U [NROWS-1][NCOLS],
					Full_V [NROWS][NCOLS-1],
					Full_P [NROWS][NCOLS];
			
			//Rank = MASTER, Impose the array
			for (int i = 0; i < local_nrows; i++)
			{
				for (int j = 0; j < local_ncols; j++)
				{
					Full_U[i][j] = u_sendback[i][j];
					Full_V[i][j] = v_sendback[i][j];
					Full_P[i][j] = p_sendback[i][j];
				}
			}
			
			//Internal_array - Except rank MASTER
			int 	collected_source = 1, count = 1;
			
			for (int ix = 0; ix < dim[0] - 1; ix ++)
			{
				for (int iy = 0; iy < dim[1] - 1; iy ++)
				{
					if ((ix == 0) && (iy == 0)) {continue;}
					MPI_Recv(&Collected_U[count][0][0], local_nrows*local_ncols, MPI_FLOAT, collected_source, trans_u, comm, MPI_STATUS_IGNORE);
					MPI_Recv(&Collected_V[count][0][0], local_nrows*local_ncols, MPI_FLOAT, collected_source, trans_v, comm, MPI_STATUS_IGNORE);
					MPI_Recv(&Collected_P[count][0][0], local_nrows*local_ncols, MPI_FLOAT, collected_source, trans_p, comm, MPI_STATUS_IGNORE);
					//int m = 0; n =0;
					for (int i = 0; i < local_nrows; i++)
					{
						for (int j = 0; j < local_ncols; j++)
						{
							Full_U[ix*local_nrows + i][iy*local_ncols + j] = Collected_U[count][i][j];
							Full_V[ix*local_nrows + i][iy*local_ncols + j] = Collected_V[count][i][j];
							Full_P[ix*local_nrows + i][iy*local_ncols + j] = Collected_P[count][i][j];
						}
					}
					count ++;
					collected_source ++;
				}
				collected_source ++;
			}
			
			//right_array
			int 	local_ncols_extra = NCOLS - (dim[1]-1)*local_ncols; 
					collected_source = dim[0]-1, count = 0;
			float 	Collected_P_col[dim[0]-1][local_nrows][local_ncols_extra],
					Collected_U_col[dim[0]-1][local_nrows][local_ncols_extra],
					Collected_V_col[dim[0]-1][local_nrows][local_ncols_extra-1];
			for (int ix = 0; ix < dim[0] - 1; ix ++)
			{
				MPI_Recv(&Collected_U_col[count][0][0], local_nrows*local_ncols_extra, MPI_FLOAT, collected_source, trans_u, comm, MPI_STATUS_IGNORE);
				MPI_Recv(&Collected_V_col[count][0][0], local_nrows*(local_ncols_extra-1), MPI_FLOAT, collected_source, trans_v, comm, MPI_STATUS_IGNORE);
				MPI_Recv(&Collected_P_col[count][0][0], local_nrows*local_ncols_extra, MPI_FLOAT, collected_source, trans_p, comm, MPI_STATUS_IGNORE);
				for (int i = 0; i < local_nrows; i++)
				{
					for (int j = 0; j < local_ncols_extra; j++)
					{
						Full_U[ix*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_U_col[count][i][j];
						Full_P[ix*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_P_col[count][i][j];
						if (j < local_ncols_extra - 1)
						{Full_V[ix*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_V_col[count][i][j];}
					}
				}
				collected_source = collected_source + dim[1];
				count ++;
			}
			
			//bottom_array
			int 	local_nrows_extra = NROWS - (dim[0]-1)*local_nrows; 
					collected_source = (dim[0]-1)*dim[1], count = 0;
			float 	Collected_P_row[dim[1]-1][local_nrows_extra][local_ncols],
					Collected_U_row[dim[1]-1][local_nrows_extra-1][local_ncols],
					Collected_V_row[dim[1]-1][local_nrows_extra][local_ncols];
			for (int iy = 0; iy < dim[1] - 1; iy ++)
			{
				MPI_Recv(&Collected_U_row[count][0][0], (local_nrows_extra-1)*local_ncols, MPI_FLOAT, collected_source, trans_u, comm, MPI_STATUS_IGNORE);
				MPI_Recv(&Collected_V_row[count][0][0], local_nrows_extra*local_ncols, MPI_FLOAT, collected_source, trans_v, comm, MPI_STATUS_IGNORE);
				MPI_Recv(&Collected_P_row[count][0][0], local_nrows_extra*local_ncols, MPI_FLOAT, collected_source, trans_p, comm, MPI_STATUS_IGNORE);
				for (int i = 0; i < local_nrows_extra; i++)
				{
					for (int j = 0; j < local_ncols; j++)
					{
						Full_V[(dim[0]-1)*local_nrows + i][iy*local_ncols + j] = Collected_V_row[count][i][j];
						Full_P[(dim[0]-1)*local_nrows + i][iy*local_ncols + j] = Collected_P_row[count][i][j];
						if (i < local_nrows_extra - 1)
						{Full_U[(dim[0]-1)*local_nrows + i][iy*local_ncols + j] = Collected_U_row[count][i][j];}
					}
				}
				collected_source ++;
				count ++;
			}
			
			//corner_array
					collected_source = dim[0]*dim[1] - 1;
			float 	Collected_P_corner[local_nrows_extra][local_ncols_extra],
					Collected_U_corner[local_nrows_extra-1][local_ncols_extra],
					Collected_V_corner[local_nrows_extra][local_ncols_extra-1];
			
			MPI_Recv(&Collected_U_corner[0][0], (local_nrows_extra-1)*local_ncols_extra, 	MPI_FLOAT, collected_source, trans_u, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&Collected_V_corner[0][0], local_nrows_extra*(local_ncols_extra-1), 	MPI_FLOAT, collected_source, trans_v, comm, MPI_STATUS_IGNORE);
			MPI_Recv(&Collected_P_corner[0][0], local_nrows_extra*local_ncols_extra, 		MPI_FLOAT, collected_source, trans_p, comm, MPI_STATUS_IGNORE);
			
			for (int i = 0; i < local_nrows_extra; i++)
			{
				for (int j = 0; j < local_ncols_extra; j++)
				{
					Full_P[(dim[0]-1)*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_P_corner[i][j];
					if (i < local_nrows_extra - 1)
					{Full_U[(dim[0]-1)*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_U_corner[i][j];}
					if (j < local_ncols_extra - 1)
					{Full_V[(dim[0]-1)*local_nrows + i][(dim[1]-1)*local_ncols + j] = Collected_V_corner[i][j];}
				}
			}
			
			//Continuity residual as error measure
			float error = 0;
			for (int i = 1; i <= NROWS-2; i++)
			{
				for (int j = 1; j <= NCOLS-2; j++)
				{
					error = error + abs(Full_U[i][j] - Full_U[i][j-1])/1 + abs(Full_V[i][j] - Full_V[i-1][j])/1;
				}
			}
			//After the converged solution, we map the staggered variables to collocated variables
			float u_final[NROWS-1][NCOLS-1], v_final[NROWS-1][NCOLS-1], p_final[NROWS-1][NCOLS-1];
			for (int i = 0; i < NROWS-1; i++)
			{
				for (int j = 0; j < NCOLS-1; j++)
				{
					u_final[i][j] = 0.5*(Full_U[i][j]+Full_U[i][j+1]);
					v_final[i][j] = 0.5*(Full_V[i][j+1]+Full_V[i][j]);
					p_final[i][j] = 0.25*(Full_P[i][j] + Full_P[i+1][j] + Full_P[i][j+1] + Full_P[i+1][j+1]);
				}
			}
			if (it  % 50 == 0) {printf("\nIteration: %d,	Error: %f", it, error);}
			if (it  == Iteration) {print_array(&u_final[0][0], NROWS-1, NCOLS-1);}
			
			// OUTPUT DATA
			FILE *fout2, *fout3;
			fout2 = fopen("UVP.plt","w+t");
			fout3 = fopen("Central_U.plt","w+t");

			if ( fout2 == NULL )
			{
				printf("\nERROR when opening file\n");
				fclose( fout2 );
			}

			else
			{
				fprintf( fout2, "VARIABLES=\"X\",\"Y\",\"U\",\"V\",\"P\"\n");
				fprintf( fout2, "ZONE  F=POINT\n");
				fprintf( fout2, "I=%d, J=%d\n", NROWS-1, NCOLS-1 );

				for ( int j = 0 ; j < NCOLS-1 ; j++ )
				{
					for ( int i = 0 ; i < NROWS-1 ; i++ )
					{
						float xpos, ypos;
						xpos = i*hx;
						ypos = j*hy;
	
						fprintf( fout2, "%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\t%5.8lf\n", xpos, ypos, u_final[i][j], v_final[i][j], p_final[i][j] );
					}
				}
			}
	
			fclose( fout2 );
	
			// CENTRAL --U
			fprintf(fout3, "VARIABLES=\"U\",\"Y\"\n");
			fprintf(fout3, "ZONE F=POINT\n");
			fprintf(fout3, "I=%d\n", NROWS-1 );
			for ( int j = 0 ; j < NCOLS-1 ; ++j )
			{
				float ypos = j*hy;
				if ((NROWS-1) % 2 != 0)fprintf( fout3, "%0.6f \t %0.6f\n", (u_final[(NROWS-1)/2][j]), ypos );
				if ((NROWS-1) % 2 == 0)fprintf( fout3, "%0.6f \t %0.6f\n", (u_final[NROWS/2][j] + u_final[(NROWS/2)-1][j])/2, ypos);
			}
			fclose( fout3 );
		}
		}
	
	MPI_Barrier(comm);
	}
	
	
	return coords;
	
}

int *Sender (int dimx, int dimy, int TAG, int *coords, int my_rank, MPI_Comm comm)
{
	float 	u[NROWS-1][NCOLS], 
			v[NROWS][NCOLS-1], 
			p[NROWS][NCOLS];
	initialize(NROWS, NCOLS, &u[0][0], &v[0][0], &p[0][0]);
	//initialize_test(NROWS-1, NCOLS, &u[0][0], 0); 
	//initialize_test(NROWS, NCOLS-1, &v[0][0], 0); 
	//initialize_test(NROWS, NCOLS, &p[0][0], 0); 
	
	int centerx, centery;
	int dim[2]; 
	dim[0] = dimx; dim[1] = dimy;
	//printf("\nRank:%d. Position:  %d, %d \n",my_rank, coords[0], coords[1]);
	centerx = coords[0];
	centery = coords[1];
	
	
	int xstart, xend, xsize,
		ystart, yend, ysize,
		local_nrows, local_ncols;
	
		local_nrows = (NROWS%dim[0] == 0) ?NROWS/dim[0] : NROWS/dim[0] + 1; //const
		local_ncols = (NCOLS%dim[1] == 0) ?NCOLS/dim[1] : NCOLS/dim[1] + 1; //const
			
			
		if(centerx==dim[0]-1)
	  	{
			if(TAG == trans_u) 
			{local_nrows = NROWS - 1 - (dim[1]-1)*local_nrows;}
			if((TAG == trans_v) || (TAG == trans_p)) 			
			{local_nrows = NROWS - (dim[1]-1)*local_nrows;}
	  	}
		if(centery==dim[1]-1)
	  	{
			if(TAG == trans_v) 
			{local_ncols = NCOLS - 1 - (dim[0]-1)*local_ncols;}
			if((TAG == trans_u) || (TAG == trans_p))			
			{local_ncols = NCOLS - (dim[0]-1)*local_ncols;}
	  	}
			
		if(centerx==0)
		{
			xstart = centerx*local_nrows;
			xend =  (centerx+1)*local_nrows;
			xsize = xend - xstart +1;
		}
		if(centery==0)
		{ 	
			ystart = centery*local_ncols;
			yend = (centery+1)*local_ncols;
			ysize = yend - ystart +1;
		}
		if(centerx>0 && centerx<dim[0]-1)
		{
			xstart = centerx*local_nrows-1;
			xend =  (centerx+1)*local_nrows;
			xsize = xend - xstart +1;
		}
		if(centery>0 && centery<dim[1]-1)
		{ 
			ystart = centery*local_ncols-1;
			yend =  (centery+1)*local_ncols;
			ysize = yend - ystart +1;
		}
		if(centerx==dim[0]-1)
		{
			if(TAG == trans_u){
			xend = NROWS - 2;}
			if((TAG == trans_v) || (TAG == trans_p)){
			xend = NROWS - 1;}
			xstart = xend - local_nrows - 1;
			xsize = xend - xstart +1;
		}
		if(centery==dim[1]-1)
		{ 
			if(TAG == trans_v){
			yend =  NCOLS - 2;}
			if((TAG == trans_u) || (TAG == trans_p)){
			yend =  NCOLS - 1;}
			ystart = yend - local_ncols - 1;
			ysize = yend - ystart + 1; 
		}
			
		int array_P[2] = {NROWS,NCOLS}, array_U[2] = {NROWS-1,NCOLS}, array_V[2] = {NROWS,NCOLS-1};
		int subsizes[2] = {xsize,ysize};
		int ar_start[2]= {0,0};
			
		MPI_Datatype ar_subsizes;
		if(TAG == trans_u) 	{MPI_Type_create_subarray(2, array_U, subsizes, ar_start, MPI_ORDER_C, MPI_FLOAT, &ar_subsizes);}
		if(TAG == trans_v) 	{MPI_Type_create_subarray(2, array_V, subsizes, ar_start, MPI_ORDER_C, MPI_FLOAT, &ar_subsizes);}
		else				{MPI_Type_create_subarray(2, array_P, subsizes, ar_start, MPI_ORDER_C, MPI_FLOAT, &ar_subsizes);}
			
		MPI_Type_commit (&ar_subsizes);
		/*Send procedure*/
		int dest = my_rank;
		if(TAG == trans_u)	{MPI_Send(&u[xstart][ystart], 1, ar_subsizes, dest, trans_u, comm);}
		if(TAG == trans_v) 	{MPI_Send(&v[xstart][ystart], 1, ar_subsizes, dest, trans_v, comm);}
		if(TAG == trans_p)	{MPI_Send(&p[xstart][ystart], 1, ar_subsizes, dest, trans_p, comm);}
			
		MPI_Type_free (&ar_subsizes);
			
		int *size = (int *)malloc(2*sizeof(int*));
		size[0] = xsize; size[1] = ysize;
		return size;
		
}

void Send_neighbor (float *uu, int nx, int ny, int left, int right, int top, int bottom, int *coords, int dimx , int dimy, int TAG, MPI_Comm comm)
{
	float u[nx][ny];
	int nei_xstart = 1, nei_ystart = 1;  
	int row_count = ny-2, col_count = nx-2;
	if(top < 0) 	{col_count = nx-1;} 
	if(left < 0) 	{row_count = ny-1;}
	if(top < 0) 	{nei_xstart = 0;} 		if(bottom < 0) {nei_xstart = 2;} 
	if(left < 0) 	{nei_ystart = 0;} 		if(right < 0) {nei_ystart = 2;}
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            u[i][j] = 0;
			u[i][j] = *(uu + i*ny + j);
        }
    }
	if (left >= 0)
	{
		MPI_Datatype left_col;
		MPI_Type_vector (col_count, 1, ny, MPI_FLOAT, &left_col);			//sửa nếu coords[0] == 0 thì nx-2 thành nx -1
		MPI_Type_commit (&left_col);
		if (right < 0){
		MPI_Send(&u[nei_xstart][2], 1, left_col, left, TAG, comm);}
		if (right >= 0){
		MPI_Send(&u[nei_xstart][1], 1, left_col, left, TAG, comm);}
		MPI_Type_free(&left_col);
	}
	if (right >= 0)
	{	
		MPI_Datatype right_col;
		if(coords[1] != dimy - 2)
		{
			MPI_Type_vector (col_count, 1, ny, MPI_FLOAT, &right_col);
			MPI_Type_commit (&right_col);
			MPI_Send(&u[nei_xstart][ny-2], 1, right_col, right, TAG, comm);
			MPI_Type_free(&right_col);
		}
		if(coords[1] == dimy - 2)
		{
			MPI_Type_vector (col_count, 2, ny, MPI_FLOAT, &right_col);
			MPI_Type_commit (&right_col);
			MPI_Send(&u[nei_xstart][ny-3], 1, right_col, right, TAG, comm);
			MPI_Type_free(&right_col);
		}
	}
	if (top >= 0)
	{
		if (bottom < 0){
		MPI_Send(&u[2][nei_ystart], row_count, MPI_FLOAT, top, TAG, comm);}
		if (bottom >= 0){
		MPI_Send(&u[1][nei_ystart], row_count, MPI_FLOAT, top, TAG, comm);} //sửa nếu coords[1] == 0 thì ny-2 thành ny -1
	}
	if (bottom >= 0)
		if(coords[0] != dimx - 2){
		MPI_Send(&u[nx-2][nei_ystart], row_count, MPI_FLOAT, bottom, TAG, comm);}
		if(coords[0] == dimx - 2)
		{
			MPI_Datatype bottom_row;
			MPI_Type_vector (2, row_count, ny, MPI_FLOAT, &bottom_row); 
			MPI_Type_commit (&bottom_row);
			MPI_Send(&u[nx-3][nei_ystart], 1, bottom_row, bottom, TAG, comm);
			MPI_Type_free(&bottom_row);
		}
}

void Recv_neighbor (float *uu, int nx, int ny, int left, int right, int top, int bottom, int *coords, int dimx , int dimy, int TAG, MPI_Comm comm) // TAG: u,v,p
{
	float rows[ny-2], rows_dim[2*(ny-2)], cols[nx-2], cols_dim[2*(nx-2)];
	int n;
	int nei_xstart = 1, nei_ystart = 1; 
	int row_count = ny-2, col_count = nx-2;
	if(top < 0) {col_count = nx-1;} 	if(left < 0) {row_count = ny-1;}	
	if(top < 0) 	{nei_xstart = 0;} 		if(left < 0) {nei_ystart = 0;}
	if(bottom < 0) 	{nei_xstart = 2;} 		if(right < 0) {nei_ystart = 2;}
	
	//Recieve procedure - left,  suy ra nhận từ source là left vì left của left của my_rank là left của sender;
	if (left >= 0)
	{
		if (right >= 0)
		{
			MPI_Recv(&cols[0], col_count, MPI_FLOAT, left, TAG, comm, MPI_STATUS_IGNORE);
			
			n = 0;
			for (int i = nei_xstart; i < col_count + nei_xstart; i++)
			{
				*(uu + i*ny)= cols[n]; n++;
			}
		}
		if (right < 0)
		{
			MPI_Recv(&cols_dim[0], 2*col_count, MPI_FLOAT, left, TAG, comm, MPI_STATUS_IGNORE);
			n = 0;			
			for (int i = nei_xstart; i < col_count + nei_xstart; i++)
			{
				for(int j = 0; j < 2; j++)
				{
					*(uu + i*ny + j) = cols_dim[n]; 
					n++;
				}
			}
		}
		
	}
	if (right >= 0)
	{	
		MPI_Recv(&cols[0], col_count, MPI_FLOAT, right, TAG, comm, MPI_STATUS_IGNORE);
		n = 0;
		for (int i = nei_xstart; i < col_count + nei_xstart; i++)
		{
			*(uu + i*ny + ny-1) = cols[n]; n++;
		}
	}
	if (top >= 0)
	{
		if(bottom >=0)
		{
			n = 0;
			MPI_Recv(&rows[0], row_count, MPI_FLOAT, top, TAG, comm, MPI_STATUS_IGNORE);
			for (int j = nei_ystart; j < row_count + nei_ystart; j++)
			{
				*(uu + j) = rows[n]; n++;
			}
		}
		if(bottom <0) //nếu là hàng cuối, recv 2 hàng
		{
			n = 0;
			MPI_Recv(&rows_dim[0], 2*row_count, MPI_FLOAT, top, TAG, comm, MPI_STATUS_IGNORE);
			for (int i = 0; i < 2; i++)
			{
				for (int j = nei_ystart; j < row_count + nei_ystart; j++)
				{
					*(uu + i*ny + j) = rows_dim[n]; n++;
				}
			}
		}
	}
	if (bottom >= 0)
	{
		MPI_Recv(&rows[0], row_count, MPI_FLOAT, bottom, TAG, comm, MPI_STATUS_IGNORE);
		n = 0;
		for (int j = nei_ystart; j < row_count + nei_ystart; j++)
		{
			*(uu + (nx-1)*ny + j) = rows[n]; n++;
		}
	}
}

void Collect_Array (int dimx, int dimy, int TAG, int *coords, int my_rank, MPI_Comm comm)
{
	
}
int isPerfectSquare(int number) {
    int squareRoot = sqrt(number);
    return (squareRoot * squareRoot == number);
	}
	
int Cal_dim_length(int num_tasks)
{
	
	int dimx, dimy;
	int tempx, tempy;
	if (isPerfectSquare(num_tasks))
	{
		tempx = sqrt(num_tasks);
		tempy = tempx;
	}		
	else
	{	

		dimx = 1; dimy = 1000000; 
		tempx = dimx; tempy = dimy;
		while (dimx <= num_tasks/2) 
		{	
			if (num_tasks%(dimx) == 0)
			{
				dimy = num_tasks/dimx;
				if (abs(dimx-dimy) < abs(tempx-tempy))
				{tempx = dimx; 
				tempy = dimy;}
			}
			dimx = dimx + 1;
		}
	}
	return tempx;
}


void Check()
{
	printf("\nCheck sucessfully\n");
}


void print_array(float *uu, int nx, int ny)
{
	printf("\n");
	for (int i =0; i< nx; i++)
	{
		for (int j = 0; j< ny; j++)
		{
			printf ("%.05f ", *(uu + i * ny + j));
		}
	printf("\n");
	}
	printf("\n");
}

void initialize_test(int nx, int ny, float *uu, float value) 
{
    int i, j;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < ny; j++) {
            *(uu + i * ny + j) = value;
            value = value + 1;
        }
    }
}

