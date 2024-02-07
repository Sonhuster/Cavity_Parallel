int *Send_back (float *array, float **sub_array, int size_0, int size_1,  int dimx, int dimy, int *coords, int my_rank, MPI_Comm comm)
{
	int centerx, centery;
	int dim[2]; 
	dim[0] = dimx; dim[1] = dimy;
	//printf("\nRank:%d. Position:  %d, %d \n",my_rank, coords[0], coords[1]);
	centerx = coords[0];
	centery = coords[1];
	
	
	int xstart, xend, xsize,
		ystart, yend, ysize;
			
			
	if(centerx==0)
	{
		xstart = 0; xend = size_0 - 2; xsize = xend - xstart + 1;
	}
	if(centery==0)
	{
		ystart = 0; yend = size_1 - 2; ysize = yend - ystart + 1;
	}
	if(centerx>0 && centerx<dim[0]-1)
	{
		xstart = 1; xend = size_0 - 2; xsize = xend - xstart + 1;
	}
	if(centery>0 && centery<dim[0]-1)
	{
		ystart = 1; yend = size_1 - 2; ysize = yend - ystart + 1;
	}
	if(centerx==dim[0]-1)
	{
		xstart = 2; xend = size_0 - 1; xsize = xend - xstart + 1;
	}
	if(centery==dim[1]-1)
	{
		ystart = 2; yend = size_1 - 1; ysize = yend - ystart + 1;
	}
	
	int m = 0, n = 0;
	for (int i = xstart; i <= xend; i++)
	{
		for (int j = ystart; j <= yend; j++)
		{
			sub_array[m][n] = *(array + i*size_1 + j);
			n++;
		}
		m++;
		n = 0;
	}
	int *size = (int *)malloc(2*sizeof(int*));
	size[0] = xsize; size[1] = ysize;
	return size;
		
}


int *Count_Send_back(int size_px, int size_py, int dimx, int dimy, int *coords, int TAG)
{
	int *size, centerx, centery;
	int xstart, ystart, xend, yend, xsize, ysize;
	size[0] = size_px; size[1] = size_py;
	if (TAG == count_u) {size [0] = size_px -1;}
	if (TAG == count_v)	{size [1] = size_py -1;}
	centerx = coords[0];
	centery = coords[1];
	
	if(centerx==0)
	{
		xstart = 0; xend = size[0] - 2; xsize = xend - xstart + 1;
	}
	if(centery==0)
	{
		ystart = 0; yend = size[1] - 2; ysize = yend - ystart + 1;
	}
	if(centerx>0 && centerx<dimx-1)
	{
		xstart = 1; xend = size[0] - 2; xsize = xend - xstart + 1;
	}
	if(centery>0 && centery<dimy-1)
	{
		ystart = 1; yend = size[1] - 2; ysize = yend - ystart + 1;
	}
	if(centerx==dimx-1)
	{
		xstart = 2; xend = size[0] - 1; xsize = xend - xstart + 1;
	}
	if(centery==dimy-1)
	{
		ystart = 2; yend = size[1] - 1; ysize = yend - ystart + 1;
	}
	return size;
}


//7-2 Collected array
		MPI_Send(&u_sendback, ux_size*uy_size, MPI_FLOAT, MASTER, trans_u, comm);
		MPI_Send(&v_sendback, vx_size*vy_size, MPI_FLOAT, MASTER, trans_v, comm);
		MPI_Send(&p_sendback, px_size*py_size, MPI_FLOAT, MASTER, trans_p, comm);
		
		if(my_rank == MASTER)
		{
			int 	local_nrows = (NROWS%dim[0] == 0) ?NROWS/dim[0] : NROWS/dim[0] + 1; //const
			int 	local_ncols = (NCOLS%dim[1] == 0) ?NCOLS/dim[1] : NCOLS/dim[1] + 1; //const
			float 	***Collected_P = createArray3D(local_nrows, local_ncols, (dim[0]-1)*(dim[1]-1), 0),
					***Collected_U = createArray3D(local_nrows, local_ncols, (dim[0]-1)*(dim[1]-1), 1), 
					***Collected_V = createArray3D(local_nrows, local_ncols, (dim[0]-1)*(dim[1]-1), 0);
			int 	collected_dest = 0, count = 0;
			
			for (int ix = 0; ix < dim[0] - 1; ix ++)
			{
				for (int iy = 0; iy < dim[1] - 1; iy ++)
				{
					collected_dest = ix*(dim[1] - 1) + iy;
					//Creat Datatype
					int Collected_recv[2] = {local_nrows,local_ncols}, 
						subsizes[2] = {local_nrows,local_ncols},
						ar_start[2]= {0,0};
					MPI_Datatype ar_subsizes;
					MPI_Type_create_subarray(2, Collected_recv, subsizes, ar_start, MPI_ORDER_C, MPI_FLOAT, &ar_subsizes);
					MPI_Type_commit (&ar_subsizes);
					
					MPI_Recv(&Collected_U[0][0][count], 1, ar_subsizes, collected_dest, trans_u, comm, MPI_STATUS_IGNORE);
					MPI_Recv(&Collected_V[0][0][count], 1, ar_subsizes, collected_dest, trans_v, comm, MPI_STATUS_IGNORE);
					MPI_Recv(&Collected_P[0][0][count], 1, ar_subsizes, collected_dest, trans_p, comm, MPI_STATUS_IGNORE);
					count ++;
					MPI_Type_free (&ar_subsizes);
				}
			}
			printf("\n");
			for (int i =0; i< local_nrows; i++)
			{
				for (int j = 0; j< local_ncols; j++)
				{
					printf ("%.0f ", Collected_U[1][j][i]);
				}
				printf("\n");
			}
			printf("\n");
		}
	
	