#include "Update.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NROWS 		11
#define NCOLS		11


void update (float *uu, float *vv, float *pp, int nx, int ny, int top, int bottom, int left, int right);
void initialize(int nx, int ny, float *uu, float *vv, float *pp);
void print_array(float *uu, int nx, int ny);

int main ()
{
	float 	u[NROWS-1][NCOLS], 
			v[NROWS][NCOLS-1], 
			p[NROWS][NCOLS];
	initialize(NROWS, NCOLS, &u[0][0], &v[0][0], &p[0][0]);
	
	
	int Iteration = 1000;
	
	for (int it = 0; it <= Iteration; it ++)
	{
	update (&u[0][0], &v[0][0], &p[0][0], NROWS, NCOLS, -1, -1, -1, -1);
	}
	print_array(&u[0][0], NROWS-1, NCOLS);
	
	return 0;
}








void update (float *uu, float *vv, float *pp, int nx, int ny, int top, int bottom, int left, int right)
{
	float 	hx = 0.01, hy = 0.01,
			delta = 4.5, dt = 0.001, Re = 100;
	int nux, nvy;
	nux = nx, nvy = ny;
	if (bottom < 0) {nux = nx -1;}
	if (right < 0) {nvy = ny -1;}
	//printf("\n nux: %d, ny: %d\n", nux,ny);
	
	float u[nux][ny], 		v[nx][nvy], 		p[nx][ny];
	float u_new[nux][ny], 	v_new[nx][nvy],		p_new[nx][ny];
	
	for (int i = 0; i < nux; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            u[i][j] = 0;
			u[i][j] = *(uu + i*ny + j);
			u_new[i][j] = u[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < nvy; j++) 
		{
            v[i][j] = 0;
			v[i][j] = *(vv + i*nvy + j);
			v_new[i][j] = v[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            p[i][j] = 0;
			p[i][j] = *(pp + i*ny + j);
			p_new[i][j] = p[i][j];
        }
    }

	float pressure = 0, advection_x = 0, advection_y = 0, diffusion = 0;
	//x-momentum eq. - Interior
	for (int i = 1; i < nux-1; i++)
	{
		for (int j = 1; j < ny-1; j++)
		{
			pressure = -(p[i+1][j] - p[i][j])/hx;
			diffusion = (1/Re)*(
				(u[i+1][j]-2*u[i][j]+u[i-1][j])/(hx*hx)+
				(u[i][j+1]-2*u[i][j]+u[i][j-1])/(hy*hy));
			advection_x = (
				pow((0.5*(u[i+1][j]+u[i][j])),2)-
				pow((0.5*(u[i-1][j]+u[i][j])),2))/hx;  
			advection_y = (
				(0.25*(u[i][j-1]+u[i][j])*(v[i][j-1]+v[i+1][j-1]))-
				(0.25*(u[i][j+1]+u[i][j])*(v[i][j]+v[i+1][j])))/hy;
			u_new[i][j]= u[i][j]+dt*(pressure + diffusion - (advection_x + advection_y));
		}
	}

	//x-momentum eq. - Boundary
	for (int j=0; j<ny; ++j)
	{
		if(top < 0) 	{u_new[0][j] = 0.0;}
		if(bottom < 0) 	{u_new[nux-1][j] = 0.0;}
	}
	
	for (int i = 0 ; i<nux; ++i)
	{
		if(left < 0) 	{u_new[i][0] = 2 - u_new[i][1];}
		if(right < 0) 	{u_new[i][ny-1] = - u_new[i][ny-2];}
	}
		
	
	//y-momentum eq. - Interior
	for (int i = 1; i < nx-1; i++)
	{
		for (int j = 1; j < nvy-1; j++)
		{
			pressure = -(p[i][j] - p[i][j+1])/hy;
			diffusion = (1/Re)*(
            (v[i][j+1]-2*v[i][j]+v[i][j-1])/(hy*hy)+
            (v[i+1][j]-2*v[i][j]+v[i-1][j])/(hx*hx));
			advection_y = (
            pow((0.5*(v[i][j-1]+v[i][j])),2)-
            pow((0.5*(v[i][j+1]+v[i][j])),2))/hy; 
			advection_x = (
            (0.25*(v[i+1][j]+v[i][j])*(u[i][j+1]+u[i][j]))-
            (0.25*(v[i-1][j]+v[i][j])*(u[i-1][j]+u[i-1][j+1])))/hx;
			v_new[i][j]= v[i][j]+dt*(pressure + diffusion - (advection_x + advection_y));
		}
	}
	
	//y-momentum eq. - Boundary
	for (int j=0; j<nvy; ++j)
	{
		if(top < 0) 		v_new[0][j] = - v_new[1][j];
		if(bottom < 0)	v_new[nx-1][j] = - v_new[nx-2][j];
	}

	for (int i=0; i<nx; ++i)
	{
		if(left < 0)		v_new[i][0] = 0;
		if(right < 0)	v_new[i][nvy-1] = 0;
	}
	
	//Continuity equation - Interior 
	for (int i = 1; i < nx-1; i++)
	{
		for (int j = 1; j < ny-1; j++)
		{
			p_new[i][j] = p[i][j]-
            dt*delta*(u[i][j]-u[i-1][j])/hx -
            dt*delta*(v[i][j-1]-v[i][j])/hy;
		}
	}
	
	//Continuity eq. - Boundary
	for (int i=1; i<nx-1; ++i)
	{
		if(left < 0)		p_new[i][0] = p_new[i][1];
		if(right < 0)	p_new[i][ny-1] = p_new[i][ny-2];
	}

	for (int j=0; j<ny; ++j)
	{
		if(top < 0)		p_new[0][j] = p_new[1][j];
		if(bottom < 0)	p_new[nx-1][j] = p_new[nx-2][j];
	}
	
	/*update matrix*/
	for (int i = 0; i < nux; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
            u[i][j] = u_new[i][j];
			*(uu + i*ny + j) = u[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < nvy; j++) 
		{
			v[i][j] = v_new[i][j];
			*(vv + i*nvy + j) = v[i][j];
        }
    }
	
	for (int i = 0; i < nx; i++) 
	{
        for (int j = 0; j < ny; j++) 
		{
			p[i][j] = p_new[i][j];
			*(pp + i*ny + j) = p[i][j];
        }
    }

	/*End update matrix*/
}

void initialize(int nx, int ny, float *uu, float *vv, float *pp)
{
    int i, j;
    for (i = 0; i < nx-1; i++) {
        for (j = 0; j < ny; j++) 
		{
            *(uu + i * ny + j) = 0;
			*(uu + i * ny) = 2;
        }
    }
	
	for (i = 0; i < nx; i++) 
	{
        for (j = 0; j < ny-1; j++) 
		{
            *(vv + i * ny + j) = 0;
        }
    }
	
	for (i = 0; i < nx; i++) 
	{
        for (j = 0; j < ny; j++) 
		{
            *(pp + i * ny + j) = 0;
        }
    }
	*(pp + (nx-1)*ny + ny -1) = 1;
	
}

void print_array(float *uu, int nx, int ny)
{
	printf("\n");
	for (int i =0; i< nx; i++)
	{
		for (int j = 0; j< ny; j++)
		{
			printf ("%.04f ", *(uu + i * ny + j));
		}
	printf("\n");
	}
	printf("\n");
}