#include "initialize.h"


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
            *(vv + i * (ny-1) + j) = 0;
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
