#ifndef DYNAMIC_ARRAY_H
#define DYNAMIC_ARRAY_H
#include <stdio.h>

float*** createArray3D(int rows, int cols, int nz, float value);
void freeArray3D(float*** array3D, int rows, int cols);

#endif