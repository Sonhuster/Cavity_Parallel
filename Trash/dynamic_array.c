#include "dynamic_array.h"
#include <stdlib.h>
#include <stdio.h>

float*** createArray3D(int nz, int rows, int cols, float value) {
    float*** array3D = (float***)malloc(nz * sizeof(float**));
    for (int i = 0; i < nz; ++i) {
        array3D[i] = (float**)malloc(rows * sizeof(float*));
        for (int j = 0; j < rows; ++j) {
            array3D[i][j] = (float*)malloc(cols * sizeof(float));
            for (int k = 0; k < cols; ++k) {
                array3D[i][j][k] = value;
				value ++;
            }
        }
    }
    return array3D;
}

void freeArray3D(float*** array3D, int nz, int rows) {
    for (int i = 0; i < nz; ++i) {
        for (int j = 0; j < rows; ++j) {
            free(array3D[i][j]);
        }
        free(array3D[i]);
    }
    free(array3D);
}