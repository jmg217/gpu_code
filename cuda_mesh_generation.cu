#include <cuda.h>
#include <iostream>
#include <math.h>
#include <ctime>
#include <cmath>
#include "enum_header.h"
#include <unistd.h>
#include <stdio.h>

/* we need these includes for CUDA's random number stuff */
#include <curand.h>
#include <curand_kernel.h>


void meshgeneration(int b, int num_assets, int m_int){

int X_N=(m_int) * b * (num_assets);
double* X_device;

}

