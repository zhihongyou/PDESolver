#ifndef FIELDCLASS_FIELDFUNCTIONSGPU_CU
#define FIELDCLASS_FIELDFUNCTIONSGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>

using namespace std;


// ----------------------------------------------------------------------
__global__ void getFFuncGPUCore(double* f_func_ptr, double* f_t, FFuncType f_func, FFuncArgs f_func_args) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*f_func_args.Nbx)*(j+f_func_args.Nby)+i+f_func_args.Nbx;

    if (i<f_func_args.Nx && j<f_func_args.Ny) {
        f_func_ptr[idx]=f_func(f_t, idx, f_func_args);
    };
};


#endif
