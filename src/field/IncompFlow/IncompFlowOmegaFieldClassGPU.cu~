#ifndef INCOMPRFLOWOMEGAFIELDCLASSGPU_CU
#define INCOMPRFLOWOMEGAFIELDCLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "IncompFlowClass.h"

using namespace std; 

//================================================================
__global__ void getIncompFlowPhi2VGPU(double* phi, double* vx, double* vy, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*f_func_args.Nbx)*(i+f_func_args.Nby)+j+f_func_args.Nbx;

    if (i<f_func_args.Ny && j<f_func_args.Nx) {
        vx[idx] = d1y(phi,idx,f_func_args);
        vy[idx] = -d1x(phi,idx,f_func_args);
    };
}


#endif
