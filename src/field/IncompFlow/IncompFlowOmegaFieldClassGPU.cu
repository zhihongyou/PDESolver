#ifndef INCOMPFLOWOMEGAFIELDCLASSGPU_CU
#define INCOMPFLOWOMEGAFIELDCLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "IncompFlowClass.h"

using namespace std; 

// ===============================================================
__global__ void getIncompFlowStreamGPU(cufftDoubleComplex* phi_complex, double* phi, double* omega, double* poisson_k2, int Nx, int Ny, int Nbx, int Nby, int getType) {
    int i=blockIdx.x;
    int j=threadIdx.x;  
    int idx1=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
    int idx = i*Nx+j;

    if (j<Nx && i<Ny) {
        if (getType==0) {
            // Get vorticity
            phi_complex[idx].x=-omega[idx1];
            phi_complex[idx].y=0;
        } else if (getType==1) {
            if (poisson_k2[idx] == 0){
                // Setting ifft of w to 0 when wavenumber is 0
                phi_complex[idx].x = 0;
                phi_complex[idx].y = 0;
            } else {
                // Get stream function in Fourier space.
                phi_complex[idx].x = -phi_complex[idx].x/(poisson_k2[idx]);
                phi_complex[idx].y = -phi_complex[idx].y/(poisson_k2[idx]);
            };
        } else if (getType==2) {
            // Get real value of stream function.
            phi[idx1]=phi_complex[idx].x/(Nx*Ny);
        };
    };
}

//================================================================
__global__ void getIncompFlowVCoreGPU(double* phi, double* vx, double* vy, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*f_func_args.Nbx)*(i+f_func_args.Nby)+j+f_func_args.Nbx;

    if (i<f_func_args.Ny && j<f_func_args.Nx) {
        vx[idx] = d1y(phi,idx,f_func_args);
        vy[idx] = -d1x(phi,idx,f_func_args);
    };
}


#endif
