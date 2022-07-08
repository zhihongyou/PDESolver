#ifndef FIELDCLASSGPU_CU
#define FIELDCLASSGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>

using namespace std;


// ----------------------------------------------------------------------
__global__ void fieldCopy(double* f_now, double* f_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;

    if (i<Nx && j<Ny) {
        f_now[idx]=f_t[idx];
    };
};


// ----------------------------------------------------------------------
void setFieldConstCPUCore(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nbx)*(Nx+2*Nbx)+i+Nbx;
            f_t[idx]=f_val;
        };                
    };
};

// ----------------------------------------------------------------------
__global__ void setFieldConstGPUCore(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;
    int idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
    
    if (i<Nx && j<Ny) {
        f_t[idx]=f_val;
    };
};


#endif
