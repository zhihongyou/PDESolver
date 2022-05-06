#ifndef FIELDCLASSGPU_FIELDFUNCTIONSGPU_CU
#define FIELDCLASSGPU_FIELDFUNCTIONSGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>

using namespace std;

// ----------------------------------------------------------------------
__global__ void getD1xGPUCore(double* d1x, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        d1x[idx]=FDM_ptrs[FDM_idx]->d1x(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void getD1yGPUCore(double* d1y, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        d1y[idx]=FDM_ptrs[FDM_idx]->d1y(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void getD2xGPUCore(double* d2x, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        d2x[idx]=FDM_ptrs[FDM_idx]->d2x(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void getD2yGPUCore(double* d2y, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        d2y[idx]=FDM_ptrs[FDM_idx]->d2y(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void getD1x1yGPUCore(double* d1x1y, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        d1x1y[idx]=FDM_ptrs[FDM_idx]->d1x1y(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void getLaplaceGPUCore(double* laplace, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        laplace[idx]=FDM_ptrs[FDM_idx]->laplace(f_t,idx,di,dj,dx,dy);
    };
};

// ----------------------------------------------------------------------
__global__ void getBiLaplaceGPUCore(double* bi_laplace, double* f_t, FiniteDifference** FDM_ptrs, int FDM_idx, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        bi_laplace[idx]=FDM_ptrs[FDM_idx]->bi_laplace(f_t,idx,di,dj,dx,dy);
    };
};


#endif
