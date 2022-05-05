#ifndef FIELDCLASSGPU_CU
#define FIELDCLASSGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include "fieldclass.h"

using namespace std;

// ----------------------------------------------------------------------
template<class FDM_class>
__global__ void getFFuncGPUCore(double* f_func_ptr, double* f_t, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        f_func_ptr[idx]=(FDM_scheme.*f_func)(f_t,idx,di,dj,dx,dy);
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

// ----------------------------------------------------------------------
__global__ void getFNowGPU(double* f_now, double* f_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;    

    if (i<Nx && j<Ny) {
        f_now[idx]=f_t[idx];
    };
};


// ----------------------------------------------------------------------
__global__ void applyBounCondPeriAnyGPU(double* f_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(Nx+2*Nbx)*(i+Nby)+j+Nbx;
    int dj=Nx+2*Nbx;
    int idx1;

    if (i<Nby) {   		// Bottom rows
        idx1=idx+Ny*dj;
        f_t[idx1]=f_t[idx];
    } else if (i>Ny-1-Nby) {	// Top rows
        idx1=idx-Ny*dj;
        f_t[idx1]=f_t[idx];
    }

    if (j<Nbx) {                 // Left columns
        idx1=idx+Nx;
        f_t[idx1]=f_t[idx];    
        if (i<Nby) {               // Left bottom corner.
            idx1=idx+(Nx+2*Nbx)*Ny+Nx;
            f_t[idx1]=f_t[idx];
        } else if (i>Ny-1-Nby) {	  // Left top corner.
            idx1=idx-(Nx+2*Nbx)*Ny+Nx;
            f_t[idx1]=f_t[idx];
        }
    } else if (j>Nx-1-Nbx) {   // Right columns
        idx1=idx-Nx;
        f_t[idx1]=f_t[idx];
        if (i<Nby) {               // Right bottom corner.
            idx1=idx+(Nx+2*Nbx)*Ny-Nx;
            f_t[idx1]=f_t[idx];
        } else if (i>Ny-1-Nby) {	 // Right top corner.
            idx1=idx-(Nx+2*Nbx)*Ny-Nx;
            f_t[idx1]=f_t[idx];
        }			       
    }    
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
