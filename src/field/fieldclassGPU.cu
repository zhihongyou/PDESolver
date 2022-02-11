#ifndef FIELDCLASSGPU_CU
#define FIELDCLASSGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include "fieldclass.h"
#include "../utility/finiteDifferenceCentralO2Isotropic.h"

using namespace std; 



// ----------------------------------------------------------------------
__global__ void getLaplaceGPUCore(double* laplace, double* f_t, int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    int di=1;
    int dj=Nx+2*Nbx;    

    if (i<Nx && j<Ny) {
        laplace[idx]=FDMCentralO2I::laplace(f_t,idx,di,dj,dx,dy);
    };
};


// ----------------------------------------------------------------------
__global__ void applyBounCondPeriAnyGPU(double* f_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;
    int idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
    int dj=Nx+2*Nbx;
    int idx1;

    if (j<Nby) {   		// Bottom rows
        idx1=idx+Ny*dj;
        f_t[idx]=f_t[idx1];
    } else if (j>Ny-1-Nby) {	// Top rows
        idx1=idx-Ny*dj;
        f_t[idx]=f_t[idx1];
    };

    if (i<Nbx) {                 // Left columns
        idx1=idx+Nx;
        f_t[idx]=f_t[idx1];    
        if (j<Nby) {               // Left bottom corner.
            idx1=idx+Ny*dj+Nx;
            f_t[idx]=f_t[idx1];
        } else if (j>Ny-1-Nby) {	  // Left top corner.
            idx1=idx-Ny*dj+Nx;
            f_t[idx]=f_t[idx1];
        };
    } else if (i>Nx-1-Nbx) {   // Right columns
        idx1=idx-Nx;
        f_t[idx]=f_t[idx1];
        if (j<Nby) {               // Right bottom corner.
            idx1=idx+Ny*j-Nx;
            f_t[idx]=f_t[idx1];
        } else if (j>Ny-1-Nby) {	 // Right top corner.
            idx1=idx-Ny*dj-Nx;
            f_t[idx]=f_t[idx1];
        };
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
