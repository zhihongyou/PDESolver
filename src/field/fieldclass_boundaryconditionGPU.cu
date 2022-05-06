#ifndef FIELDCLASSGPU_BOUNDARYCONDITIONGPU_CU
#define FIELDCLASSGPU_BOUNDARYCONDITIONGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>

using namespace std;


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


#endif
