#ifndef POISSONEQFIELDCLASSGPU_CU
#define POISSONEQFIELDCLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "poissonEqFieldclass.h"

using namespace std; 

// ===============================================================
__global__ void getPhiGPU(cufftDoubleComplex* phi_complex, double* phi, double* f, double* poisson_k2, int Nx, int Ny, int Nbx, int Nby, int getType) {
    // This is used to solve the Poisson equation:
    //   \nabla^2 \phi = f
    int i=blockIdx.x;
    int j=threadIdx.x;  
    int idx1=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;
    int idx = i*Nx+j;

    if (j<Nx && i<Ny) {
        if (getType==0) {
            phi_complex[idx].x=f[idx1];
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


#endif
