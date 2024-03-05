#ifndef LIVINGLCPOLARFIELDCLASSGPU_CU
#define LIVINGLCPOLARFIELDCLASSGPU_CU
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "LivingLCPolarFieldClass.h"

using namespace std; 

//================================================================
__global__ void getLivingLCPxPyThetaGPU(double* Pxx, double* Pxy, double* px, double* py, double* theta, double* theta_old, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;

    if (i<Ny && j<Nx) {
      theta_old[idx] = theta[idx];
      theta[idx] = atan2(Pxy[idx],Pxx[idx]);
      double p = 2*sqrt(Pxx[idx]*Pxx[idx]+Pxy[idx]*Pxy[idx]);
      px[idx] = p*cos(theta[idx]);
      py[idx] = p*sin(theta[idx]);
    };
}


//================================================================
__global__ void getLivingLCFlipGPU(double* theta_old, double* theta, double* flip, int Nx, int Ny, int Nbx, int Nby) {
    int i=blockIdx.x;
    int j=threadIdx.x;
    int idx=(blockDim.x+2*Nbx)*(i+Nby)+j+Nbx;

    if (i<Ny && j<Nx) {
      if (abs(theta_old[idx]-theta[idx]) > 1.57) {
          flip[idx] = -1*flip[idx];
      };
    };
}


#endif
