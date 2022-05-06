#ifndef FIELDCLASS_FIELDFUNCTIONS_GENERICGPU_CU
#define FIELDCLASS_FIELDFUNCTIONS_GENERICGPU_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>

using namespace std;


// ----------------------------------------------------------------------
template<class FDM_class>
__global__ void getFFuncGPUCore1(double* f_func_ptr, double* f_t, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), int Nx, int Ny, int Nbx, int Nby, double dx, double dy) {
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
template<typename T>
__global__ void getFFuncGPUCore(T* f_func_ptr, double* f_t, T f_func(double*,int), int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;

    if (i<Nx && j<Ny) {
        f_func_ptr[idx]=f_func(f_t,idx);
        // f_func_ptr[idx]=FieldFunction::sinF(f_t,idx);
    };
};


#endif
