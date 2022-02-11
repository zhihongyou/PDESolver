#ifndef EVOLVERCLASS_CU
#define EVOLVERCLASS_CU

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream> 
#include <vector>
#include "evolverclass.h"

using namespace std;


// ----------------------------------------------------------------------
__global__ void addRHSTermGPU(rhs_term rhs_term_t, double* f_rhs_temp_ptr, double** f_func_ptrs_dev, int N_funcs, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;

    if (i<Nx && j<Ny) {
        double temp=rhs_term_t.prefactor;
        for (int i_func=0; i_func<N_funcs;i_func++) {
            temp=temp*f_func_ptrs_dev[i_func][idx];
        };
        f_rhs_temp_ptr[idx]+=temp;
    };
};


// ----------------------------------------------------------------------
__global__ void fieldUpdateGPUCore(double* f_new, double* f_old, double* f_rhs, double time_step_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;

    if (i<Nx && j<Ny) {
        f_new[idx]=f_old[idx]+f_rhs[idx]*time_step_t;
    };
};



#endif
