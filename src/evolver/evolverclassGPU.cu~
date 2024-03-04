#ifndef EVOLVERCLASSGPU_CU
#define EVOLVERCLASSGPU_CU

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include <iostream> 
#include <vector>

using namespace std;


// ----------------------------------------------------------------------
__global__ void updateRHSCoreGPU(rhsPtrs rhs_ptrs, double* rhs_temp, double* lhs_temp, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
    
    if (i<Nx && j<Ny) {
        rhs_temp[idx]=0;
        lhs_temp[idx]=0;
        int i_func=0;
        for (int i_term=0; i_term<rhs_ptrs.num_terms[0]; i_term++) {
            double temp=rhs_ptrs.prefactors[i_term];
            for (int i_func1=0; i_func1<rhs_ptrs.num_funcs_1term[i_term]; i_func1++) {
                temp=temp*rhs_ptrs.f_func_ptrs[i_func][idx];
                // temp=1;
                i_func+=1;
            };
            rhs_temp[idx]+=temp;                
        };

        for (int i_term=rhs_ptrs.num_terms[0]; i_term<rhs_ptrs.num_terms[0]+rhs_ptrs.num_terms[1]; i_term++) {
            double temp=rhs_ptrs.prefactors[i_term];
            for (int i_func1=0; i_func1<rhs_ptrs.num_funcs_1term[i_term]-1; i_func1++) {
                temp=temp*rhs_ptrs.f_func_ptrs[i_func][idx];
                i_func+=1;
            };
            lhs_temp[idx]+=temp;
        };

    };
};


// ----------------------------------------------------------------------
__global__ void fieldUpdateGPUCore(double* f_new, double* f_old, double* f_rhs, double* f_lhs, double time_step_t, int Nx, int Ny, int Nbx, int Nby) {
    int i=threadIdx.x;
    int j=blockIdx.x;    
    int idx=(blockDim.x+2*Nbx)*(j+Nby)+i+Nbx;
    
    if (i<Nx && j<Ny) {
        f_new[idx]=(f_old[idx]+f_rhs[idx]*time_step_t)/(1+f_lhs[idx]*time_step_t);
    };
};


// ----------------------------------------------------------------------

#endif
