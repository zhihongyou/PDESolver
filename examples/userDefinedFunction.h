#ifndef USERDEFINEDFUNCTION_H
#define USERDEFINEDFUNCTION_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "/home/you/Research/codes/PDESolver/src/field/fieldFunction.h"

// ===============================================================
// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double getPressure(double * f, int idx, FFuncArgs f_func_args) {
    double func=0;
    if (f[idx]>f_func_args.f_func_arg2) {
        func=f_func_args.f_func_arg1*(f[idx]-f_func_args.f_func_arg2);
    }
    return func;
};

// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double oneOverPhi(double * f, int idx, FFuncArgs f_func_args) {
    double func=f[idx];
    if (f[idx]>f_func_args.f_func_arg1) {
        func=f_func_args.f_func_arg1;
    }
    return func;
};

// Define device functions========================================
// ---------------------------------------------------------------
__device__ FFuncType getPressure_dev=getPressure;
__device__ FFuncType oneOverPhi_dev=oneOverPhi;


// ---------------------------------------------------------------
void addUserDefinedFuncs () {
    f_func_map_all[{"getPressure",""}]=getPressure;
    f_func_map_all[{"oneOverPhi",""}]=oneOverPhi;
    f_func_map_all_dev[{"getPressure",""}]=getFFuncDevPtr(&getPressure_dev);
    f_func_map_all_dev[{"oneOverPhi",""}]=getFFuncDevPtr(&oneOverPhi_dev);
};


#endif
