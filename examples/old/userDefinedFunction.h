#ifndef USERDEFINEDFUNCTION_H
#define USERDEFINEDFUNCTION_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "/home/you/Research/codes/PDESolver/src/field/fieldFunction.h"


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double myFunc(double * f, int idx, FFuncArgs f_func_args) {
    double func=0;
    if (f[idx]>1) {
        func=f[idx]-1;
    }
    return func;
};

__device__ FFuncType myFunc_dev=myFunc;


// ----------------------------------------------------------------------
void addUserDefinedFuncs () {
    f_func_map_all[{"myFunc",""}]=myFunc;
    f_func_map_all_dev[{"myFunc",""}]=getFFuncDevPtr(&myFunc_dev);
};


#endif
