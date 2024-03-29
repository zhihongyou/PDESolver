#ifndef USERDEFINEDFUNCTION_H
#define USERDEFINEDFUNCTION_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "/home/you/PDESolver/src/field/fieldFunction.h"

// ===============================================================
// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double getPressure(double * f, int idx, FFuncArgs f_func_args) {
    double func=0;
    if (f[idx]>f_func_args.f_func_arg2) {
        func=f_func_args.f_func_arg1*(f[idx]-f_func_args.f_func_arg2);
    }
    return func;
};
__device__ FFuncType getPressure_dev=getPressure;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double oneOverPhi(double * f, int idx, FFuncArgs f_func_args) {
    double func=f[idx];
    if (f[idx]>f_func_args.f_func_arg1) {
        func=f_func_args.f_func_arg1;
    }
    return func;
};
__device__ FFuncType oneOverPhi_dev=oneOverPhi;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double atan2F(double * f, int idx, FFuncArgs f_func_args) {
   
    return atan2(f[idx],f_func_args.field2[idx]);
};
__device__ FFuncType atan2F_dev=atan2F;


// ---------------------------------------------------------------
CUDA_CALLABLE_MEMBER double toatan2F(double * f, int idx, FFuncArgs f_func_args) {
    double theta;
    if (f[idx] == 4){
        theta = 0;
    } else if(f[idx] == 2){
        theta = 0;
    } else if(f[idx] == -2){
        theta = 3.14159265358979323846;
    } else if(f[idx] == -4){
        theta = -3.14159265358979323846;
    } else if(f[idx] == 3){
        theta = 0;
    } else if(f[idx] == -3){
        theta = 0;
    }
    return theta;
};
__device__ FFuncType toatan2F_dev=toatan2F;


// ---------------------------------------------------------------
void addUserDefinedFuncs () {
    f_func_map_all[{"getPressure",""}]=getPressure;
    f_func_map_all[{"oneOverPhi",""}]=oneOverPhi;
    f_func_map_all[{"atan2F",""}]=atan2F;
    f_func_map_all[{"toatan2F",""}]=toatan2F;    
    f_func_map_all_dev[{"getPressure",""}]=getFFuncDevPtr(&getPressure_dev);
    f_func_map_all_dev[{"oneOverPhi",""}]=getFFuncDevPtr(&oneOverPhi_dev);
    f_func_map_all_dev[{"atan2F",""}]=getFFuncDevPtr(&atan2F_dev);
    f_func_map_all_dev[{"toatan2F",""}]=getFFuncDevPtr(&toatan2F_dev);    
};


#endif
