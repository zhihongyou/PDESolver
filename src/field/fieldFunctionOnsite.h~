#ifndef FIELDFUNCTIONONSITE_H
#define FIELDFUNCTIONONSITE_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "fieldFunction.h"

// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double oneF(double * f, int idx, FFuncArgs f_func_args) {
    return f[idx];
};

__device__ FFuncType oneF_dev=oneF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double oneOverF(double * f, int idx, FFuncArgs f_func_args) {
    return 1/(f[idx]+0.00000000000001);
};

__device__ FFuncType oneOverF_dev=oneOverF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double sinF(double * f, int idx, FFuncArgs f_func_args) {
    return sin(f_func_args.f_func_arg1*f[idx]+f_func_args.f_func_arg2);
};

__device__ FFuncType sinF_dev=sinF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double cosF(double * f, int idx, FFuncArgs f_func_args) {
  //return cos(f_func_args.f_func_arg1*f[idx]+f_func_args.f_func_arg2);
    return cos(f[idx]);
};

__device__ FFuncType cosF_dev=cosF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double expF(double * f, int idx, FFuncArgs f_func_args) {
    return exp(f_func_args.f_func_arg1*f[idx]+f_func_args.f_func_arg2);
};

__device__ FFuncType expF_dev=expF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double absF(double * f, int idx, FFuncArgs f_func_args) {
    return abs(f[idx]);
};

__device__ FFuncType absF_dev=absF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double signF(double * f, int idx, FFuncArgs f_func_args) {
    return (double(0) < f[idx]) - (f[idx] < double(0));
};

__device__ FFuncType signF_dev=signF;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double powF(double * f, int idx, FFuncArgs f_func_args) {
    return pow(f[idx],f_func_args.f_func_arg1);
};

__device__ FFuncType powF_dev=powF;


#endif
