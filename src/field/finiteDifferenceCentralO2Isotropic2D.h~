#ifndef FINITEDIFFERENCECENTRALO2ISOTROPIC2D_H
#define FINITEDIFFERENCECENTRALO2ISOTROPIC2D_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "fieldFunction.h"

/* The name "CO2I2D" stands for "centeral difference,
    2nd order, Isotropic, 2D".*/
    
// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1xCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    /* di-y, dj-x */
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(12.0*dx)*(
        (f[idx+di+dj] - f[idx+di-dj])
        +4*(f[idx+dj] - f[idx-dj])
        +(f[idx-di+dj] - f[idx-di-dj]));
};

__device__ FFuncType d1xCO2I2D_dev=d1xCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1yCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(12.0*dy)*(
        (f[idx+di+dj] - f[idx-di+dj])
        +4*(f[idx+di] - f[idx-di])
        +(f[idx+di-dj] - f[idx-di-dj]));
};

__device__ FFuncType d1yCO2I2D_dev=d1yCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1x1yCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(4.0*dx*dy)*(f[idx+di+dj]-f[idx-di+dj]
    -f[idx+di-dj]+f[idx-di-dj]);
};

__device__ FFuncType d1x1yCO2I2D_dev=d1x1yCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d2xCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(12.0*dx*dx)*(
        (f[idx+di+dj]-2*f[idx+di]+f[idx+di-dj])
        +10*(f[idx+dj]-2*f[idx]+f[idx-dj])
        +(f[idx-di+dj]-2*f[idx-di]+f[idx-di-dj]));
};

__device__ FFuncType d2xCO2I2D_dev=d2xCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d2yCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(12.0*dy*dy)*(
        (f[idx+di+dj]-2*f[idx+dj]+f[idx-di+dj])
        +10*(f[idx+di]-2*f[idx]+f[idx-di])
        +(f[idx+di-dj]-2*f[idx-dj]+f[idx-di-dj]));
};

__device__ FFuncType d2yCO2I2D_dev=d2yCO2I2D;


// ------------------------------------------------------------------
// This actually is O4
CUDA_CALLABLE_MEMBER double d2x2yCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(16.0*dx*dx*dy*dy)*( f[idx+2*dj+2*di] - 2*f[idx+2*dj] + f[idx-2*di+2*dj] - 2*f[idx+2*di] + 4*f[idx] - 2*f[idx-2*di] + f[idx+2*di-2*dj] - 2*f[idx-2*dj] + f[idx-2*di-2*dj] ); 
};

__device__ FFuncType d2x2yCO2I2D_dev=d2x2yCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double laplaceCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dy)*( -10.0/3.0*f[idx]
    +2.0/3.0*( f[idx+di] + f[idx-di]
    +f[idx+dj] + f[idx-dj] )
    +1.0/6.0*( f[idx+di+dj] + f[idx-di+dj]
    + f[idx+di-dj] + f[idx-di-dj] ));
};

__device__ FFuncType laplaceCO2I2D_dev=laplaceCO2I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double biLaplaceCO2I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dx*dy*dy)*( 12*f[idx]
    -10.0/3.0*( f[idx+di] + f[idx-di]
    + f[idx+dj] + f[idx-dj] )
    -2.0/3.0*( f[idx+di+dj] + f[idx-di+dj]
    + f[idx+di-dj] + f[idx-di-dj] )
    +1.0/3.0*( f[idx+2*di] + f[idx-2*di]
    + f[idx+2*dj] + f[idx-2*dj] )
    +1.0/3.0*( f[idx+di+2*dj] + f[idx-di+2*dj]
    + f[idx+di-2*dj] + f[idx-di-2*dj]
    + f[idx+2*di+dj] + f[idx-2*di+dj]
    + f[idx+2*di-dj] + f[idx-2*di-dj] )
    );
};

__device__ FFuncType biLaplaceCO2I2D_dev=biLaplaceCO2I2D;


#endif
