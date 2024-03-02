#ifndef FINITEDIFFERENCECENTRALO4ISOTROPIC2D_H
#define FINITEDIFFERENCECENTRALO4ISOTROPIC2D_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif


/* The name "CO4I2D" stands for "centeral difference,
    4th order, Isotropic, 2D".*/
    
// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1xCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/dx*(
        13.0/30.0*( f[idx+dj] - f[idx-dj] )
        +2.0/15.0*( f[idx+di+dj] + f[idx-di+dj] - f[idx+di-dj] - f[idx-di-dj] )
        -1.0/60.0*( f[idx+2*dj] - f[idx-2*dj] )
        -1.0/60.0*( f[idx+2*di+dj] + f[idx-2*di+dj] - f[idx+2*di-dj] - f[idx-2*di-dj] )
        -1.0/30.0*( f[idx+di+2*dj] + f[idx-di+2*dj] - f[idx+di-2*dj] - f[idx-di-2*dj] )
    );
};

__device__ FFuncType d1xCO4I2D_dev=d1xCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1yCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/dy*(
        13.0/30.0*( f[idx+di] - f[idx-di] )
        +2.0/15.0*( f[idx+dj+di] + f[idx-dj+di] - f[idx+dj-di] - f[idx-dj-di] )
        -1.0/60.0*( f[idx+2*di] - f[idx-2*di] )
        -1.0/60.0*( f[idx+2*dj+di] + f[idx-2*dj+di] - f[idx+2*dj-di] - f[idx-2*dj-di] )
        -1.0/30.0*( f[idx+dj+2*di] + f[idx-dj+2*di] - f[idx+dj-2*di] - f[idx-dj-2*di] )
    );
};

__device__ FFuncType d1yCO4I2D_dev=d1yCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d1x1yCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dy)*(17.0/45.0*(f[idx-di-dj]+f[idx+di+dj]
    -f[idx+di-dj]-f[idx-di+dj])
    -1.0/45.0*(f[idx+2*di+dj]+f[idx+di+2*dj]
    +f[idx-2*di-dj]+f[idx-di-2*dj]
    -f[idx-di+2*dj]-f[idx-2*di+dj]
    -f[idx+2*di-dj]-f[idx+di-2*dj])
    -7.0/720.0*(f[idx-2*di-2*dj]+f[idx+2*di+2*dj]
    -f[idx-2*di+2*dj]-f[idx+2*di-2*dj])
    );
};

__device__ FFuncType d1x1yCO4I2D_dev=d1x1yCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d2xCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dx)*( 5.0*f[idx]
    -164.0/45.0*( f[idx+dj] + f[idx-dj] )
    +103.0/90.0*( f[idx+2*dj] + f[idx-2*dj] )
    -223.0/45.0*( f[idx+di] + f[idx-di])  
    +148.0/45.0*( f[idx+di+dj] + f[idx-di+dj] +
    f[idx+di-dj] + f[idx-di-dj] )
    -73.0/90.0* ( f[idx+di+2*dj] + f[idx-di+2*dj] +
    f[idx+di-2*dj] + f[idx-di-2*dj] )
    +217.0/180.0*(f[idx+2*di] + f[idx-2*di])
    -4.0/5.0*(  + f[idx+2*di+dj] + f[idx-2*di+dj] +
    f[idx+2*di-dj] + f[idx-2*di-dj])
    +71.0/360.0*(+f[idx+2*di+2*dj] + f[idx-2*di+2*dj] +
    f[idx+2*di-2*dj] + f[idx-2*di-2*dj])
    );
};

__device__ FFuncType d2xCO4I2D_dev=d2xCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double d2yCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dy*dy)*( 5.0*f[idx]
    -164.0/45.0*( f[idx+di] + f[idx-di] )
    +103.0/90.0*( f[idx+2*di] + f[idx-2*di] )
    -223.0/45.0*( f[idx+dj] + f[idx-dj])  
    +148.0/45.0*( f[idx+dj+di] + f[idx-dj+di] +
    f[idx+dj-di] + f[idx-dj-di] )
    -73.0/90.0* ( f[idx+dj+2*di] + f[idx-dj+2*di] +
    f[idx+dj-2*di] + f[idx-dj-2*di] )
    +217.0/180.0*(f[idx+2*dj] + f[idx-2*dj])
    -4.0/5.0*(  + f[idx+2*dj+di] + f[idx-2*dj+di] +
    f[idx+2*dj-di] + f[idx-2*dj-di])
    +71.0/360.0*(+f[idx+2*dj+2*di] + f[idx-2*dj+2*di] +
    f[idx+2*dj-2*di] + f[idx-2*dj-2*di])
    );
};

__device__ FFuncType d2yCO4I2D_dev=d2yCO4I2D;


// ------------------------------------------------------------------
// This actually is O4
CUDA_CALLABLE_MEMBER double d2x2yCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(16.0*dx*dx*dy*dy)*( f[idx+2*dj+2*di] - 2*f[idx+2*dj] + f[idx-2*di+2*dj] - 2*f[idx+2*di] + 4*f[idx] - 2*f[idx-2*di] + f[idx+2*di-2*dj] - 2*f[idx-2*dj] + f[idx-2*di-2*dj] ); 
};

__device__ FFuncType d2x2yCO4I2D_dev=d2x2yCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double laplaceCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dy)*( -21.0/5.0*f[idx]
    +13.0/15.0*( f[idx+di] + f[idx-di] + f[idx+dj] + f[idx-dj] )
    +4.0/15.0*( f[idx+di+dj] + f[idx-di+dj] + f[idx+di-dj] + f[idx-di-dj] )
    -1.0/60.0*( f[idx+2*di] + f[idx-2*di] + f[idx+2*dj] + f[idx-2*dj] )
    -1.0/30.0*( f[idx+di+2*dj] + f[idx-di+2*dj] + f[idx+di-2*dj] + f[idx-di-2*dj]
    + f[idx+2*di+dj] + f[idx-2*di+dj] + f[idx+2*di-dj] + f[idx-2*di-dj] )
    );
};

__device__ FFuncType laplaceCO4I2D_dev=laplaceCO4I2D;


// ------------------------------------------------------------------
CUDA_CALLABLE_MEMBER double biLaplaceCO4I2D(double * f, int idx, FFuncArgs f_func_args) {
    double dx=f_func_args.dx;
    double dy=f_func_args.dy;
    int di=f_func_args.di;
    int dj=f_func_args.dj;
    return 1.0/(dx*dx*dy*dy)*( 779.0/45.0*f[idx]
    -191.0/45.0*( f[idx+di] + f[idx-di] + f[idx+dj] + f[idx-dj] )
    -187.0/90.0*( f[idx+di+dj] + f[idx-di+dj] + f[idx+di-dj] + f[idx-di-dj] )
    +7.0/30.0*( f[idx+2*di] + f[idx-2*di] + f[idx+2*dj] + f[idx-2*dj] )
    +47.0/45.0*( f[idx+di+2*dj] + f[idx-di+2*dj] + f[idx+di-2*dj] + f[idx-di-2*dj]
    + f[idx+2*di+dj] + f[idx-2*di+dj] + f[idx+2*di-dj] + f[idx-2*di-dj] )
    -29.0/180.0*( f[idx+2*di+2*dj] + f[idx-2*di+2*dj] + f[idx+2*di-2*dj] + f[idx-2*di-2*dj] )
    +1.0/45.0*( f[idx+3*di] + f[idx-3*di] + f[idx+3*dj] + f[idx-3*dj] )
    -17.0/180.0*( f[idx+di+3*dj] + f[idx-di+3*dj] + f[idx+di-3*dj] + f[idx-di-3*dj]
    + f[idx+3*di+dj] + f[idx-3*di+dj] + f[idx+3*di-dj] + f[idx-3*di-dj] )
    );
};
        
__device__ FFuncType biLaplaceCO4I2D_dev=biLaplaceCO4I2D;


#endif
