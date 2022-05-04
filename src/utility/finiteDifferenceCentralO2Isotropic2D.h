#ifndef FINITEDIFFERENCECENTRALO2ISOTROPIC2D_H
#define FINITEDIFFERENCECENTRALO2ISOTROPIC2D_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "vector3class.h"
#include "finiteDifference.h"

class FDMCentralO2Iso2D: public FiniteDifference {
    
public:

    string name() {
        string name_t="CentralDifferenceO2Iso2D";
        return name_t;
    };
    
    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1x(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(12.0*dx)*(
            (f[idx+di+dj] - f[idx+di-dj])
            +4*(f[idx+dj] - f[idx-dj])
            +(f[idx-di+dj] - f[idx-di-dj]));
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(12.0*dy)*(
            (f[idx+di+dj] - f[idx-di+dj])
            +4*(f[idx+di] - f[idx-di])
            +(f[idx+di-dj] - f[idx-di-dj]));
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1x1y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(4.0*dx*dy)*(f[idx+di+dj]-f[idx-di+dj]
        -f[idx+di-dj]+f[idx-di-dj]);
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d2x(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(12.0*dx*dx)*(
            (f[idx+di+dj]-2*f[idx+di]+f[idx+di-dj])
            +10*(f[idx+dj]-2*f[idx]+f[idx-dj])
            +(f[idx-di+dj]-2*f[idx-di]+f[idx-di-dj]));
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d2y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(12.0*dy*dy)*(
            (f[idx+di+dj]-2*f[idx+dj]+f[idx-di+dj])
            +10*(f[idx+di]-2*f[idx]+f[idx-di])
            +(f[idx+di-dj]-2*f[idx-dj]+f[idx-di-dj]));
    };

    // ------------------------------------------------------------------
    // This actually is O4
    CUDA_CALLABLE_MEMBER double d2x2y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(16.0*dx*dx*dy*dy)*( f[idx+2*dj+2*di] - 2*f[idx+2*dj] + f[idx-2*di+2*dj] - 2*f[idx+2*di] + 4*f[idx] - 2*f[idx-2*di] + f[idx+2*di-2*dj] - 2*f[idx-2*dj] + f[idx-2*di-2*dj] ); 
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double laplace(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(dx*dy)*( -10.0/3.0*f[idx]
        +2.0/3.0*( f[idx+di] + f[idx-di]
        +f[idx+dj] + f[idx-dj] )
        +1.0/6.0*( f[idx+di+dj] + f[idx-di+dj]
        + f[idx+di-dj] + f[idx-di-dj] ));
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double bi_laplace(double * f, int idx, int di, int dj, double dx, double dy) {        
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
        
    
};


#endif
