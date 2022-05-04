#ifndef FINITEDIFFERENCECENTRALO4ISOTROPIC_H
#define FINITEDIFFERENCECENTRALO4ISOTROPIC_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "vector3class.h"
#include "finiteDifference.h"

class FDMCentralO4I: public FiniteDifference {
    
public:

    string name() {
        string name_t="CentralDifferenceO4I";
        return name_t;
    };
    
    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1x(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/dx*(
		  13.0/30.0*( f[idx+dj] - f[idx-dj] )
		  +2.0/15.0*( f[idx+di+dj] + f[idx-di+dj] - f[idx+di-dj] - f[idx-di-dj] )
		  -1.0/60.0*( f[idx+2*dj] - f[idx-2*dj] )
		  -1.0/60.0*( f[idx+2*di+dj] + f[idx-2*di+dj] - f[idx+2*di-dj] - f[idx-2*di-dj] )
		  -1.0/30.0*( f[idx+di+2*dj] + f[idx-di+2*dj] - f[idx+di-2*dj] - f[idx-di-2*dj] )
		  );
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/dy*(
		  13.0/30.0*( f[idx+dj] - f[idx-dj] )
		  +2.0/15.0*( f[idx+di+dj] + f[idx-di+dj] - f[idx+di-dj] - f[idx-di-dj] )
		  -1.0/60.0*( f[idx+2*dj] - f[idx-2*dj] )
		  -1.0/60.0*( f[idx+2*di+dj] + f[idx-2*di+dj] - f[idx+2*di-dj] - f[idx-2*di-dj] )
		  -1.0/30.0*( f[idx+di+2*dj] + f[idx-di+2*dj] - f[idx+di-2*dj] - f[idx-di-2*dj] )
		  );
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d1x1y(double * f, int idx, int di, int dj, double dx, double dy) {        
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

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d2x(double * f, int idx, int di, int dj, double dx, double dy) {        
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

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double d2y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(dy*dy)*( 5.0*f[idx]
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

    // ------------------------------------------------------------------
    // This actually is O4
    CUDA_CALLABLE_MEMBER double d2x2y(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(16.0*dx*dx*dy*dy)*( f[idx+2*dj+2*di] - 2*f[idx+2*dj] + f[idx-2*di+2*dj] - 2*f[idx+2*di] + 4*f[idx] - 2*f[idx-2*di] + f[idx+2*di-2*dj] - 2*f[idx-2*dj] + f[idx-2*di-2*dj] ); 
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double laplace(double * f, int idx, int di, int dj, double dx, double dy) {        
        return 1.0/(dx*dy)*( -21.0/5.0*f[idx]
			  +13.0/15.0*( f[idx+di] + f[idx-di] + f[idx+dj] + f[idx-dj] )
			  +4.0/15.0*( f[idx+di+dj] + f[idx-di+dj] + f[idx+di-dj] + f[idx-di-dj] )
			  -1.0/60.0*( f[idx+2*di] + f[idx-2*di] + f[idx+2*dj] + f[idx-2*dj] )
			  -1.0/30.0*( f[idx+di+2*dj] + f[idx-di+2*dj] + f[idx+di-2*dj] + f[idx-di-2*dj]
				     + f[idx+2*di+dj] + f[idx-2*di+dj] + f[idx+2*di-dj] + f[idx-2*di-dj] )
			  );
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER double bi_laplace(double * f, int idx, int di, int dj, double dx, double dy) {        
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
        
    
};


#endif
