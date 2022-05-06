#ifndef FIELDFUNCTIONCLASS_H
#define FIELDFUNCTIONCLASS_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif


class FieldFunction {
    
public:

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER static double oneF(double * f, int idx) {
        return f[idx];
    };
    
    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER static double oneOverF(double * f, int idx) {
        // return 1.0/f[idx];
        return 4*f[idx];
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER static double sinF(double * f, int idx) {
        return sin(f[idx]);
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER static double cosF(double * f, int idx) {
        return cos(f[idx]);
    };

    // ------------------------------------------------------------------
    CUDA_CALLABLE_MEMBER static double absF(double * f, int idx) {
        return abs(f[idx]);
    };
};


#endif
