#ifndef FINITEDIFFERENCE_H
#define FINITEDIFFERENCE_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif


struct FDMScheme {
    std::string scheme_name;
};
FDMScheme central_O2;
FDMScheme central_O2I;
FDMScheme central_O4;
FDMScheme central_O4I;

template <FDMScheme*> struct FDM;

#endif
