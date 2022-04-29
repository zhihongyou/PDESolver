#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "cuda_runtime.h"
#include "cuda.h"
#include "device_launch_parameters.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
// #include "src/mesh/meshclass.cuh"
// #include "src/field/fieldclass.cu"
// #include "src/system/systemclass.cpp"
#include "src/evolver/evolverclass.h"
#include "src/evolver/evolverclass.cu"
// #include "src/utility/finiteDifference.h"
// #include "src/utility/finiteDifferenceCentralO2Isotropic.h"

using namespace std;

__global__ void testGPU (rhsPtrs rhs_ptrs_dev, int a) {
    rhs_ptrs_dev.num_terms[0]+=a;
    // rhs_ptrs_dev.
    // f_dev[0][10]=1000;
};

__global__ void testGPU1 (double* prefactors) {
    // rhs_ptrs_dev.num_terms[0]+=0;
    prefactors[0]+=11;
    // rhs_ptrs_dev.num_funcs_1term[0]+=12;
};


// ======================================================================
int main() {
    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field mu(&mesh, "mu", 0, 1, "periodic", "Gaussian", "on");
    Field phi(&mesh, "phi", 0, 0, "periodic", "sin", "on");

    // Set field equations.
    mu.setRhsTerms({
        {-1,{{"laplace",&phi}},"explicit"},
        {-0.2,{{"1",&phi}},"explicit"},
        {1,{{"1",&phi},{"1",&phi},{"1",&phi}},"explicit"}
    });
    
    phi.setRhsTerms({
        {2.5,{{"laplace",&mu}},"explicit"}
    });
    
    // Add fields to the system.
    mySys.field_ptrs.push_back(&mu);
    mySys.field_ptrs.push_back(&phi);
    // Print system information.
    // mySys.printSysInfo();
    
    // Creating an evolver:
    string device="gpu";
    Evolver evolver(&mySys,0,10000,0.01,100,"EulerForward",device);    
    evolver.run();
    
    
    return 0;
};

