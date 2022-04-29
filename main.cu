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


// ======================================================================
int main() {
    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field mua(&mesh, "mua", 0, 1, "periodic", "Gaussian", "on");
    Field phia(&mesh, "phia", 0, 0, "periodic", "sin", "on");
    Field mub(&mesh, "mub", 0, 1, "periodic", "Gaussian", "on");
    Field phib(&mesh, "phib", 0, 0, "periodic", "sin", "on");

    // Set field equations.
    mua.setRhsTerms({
        {-1,{{"laplace",&phia}},"explicit"},
        {-0.05,{{"1",&phia}},"explicit"},
        {1,{{"1",&phia},{"1",&phia},{"1",&phia}},"explicit"}
    });

    mub.setRhsTerms({
        {-1,{{"laplace",&phib}},"explicit"},
        {-0.05,{{"1",&phib}},"explicit"},
        {1,{{"1",&phib},{"1",&phib},{"1",&phib}},"explicit"}
    });
    
    phia.setRhsTerms({
        {1,{{"laplace",&mua}},"explicit"},
        {-0.1,{{"laplace",&phib}},"explicit"}
    });

    phib.setRhsTerms({
        {1,{{"laplace",&mub}},"explicit"},
        {0.1,{{"laplace",&phia}},"explicit"}
    });
    
    // Add fields to the system.
    mySys.field_ptrs.push_back(&mua);
    mySys.field_ptrs.push_back(&mub);
    mySys.field_ptrs.push_back(&phia);
    mySys.field_ptrs.push_back(&phib);
    // Print system information.
    // mySys.printSysInfo();
    
    // Creating an evolver:
    string device="gpu";
    Evolver evolver(&mySys,0,30000,0.02,100,"EulerForward",device);    
    evolver.run();
    
    
    return 0;
};

