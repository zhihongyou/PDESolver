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



// ======================================================================
int main() {
    
    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field phi(&mesh, "phi", 0, 0, "periodic", "Gaussian", "on");
    // Field phi(&mesh, "phi", 0, 0, "periodic", "sin", "on");
    // Set field equations.
    phi.setRhsTerms({
        {2.5,{{"laplace",&phi}},"explicit"}
    });

    // Add fields to the system.
    mySys.field_ptrs.push_back(&phi);
    // Print system information.
    // mySys.printSysInfo();
    
    // Creating an evolver:
    string device="gpu";
    Evolver evolver(&mySys,0,100,0.01,1,"EulerForward",device);
    // evolver.EulerForward();
    
    // Run simulations.
    evolver.run();
    
    return 0;
};
