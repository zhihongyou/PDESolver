#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "cuda.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
#include "../src/evolver/evolverclass.cu"

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
    Field mu(&mesh, "mu", 1);
    Field phi(&mesh, "phi");

    mu.setRhsTerms({
        {-1,{{&phi}}},
        {1,{{&phi},{&phi},{&phi}}},
        {-1,{{"laplace",&phi}}},
    });
    
    phi.setRhsTerms({
        {{{"laplace",&mu}}}
    });
    
    // Add fields to the system.
    mySys.field_ptrs.push_back(&mu);
    mySys.field_ptrs.push_back(&phi);
    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:
    string device="gpu";
    string FDMScheme="CentralDifferenceO2Iso2D";
    // Low activity
    Evolver evolver(&mySys,0,500,0.01,5,device,"EulerForward",FDMScheme);    
    // Evolver evolver(&mySys,0,5000,0.01,50,device,"EulerForward",FDMScheme);

    // Running simulations
    evolver.run();    

    
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};

