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
    Field mu1(&mesh, "mu1", 1);
    Field phi1(&mesh, "phi1");

    mu.setRhsTerms({
        {-1,{{&phi}}},
        {1,{{&phi},{&phi},{&phi}}},
        {-1,{{"laplace",&phi}}},
    });
    
    phi.setRhsTerms({
        {{{"laplace",&mu}}},
        {-0.2, {{"laplace",&phi1}}}
    });

    mu1.setRhsTerms({
        {-1,{{&phi1}}},
        {1,{{&phi1},{&phi1},{&phi1}}},
        {-1,{{"laplace",&phi1}}},
    });
    
    phi1.setRhsTerms({
        {{{"laplace",&mu1}}},
        {0.2, {{"laplace",&phi}}}
    });
    
    // Add fields to the system.
    mySys.field_ptrs.push_back(&mu);
    mySys.field_ptrs.push_back(&phi);
    mySys.field_ptrs.push_back(&mu1);
    mySys.field_ptrs.push_back(&phi1);
    
    // Print system information.
    // mySys.printSysInfo();    
    // Creating an evolver:
    string device="gpu";
    string FDMScheme="CentralDifferenceO2Iso2D";
    // Low activity
    Evolver evolver(&mySys,0,10000,0.01,50,device,"EulerForward",FDMScheme);
    evolver.run();    
    // Evolver evolver(&mySys,0,5000,0.01,50,device,"EulerForward",FDMScheme);

    // Running simulations
    

    
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};

