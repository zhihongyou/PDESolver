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
#include "src/evolver/evolverclass.cu"

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
    Field mua(&mesh, "mua", 0, 1, "periodic", "Gaussian", "on");
    Field phia(&mesh, "phia", 0, 0, "periodic", "sin", "on");
    Field mub(&mesh, "mub", 0, 1, "periodic", "Gaussian", "on");
    Field phib(&mesh, "phib", 0, 0, "periodic", "sin", "on");

    // Set equations
    mua.setRhsTerms({
        {-1, {{"laplace",&phia}}},
        {-0.2, {{&phia}}},
        {{ {&phia}, {&phia}, {&phia} }}
    });

    mub.setRhsTerms({
        {-1,{{"laplace",&phib}}},
        {-0.2,{{"1",&phib}}},
        {{ {&phib}, {&phib}, {&phib} }}
    });
    
    phia.setRhsTerms({
        {{{"laplace",&mua}}},
        {-0.1, {{"laplace",&phib}}}
    });

    phib.setRhsTerms({
        {{{"laplace",&mub}}},
        {0.1, {{"laplace",&phia}}}
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
    string FDScheme="CentralDifferenceO2Iso2D";
    Evolver evolver(&mySys,0,20000,0.02,100,device,"EulerForward",FDScheme);

    // Running simulations
    evolver.run();    

    // -----------------------------------------------------------
    // Testing

    // -----------------------------------------------------------
    
    return 0;
};

