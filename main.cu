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

__global__ void testGPU (FiniteDifference** FDM_ptrs, int *FDM_idx) {
    FDM_idx[0]=FDM_idx[0]+1;
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
        {-0.2,{{"1",&phia}},"explicit"},
        {1,{{"1",&phia},{"1",&phia},{"1",&phia}},"explicit"}
    });

    mub.setRhsTerms({
        {-1,{{"laplace",&phib}},"explicit"},
        {-0.2,{{"1",&phib}},"explicit"},
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
    string FDScheme="CentralDifferenceO2I";
    Evolver evolver(&mySys,0,30000,0.02,100,device,"EulerForward",FDScheme);
    // Evolver evolver(&mySys,0,10000,0.02,100,device,"RK4",FDScheme);
    evolver.run();

    // -----------------------------------------------------------
    // evolver.initEvolver();
    // evolver.evalFieldFuncs(&mua,0);
    // evolver.updateRHS(&mua,0);
    // mua.applyBounCondPeriGPU(mua.f[0]);
    // evolver.evalFieldFuncs(&phia,0);
    // evolver.updateRHS(&phia,0);
    // phia.export_conf_any(phia.f[0],"phia","0", device, 1);
    // phia.export_conf_any(phia.laplace,"phia_laplace","0", device, 1);
    // phia.export_conf_any(phia.rhs[0],"phia_rhs","0", device, 1);
    // mua.export_conf_any(mua.f[0],"mua","0", device, 1);
    // mua.export_conf_any(mua.laplace,"mua_laplace","0", device, 1);        

    // -----------------------------------------------------------
    
    return 0;
};

