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
#include "userDefinedFunction.h"
#include "../src/evolver/evolverclass.cu"

using namespace std;


// ======================================================================
int main() {

    // Simulation parameters
    string direExpo="data/";
    string device="gpu";
    string FDMScheme="CentralDifferenceO4Iso2D";
    string timeScheme="RK4";
    double dt=0.001;
    double T=10;
    double dtExpo=1;
    int    NGrid=128;
    double L=NGrid*1.0;       
    
    // Generating a new system.
    System mySys(direExpo);
    // Generating a new mesh.
    Mesh mesh(2);
    mesh.setGridNumber(NGrid,NGrid,1);
    mesh.setBoxSize(L,L,1);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    LivingLC livingLC(&mesh, "livingLC", 0);    
    Field phi(&mesh, "phi",1);
    livingLC.Pxx.initFieldSin(0.5, 2, 0);
    livingLC.Pxy.initFieldSin(0.3, 4, 0);
    
    livingLC.Pxx.setRhsTerms({
        {-1,{{&livingLC.Pxx}}}
    });

    livingLC.Pxy.setRhsTerms({
        {-1,{{&livingLC.Pxy}}}
    });
    
    phi.setRhsTerms({
        {{{"sin",&livingLC.Pxx,{2,0}}}}
    });
    
    // phi.initFieldGaus(L/2, 0.1*L, 1.5);
    mySys.addLivingLC(&livingLC);
    mySys.addField(&phi);
    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:    
    Evolver evolver(&mySys,0,T,dt,dtExpo,device,timeScheme,FDMScheme);
    // evolver.initEvolver();
    // Running simulations
    evolver.run();

    // f.export_conf_any(f.rhs[0], "f_rhs", "", device);
    // phi.export_conf_any(phi.rhs[0], "phi_rhs", "", device);
    // psi.export_conf_any(psi.rhs[0], "psi_rhs", "", device);
    
    // -----------------------------------------------------------
    
    return 0;
};

