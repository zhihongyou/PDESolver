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
    Field f(&mesh, "f",0);
    // Field phi(&mesh, "phi",1);
    LaplaceNFEqField phi(&mesh, "phi");
    double prefactors[2]={2,10};
    phi.setLaplaceNFEq(1,prefactors);
    cout <<phi.LaplNFEqSolver.max_power<<endl;
    cout <<phi.LaplNFEqSolver.prefactors[0]<<", "<<phi.LaplNFEqSolver.prefactors[1]<<endl;
    cout <<phi.specialty<<endl;
    cout <<phi.traits_host.priority<<endl;
    // LaplaceNFEqSolver LaplNFEqSolver;
    // LaplNFEqSolver.initLaplaceNFEqSolver(&mesh);
    
    f.setRhsTerms({
        {{{&f}}}
    });
    
    phi.setRhsTerms({
        {{{&f}}}
    });
    
    f.initFieldGaus(L/2, 0.1*L, 1.5);
    
    // phi.initFieldGaus(L/2, 0.1*L, 1.5);
    mySys.addField(&f);
    mySys.addField(&phi);
    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:    
    Evolver evolver(&mySys,0,T,dt,dtExpo,device,timeScheme,FDMScheme);
    // evolver.initEvolver();
    // Running simulations
    evolver.run();

    phi.export_conf("1000", device);
    f.export_conf("1000", device);
    phi.LaplNFEqSolver.solveLaplaceNFEq(phi.f[0],f.f[0]);
    phi.export_conf("1001", device);
    f.export_conf("1001", device);
    // phi.export_conf_any(phi.phi_complex[3], "phi_rhs", "3", "gpu", 1);
    // evolver.getRHS(0);
    // f.export_conf("0","gpu");
    // phi.export_conf("0","gpu");
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};

