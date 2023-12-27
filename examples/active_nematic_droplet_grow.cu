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
#include "/home/you/Research/codes/PDESolver/src/evolver/evolverclass.cu"

using namespace std;


// ======================================================================
int main() {

    // Simulation parameters
    string direExpo="data/";
    string device="gpu";
    string FDMScheme="CentralDifferenceO4Iso2D";
    string timeScheme="RK4";
    double dt=0.001;
    double T=1000;
    double dtExpo=1;
    int    NGrid=128;
    double L=NGrid*1.0;
    
    // Model parameters
    // Nematics
    double GammaQ=1;
    double KQ=1;
    double aQ=-1;
    double bQ=2.5;
    double lambda=0.1;
    double alpha=-20;
    // Flow
    double eta=10;
    double Gamma=0.1;
    // Phase field
    double GammaPhi=0.01;
    double KPhi=1000;
    double aPhi=KPhi;
    double phi0=0.3;
    double gPhi=0.01;
    
    // Generating a new system.
    System mySys(direExpo);
    // Generating a new mesh.
    Mesh mesh(2);
    mesh.setGridNumber(NGrid,NGrid,1);
    mesh.setBoxSize(L,L,1);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field S2(&mesh, "S2",1);
    Field mu(&mesh, "mu",1);
    Field Qxx(&mesh, "Qxx",0,"zeros");
    Field Qxy(&mesh, "Qxy",0,"zeros","periodic","on");    
    Field phi(&mesh, "phi",0,"zeros","periodic","on");
    IncompressibleFlow incomFlow(&mesh, "incomFlow");   
    
    S2.setRhsTerms({
        {4,{{&Qxx},{&Qxx}}},
        {4,{{&Qxy},{&Qxy}}}
    });

    mu.setRhsTerms({
        {2*aPhi,{{&phi}}},
        {-6*aPhi,{{&phi},{&phi}}},
        {4*aPhi,{{&phi},{&phi},{&phi}}},
        {-KPhi,{{"laplace",&phi}}},
    });
    
    Qxx.setRhsTerms({
        {-1, { {&incomFlow.vx}, {"d1x",&Qxx} }},
        {-1,{ {&incomFlow.vy}, {"d1y",&Qxx} }},
        {lambda, { {"d1x",&incomFlow.vx}, {&phi} }},
        {-1, { {"d1x",&incomFlow.vy}, {&Qxy}, {&phi} }},
        {{ {"d1y",&incomFlow.vx},{&Qxy} }},        
        {-GammaQ*aQ,{{&Qxx}, {&phi}}},
        {-0.5*GammaQ*bQ, {{&S2},{&Qxx}}},
        {GammaQ*KQ,{{"laplace",&Qxx}}}
    });

    Qxy.setRhsTerms({
        {-1, { {&incomFlow.vx}, {"d1x",&Qxy} }},
        {-1,{ {&incomFlow.vy}, {"d1y",&Qxy} }},
        {0.5*lambda, { {"d1x",&incomFlow.vy}, {&phi} }},
        {0.5*lambda, { {"d1y",&incomFlow.vx}, {&phi} }},
        {{ {"d1x",&incomFlow.vy},{&Qxx} }},
        {-1, { {"d1y",&incomFlow.vx},{&Qxx} }},        
        {-GammaQ*aQ,{{&Qxy}, {&phi} }},
        {-0.5*GammaQ*bQ, {{&S2},{&Qxy}}},
        {GammaQ*KQ,{{"laplace",&Qxy}}}
    });

    phi.setRhsTerms({
        {-1, { {&incomFlow.vx}, {"d1x",&phi} }},
        {-1,{ {&incomFlow.vy}, {"d1y",&phi} }},
        {GammaPhi,{{"laplace",&mu}}},
        {gPhi, {{&phi}}}
    });
    
    incomFlow.omega.setRhsTerms({
        {-0, { {&incomFlow.vx}, {"d1x",&incomFlow.omega} }},
        {-0, { {&incomFlow.vy}, {"d1y",&incomFlow.omega} }},
        {eta, {{"laplace",&incomFlow.omega}}},
        {-Gamma, {{&incomFlow.omega}}},
        {alpha, {{"d2x",&Qxy}, {&phi}}},
        {-alpha, {{"d2y",&Qxy}, {&phi}}},
        {-2*alpha, {{"d1x1y",&Qxx}, {&phi}}},
        {-1, {{"d1x",&phi}, {"d1y",&mu}}},
        { 1, {{"d1y",&phi}, {"d1x",&mu}}}
    });
    
    // Add fields to the system.
    mySys.addField(&S2);
    mySys.addField(&mu);
    mySys.addField(&Qxx);
    mySys.addField(&Qxy);
    mySys.addField(&phi);
    mySys.addIncompressibleFlow(&incomFlow);
    // phi.initFieldConst(phi0);
    phi.initFieldGaus(L/2, 0.1*L, 1.5);
    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:    
    Evolver evolver(&mySys,0,T,dt,dtExpo,device,timeScheme,FDMScheme);

    // Running simulations
    evolver.run();
    
    
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};

