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
    double gammav=0.1;
    double kp=1;
    // Phase field
    double GammaPhi=0.01;
    double DPhi=0.01;
    double aPhi=DPhi;
    double phi0=1;
    double gPhi=0.01;
    
    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    mesh.setGridNumber(NGrid,NGrid,1);
    mesh.setBoxSize(L,L,1);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    ConstantField ConstF(&mesh,"ConstF");
    Field S2(&mesh, "S2",1);
    Field sigxx(&mesh, "sigxx",1);
    Field sigxy(&mesh, "sigxy",1);
    Field sigyx(&mesh, "sigyx",1);
    Field sigyy(&mesh, "sigyy",1);
    Field vx(&mesh, "vx",1);
    Field vy(&mesh, "vy",1);
    Field phiCorr(&mesh, "phiCorr",1);
    Field phivx(&mesh, "phivx",1);
    Field phivy(&mesh, "phivy",1);
    
    Field Qxx(&mesh, "Qxx",0,"zeros");
    Field Qxy(&mesh, "Qxy",0,"zeros","periodic","on");    
    Field phi(&mesh, "phi",0,"zeros","periodic","on");
    
    
    S2.setRhsTerms({
        {4,{{&Qxx},{&Qxx}}},
        {4,{{&Qxy},{&Qxy}}}
    });
    // sigxx.setRhsTerms({
        // {kp,{{&phi}}}, {-kp*phi0,{{&ConstF.one}}}
    // });
    sigxx.setRhsTerms({
        {{{"myFunc",&phi}}}
    });
    sigyy.setRhsTerms({
        {{{"myFunc",&phi}}}
    });
    sigxy.setRhsTerms({
        {{{&ConstF.zero}}}
    });
    sigyx.setRhsTerms({
        {{{&ConstF.zero}}}
    });
    // sigyy.setRhsTerms({
        // {kp,{{&phi}}}, {-kp*phi0,{{&ConstF.one}}}
    // });
    
    phivx.setRhsTerms({
        {-1/gammav, {{"d1x",&sigxx},{&phi}}},
        {-1/gammav, {{"d1y",&sigxy},{&phi}}}
    });
    phivy.setRhsTerms({
        {-1/gammav, {{"d1x",&sigyx},{&phi}}},
        {-1/gammav, {{"d1y",&sigyy},{&phi}}}
    });
    
    
    // Qxx.setRhsTerms({
    //     {-1, { {&incomFlow.vx}, {"d1x",&Qxx} }},
    //     {-1,{ {&incomFlow.vy}, {"d1y",&Qxx} }},
    //     {lambda, { {"d1x",&incomFlow.vx}, {&phi} }},
    //     {-1, { {"d1x",&incomFlow.vy}, {&Qxy}, {&phi} }},
    //     {{ {"d1y",&incomFlow.vx},{&Qxy} }},        
    //     {-GammaQ*aQ,{{&Qxx}, {&phi}}},
    //     {-0.5*GammaQ*bQ, {{&S2},{&Qxx}}},
    //     {GammaQ*KQ,{{"laplace",&Qxx}}}
    // });
    
    // Qxy.setRhsTerms({
    //     {-1, { {&incomFlow.vx}, {"d1x",&Qxy} }},
    //     {-1,{ {&incomFlow.vy}, {"d1y",&Qxy} }},
    //     {0.5*lambda, { {"d1x",&incomFlow.vy}, {&phi} }},
    //     {0.5*lambda, { {"d1y",&incomFlow.vx}, {&phi} }},
    //     {{ {"d1x",&incomFlow.vy},{&Qxx} }},
    //     {-1, { {"d1y",&incomFlow.vx},{&Qxx} }},        
    //     {-GammaQ*aQ,{{&Qxy}, {&phi} }},
    //     {-0.5*GammaQ*bQ, {{&S2},{&Qxy}}},
    //     {GammaQ*KQ,{{"laplace",&Qxy}}}
    // });
    
    phi.setRhsTerms({
        {-1, { {"d1x",&phivx} }},
        {-1, { {"d1y",&phivy} }},
        {DPhi,{{"laplace",&phi}}},
        {gPhi, {{&phi}}}
    });    
    
    // Add fields to the system.
    // mySys.field_ptrs.push_back(&S2);
    // mySys.field_ptrs.push_back(&mu);
    // mySys.field_ptrs.push_back(&Qxx);
    // mySys.field_ptrs.push_back(&Qxy);
    mySys.addConstantField(&ConstF);
    mySys.field_ptrs.push_back(&sigxx);
    mySys.field_ptrs.push_back(&sigxy);
    mySys.field_ptrs.push_back(&sigyx);
    mySys.field_ptrs.push_back(&sigyy);    
    mySys.field_ptrs.push_back(&phi);
    mySys.field_ptrs.push_back(&phivx);
    mySys.field_ptrs.push_back(&phivy);
    
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

