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
    Field S2(&mesh, "S2",1);
    Field Qxx(&mesh, "Qxx",0,"sin");
    Field Qxy(&mesh, "Qxy",0,"sin","periodic","on");
    IncompressibleFlow incomFlow(&mesh, "incomFlow");

    // Set equations
    double gamma=1;
    double K=1;
    double A2=-0.25;
    double A4=1;
    double lambda=0.7;
    double eta=10;
    double Gamma=0.1;
    double alpha=-40;
    // double alpha=-40;
    
    S2.setRhsTerms({
        {4,{{&Qxx},{&Qxx}}},
        {4,{{&Qxy},{&Qxy}}}
    });
    
    Qxx.setRhsTerms({
        {-1, { {&incomFlow.vx}, {"d1x",&Qxx} }}, {-1,{ {&incomFlow.vy}, {"d1y",&Qxx} }},
        {lambda, { {"d1x",&incomFlow.vx} }},
        {-1, { {"d1x",&incomFlow.vy},{&Qxy} }}, {{ {"d1y",&incomFlow.vx},{&Qxy} }},
        {K/gamma,{{"laplace",&Qxx}}}, {-A2/gamma,{{&Qxx}}}, {-0.5*A4/gamma, {{&S2},{&Qxx}}}
    });

    Qxy.setRhsTerms({
        {-1, { {&incomFlow.vx}, {"d1x",&Qxy} }}, {-1,{ {&incomFlow.vy}, {"d1y",&Qxy} }},
        {0.5*lambda, { {"d1x",&incomFlow.vy} }}, {0.5*lambda, { {"d1y",&incomFlow.vx} }},
        {{ {"d1x",&incomFlow.vy},{&Qxx} }}, {-1, { {"d1y",&incomFlow.vx},{&Qxx} }},
        {K/gamma,{{"laplace",&Qxy}}}, {-A2/gamma,{{&Qxy}}}, {-0.5*A4/gamma, {{&S2},{&Qxy}}}
    });
    
    incomFlow.omega.setRhsTerms({
        {-0, { {&incomFlow.vx}, {"d1x",&incomFlow.omega} }},
        {-0, { {&incomFlow.vy}, {"d1y",&incomFlow.omega} }},
        {eta, {{"laplace",&incomFlow.omega}}}, {-Gamma, {{&incomFlow.omega}}},
        {alpha, {{"d2x",&Qxy}}}, {-alpha, {{"d2y",&Qxy}}}, {-2*alpha, {{"d1x1y",&Qxx}}}
    });
    
    // Add fields to the system.
    mySys.field_ptrs.push_back(&S2);
    mySys.field_ptrs.push_back(&Qxx);
    mySys.field_ptrs.push_back(&Qxy);
    mySys.addIncompressibleFlow(&incomFlow);
    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:
    string device="gpu";
    string FDMScheme="CentralDifferenceO2Iso2D";
    // Low activity
    // Evolver evolver(&mySys,0,500,0.005,2,device,"EulerForward",FDMScheme);
    // Evolver evolver(&mySys,0,100,0.001,0.2,device,"EulerForward",FDMScheme);
    Evolver evolver(&mySys,0,100,0.01,0.2,device,"RK4",FDMScheme);

    // Running simulations
    evolver.run(); 

    
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};

