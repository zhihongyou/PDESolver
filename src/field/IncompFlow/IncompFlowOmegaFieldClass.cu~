#ifndef INCOMPFLOWOMEGAFIELDCLASS_CU
#define INCOMPFLOWOMEGAFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "IncompFlowOmegaFieldClass.h"
#include "IncompFlowOmegaFieldClassGPU.cu"

using namespace std;

// =============================================================
// Constructors
// -------------------------------------------------------------
IncompFlowOmegaField::IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=0;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.location = "both";
    traits_host.expo_data = "on";    
    initOmegaField ();
};


// -------------------------------------------------------------
IncompFlowOmegaField::IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initOmegaField ();
};


// -------------------------------------------------------------
IncompFlowOmegaField::IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initOmegaField ();
};


// -------------------------------------------------------------
IncompFlowOmegaField::IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initOmegaField ();
};


// -------------------------------------------------------------
void IncompFlowOmegaField::initOmegaField() {
    // Empty constructor requires specifying mesh before initialization.
    LaplSolverW2Phi.mesh_ptr=traits_host.mesh_ptr;
    // Initiate a Poisson solver
    LaplSolverW2Phi.initLaplaceNFEqSolver(traits_host.mesh_ptr);
    double prefactors_t[2]={0, -1};
    LaplSolverW2Phi.setLaplaceNFEqSolver(1, prefactors_t);
    
    // Customize solver for Frictional-Stokes equation
    // The solver solves the Stokes equation by default
    LaplSolverW.mesh_ptr=traits_host.mesh_ptr;
    // Initiate a Poisson solver
    LaplSolverW.initLaplaceNFEqSolver(traits_host.mesh_ptr);
    setOmegaEq(1, 0, 1);
};


// -------------------------------------------------------------
void IncompFlowOmegaField::setOmegaEq (int omegaEqType_t, double gamma_t, double eta_t) {
    // Initiate a Poisson solver
    omegaEqType=omegaEqType_t;
    gamma=gamma_t;
    eta=eta_t;
    double prefactors_t[2]={gamma, -1*eta};
    LaplSolverW.setLaplaceNFEqSolver(1, prefactors_t);
};


// -------------------------------------------------------------
void IncompFlowOmegaField::getOmegaAddi(int i_field) {
    if (omegaEqType>1) {
        LaplSolverW.solveLaplaceNFEq(f[i_field],f[i_field]);
    };
    applyBounCondPeriGPU(f[i_field]);
}


// -------------------------------------------------------------
void IncompFlowOmegaField::getVelocity(int i_field) {
    // Get velocity from vorticity field
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;

    // Solve the Poission Eq to get the stream function
    LaplSolverW2Phi.solveLaplaceNFEq((*ptr_phi).f[i_field], f[i_field]);
    (*ptr_phi).applyBounCondPeriGPU((*ptr_phi).f[i_field]);

    FFuncType d1x=f_func_map_all_dev[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all_dev[{"d1y","CentralDifferenceO4Iso2D"}];
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    // Get velocity from stream function
    getIncompFlowPhi2VGPU<<<Ny,Nx>>>((*ptr_phi).f[i_field], (*ptr_vx).f[i_field], (*ptr_vy).f[i_field], d1x, d1y, f_func_args);
    (*ptr_vx).applyBounCondPeriGPU((*ptr_vx).f[i_field]);
    (*ptr_vy).applyBounCondPeriGPU((*ptr_vy).f[i_field]);    
};


//===============================================================
void IncompFlowOmegaField::getIncompFlowVCoreCPU(int i_field, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args) {
    // Get velocity from stream function, cpu version
    int Nx=f_func_args.Nx;
    int Ny=f_func_args.Ny;
    int Nbx=f_func_args.Nbx;
    int Nby=f_func_args.Nby;
    double* vx=(*ptr_vx).f[i_field];
    double* vy=(*ptr_vy).f[i_field];
    double* phi=(*ptr_phi).f[i_field];
    
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            int idx=(i+Nby)*(Nx+2*Nbx)+j+Nbx;
            vx[idx]=d1y(phi,idx,f_func_args);
            vy[idx]=-d1x(phi,idx,f_func_args);
        };
    };
}


// =============================================================

#endif
