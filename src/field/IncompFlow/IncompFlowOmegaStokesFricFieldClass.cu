#ifndef INCOMPFLOWOMEGASTOKESFRICFIELDCLASS_CU
#define INCOMPFLOWOMEGASTOKESFRICFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "IncompFlowOmegaStokesFricFieldClass.h"
#include "IncompFlowOmegaFieldClassGPU.cu"


using namespace std;

// =============================================================
// Constructors
IncompFlowOmegaStokesFricField::IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=0;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.location = "both";
    traits_host.expo_data = "on";    
    initFieldAddi();
};


// -------------------------------------------------------------
IncompFlowOmegaStokesFricField::IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initFieldAddi();
};


// -------------------------------------------------------------
IncompFlowOmegaStokesFricField::IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initFieldAddi ();
};

// -------------------------------------------------------------
IncompFlowOmegaStokesFricField::IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initFieldAddi ();
};

// -------------------------------------------------------------
void IncompFlowOmegaStokesFricField::getVelocity(int i_field) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;

    getIncompFlowStreamGPU<<<Ny,Nx>>>(phi_complex, (*ptr_phi).f[i_field], f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 0);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_FORWARD);
    getIncompFlowStreamGPU<<<Ny,Nx>>> (phi_complex, (*ptr_phi).f[i_field], f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 1);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_INVERSE);
    getIncompFlowStreamGPU<<<Ny,Nx>>> (phi_complex, (*ptr_phi).f[i_field], f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 2);
    (*ptr_phi).applyBounCondPeriGPU((*ptr_phi).f[i_field]);

    FFuncType d1x=f_func_map_all_dev[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all_dev[{"d1y","CentralDifferenceO4Iso2D"}];
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    getIncompFlowVCoreGPU<<<Ny,Nx>>>((*ptr_phi).f[i_field], (*ptr_vx).f[i_field], (*ptr_vy).f[i_field], d1x, d1y, f_func_args);
    (*ptr_vx).applyBounCondPeriGPU((*ptr_vx).f[i_field]);
    (*ptr_vy).applyBounCondPeriGPU((*ptr_vy).f[i_field]);    
};


//===============================================================
void IncompFlowOmegaStokesFricField::getIncompFlowVCoreCPU(int i_field, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args) {
    
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


//==============================================================
void IncompFlowOmegaStokesFricField::initPoissonSolver() {
  // Creating wavenumber array
    double kx,ky;
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    poisson_k2_host=new real[Nx*Ny];
    cudaMalloc((void **)&poisson_k2_dev, (Nx*Ny)*sizeof(double));
    cudaMalloc((void **)&phi_complex, sizeof(cufftDoubleComplex)*Nx*Ny);
    
    for (int i=0; i<Ny; i++){
        ky = 2*Pi*i/(Ny*dy+0.0);
        if (i>=Ny/2) {
            ky=2*Pi*(i-Ny)/(Ny*dy+0.0);
        }
        for (int j=0; j<Nx; j++){
            kx=2*Pi*j/(Nx*dx+0.0);
            if (j>=Nx/2) {
                kx=2*Pi*(j-Nx)/(Nx*dx+0.0);
            }
            int idx=i*Nx+j;
            poisson_k2_host[idx]=kx*kx+ky*ky;
        }
    }
    poisson_k2_host[0]=1;
    cudaMemcpy(poisson_k2_dev,poisson_k2_host,sizeof(double)*Nx*Ny,cudaMemcpyHostToDevice);  
    cufftPlan2d(&cufftPlan, Ny, Nx, CUFFT_Z2Z);
}


// =============================================================

#endif