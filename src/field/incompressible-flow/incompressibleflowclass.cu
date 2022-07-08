#ifndef INCOMPRESSIBLEFLOWCLASS_CU
#define INCOMPRESSIBLEFLOWCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "incompressibleflowclass.h"
#include "incompressibleflowclassGPU.cu"


using namespace std; 

// ===================================================================
// Constructors
IncompressibleFlow::IncompressibleFlow (Mesh* mesh_ptr_t, string name_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=0;
    boun_cond = "periodic";
    init_cond = "sin";
    location = "both";
    expo_data = "on";
    initFields ();
};


// -------------------------------------------------------------------
IncompressibleFlow::IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = "sin";
    expo_data = "on";
    initFields ();    
};


// -------------------------------------------------------------------
IncompressibleFlow::IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = init_cond_t;
    expo_data = "on";
    initFields ();
};

// -------------------------------------------------------------------
IncompressibleFlow::IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = boun_cond_t;
    init_cond = init_cond_t;
    expo_data = expo_data_t;
    initFields ();
};


// -------------------------------------------------------------------
void IncompressibleFlow::initFields () {
    setFieldProperties(&vx, name+".vx",-1);
    setFieldProperties(&vy, name+".vy",-1);
    setFieldProperties(&omega, name+".omega",priority);
    setFieldProperties(&phi, name+".phi",-1);
    // This can trigger calculating stream function and velocity after
    //   a new omega has been obtained
    omega.specialty="IncompressibleFlow.omega";
    // No RHS for vx, vy, phi, as they need special functions to get values
    vx.setRhsTerms({});
    vy.setRhsTerms({});
    phi.setRhsTerms({});
    
    vx.allocField<double>(vx.f_host[0], "cpu");
    vy.allocField<double>(vy.f_host[0], "cpu");
    omega.allocField<double>(omega.f_host[0], "cpu");
    phi.allocField<double>(phi.f_host[0], "cpu");
    setFieldValues(init_cond);
    
    initPoissonSolver();
};


// -------------------------------------------------------------------
void IncompressibleFlow::setFieldProperties (Field* field_ptr, string field_name, int priority_t) {
    (*field_ptr).traits_host.mesh_ptr=mesh_ptr;
    (*field_ptr).traits_host.name=field_name;
    (*field_ptr).traits_host.priority=priority_t;
    (*field_ptr).traits_host.boun_cond = boun_cond;
    (*field_ptr).traits_host.init_cond = init_cond;
    (*field_ptr).traits_host.expo_data = expo_data;
    (*field_ptr).num_f_funcs=0;
    for (int i=0; i<200; i++) {
        (*field_ptr).f_funcs_host[i]=NULL;
    };
};


// -------------------------------------------------------------------
void IncompressibleFlow::setFieldValues (string init_cond_t) {
    int Nx=vx.gridNumber().x;
    int Ny=vx.gridNumber().y;
    int Nbx=vx.gridNumberBoun().x;
    int Nby=vx.gridNumberBoun().y;
    double dx=vx.gridSize().x;
    double dy=vx.gridSize().y;   
    
    if (vx.traits_host.init_cond=="Gaussian") {
        phi.initFieldGaus(0,10,1);
    } else if (vx.traits_host.init_cond=="ones") {
        phi.initFieldConst(1);
    } else if (vx.traits_host.init_cond=="sin") {
        phi.initFieldSin(0.01,1,0);
    };

    FFuncType d1x=f_func_map_all[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all[{"d1y","CentralDifferenceO4Iso2D"}];
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    getVCoreCPU(phi.f_host[0], vx.f_host[0], vy.f_host[0], d1x, d1y, f_func_args);
    vx.applyBounCondPeriCPU(vx.f_host[0]);
    vy.applyBounCondPeriCPU(vy.f_host[0]);
    setOmegaValues();
};


// --------------------------------------------------------------------
void IncompressibleFlow::setOmegaValues () {
    int Nx=vx.gridNumber().x;
    int Ny=vx.gridNumber().y;
    int Nbx=vx.gridNumberBoun().x;
    int Nby=vx.gridNumberBoun().y;
    double dx=vx.gridSize().x;
    double dy=vx.gridSize().y;
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    FFuncType d1x=f_func_map_all[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all[{"d1y","CentralDifferenceO4Iso2D"}];
    
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            int idx=(i+Nby)*(Nx+2*Nbx)+j+Nbx;
            omega.f_host[0][idx]=d1x(vy.f_host[0],idx,f_func_args)-d1y(vx.f_host[0],idx,f_func_args);
        };
    };
    omega.applyBounCondPeriCPU(omega.f_host[0]);
};


// -------------------------------------------------------------------
void IncompressibleFlow::getVelocity (int i_field) {
    int Nx=vx.gridNumber().x;
    int Ny=vx.gridNumber().y;
    int Nbx=vx.gridNumberBoun().x;
    int Nby=vx.gridNumberBoun().y;
    double dx=vx.gridSize().x;
    double dy=vx.gridSize().y;

    getStreamGPU<<<Ny,Nx>>>(phi_complex, phi.f[i_field], omega.f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 0);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_FORWARD);
    getStreamGPU<<<Ny,Nx>>> (phi_complex, phi.f[i_field], omega.f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 1);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_INVERSE);
    getStreamGPU<<<Ny,Nx>>> (phi_complex, phi.f[i_field], omega.f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 2);
    phi.applyBounCondPeriGPU(phi.f[i_field]);

    FFuncType d1x=f_func_map_all_dev[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all_dev[{"d1y","CentralDifferenceO4Iso2D"}];
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    getVCoreGPU<<<Ny,Nx>>>(phi.f[i_field], vx.f[i_field], vy.f[i_field], d1x, d1y, f_func_args);
    vx.applyBounCondPeriGPU(vx.f[i_field]);
    vy.applyBounCondPeriGPU(vy.f[i_field]);    
};


//=======================================================================
void IncompressibleFlow::getVCoreCPU(double* phi, double* vx, double* vy, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args) {

    int Nx=f_func_args.Nx;
    int Ny=f_func_args.Ny;
    int Nbx=f_func_args.Nbx;
    int Nby=f_func_args.Nby;
    
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            int idx=(i+Nby)*(Nx+2*Nbx)+j+Nbx;
            vx[idx]=d1y(phi,idx,f_func_args);
            vy[idx]=-d1x(phi,idx,f_func_args);
        };
    };
}


//=======================================================================
void IncompressibleFlow::initPoissonSolver() {
  // Creating wavenumber array
    double kx,ky;
    int Nx=vx.gridNumber().x;
    int Ny=vx.gridNumber().y;
    double dx=vx.gridSize().x;
    double dy=vx.gridSize().y;
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


// =====================================================================

#endif
