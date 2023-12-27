#ifndef POISSONEQFIELDCLASS_CU
#define POISSONEQFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "poissonEqFieldclass.h"
#include "poissonEqFieldclassGPU.cu"


using namespace std;

// ==============================================================
// Constructors
// =============================================================
// Constructors
poissonEqField::poissonEqField (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=1;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.location = "both";
    traits_host.expo_data = "on";
    initFieldAddi ();
};


// -------------------------------------------------------------
poissonEqField::poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initFieldAddi ();    
};


// -------------------------------------------------------------
poissonEqField::poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initFieldAddi ();
};

// -------------------------------------------------------------
poissonEqField::poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initFieldAddi ();
};


// -------------------------------------------------------------
void poissonEqField::initFieldAddi () {    
    allocField<double>(f_host[0], "cpu");
    num_f_funcs=0;
    for (int i=0; i<200; i++) {
        f_funcs_host[i]=NULL;
    };
    specialty="poissonEqField";
    initPoissonSolver();
    cout << "Initiating a poissonEqField: " << name()<<"."<<endl;
};


// --------------------------------------------------------------
// void poissonEqField::getRHSAddi (int i_field) {
    // solvePoissonEq(i_field);
    // cout << "Solving the poisson equation." <<endl;
// };


// --------------------------------------------------------------
void poissonEqField::solvePoissonEq (int i_field) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    
    getPhiGPU<<<Ny,Nx>>>(phi_complex, f[i_field], f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 0);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_FORWARD);
    getPhiGPU<<<Ny,Nx>>>(phi_complex, f[i_field], rhs[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 1);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_INVERSE);
    getPhiGPU<<<Ny,Nx>>>(phi_complex, f[i_field], f[i_field], poisson_k2_dev, Nx, Ny, Nbx, Nby, 2);
    applyBounCondPeriGPU(f[i_field]);    
};


//===============================================================
void poissonEqField::initPoissonSolver() {
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


// ==============================================================

#endif
