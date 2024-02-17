#ifndef LAPLACENFEQFIELDCLASS_CU
#define LAPLACENFEQFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "LaplaceNFEqFieldClass.h"
#include "LaplaceNFEqFieldClassGPU.cu"


using namespace std;

// ==============================================================
// Constructors
// =============================================================
// Constructors
LaplaceNFEqField::LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t) {
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
LaplaceNFEqField::LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initFieldAddi();
};


// -------------------------------------------------------------
LaplaceNFEqField::LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initFieldAddi ();
};

// -------------------------------------------------------------
LaplaceNFEqField::LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initFieldAddi ();
};


// -------------------------------------------------------------
void LaplaceNFEqField::initFieldAddi () {    
    allocField<double>(f_host[0], "cpu");
    num_f_funcs=0;
    for (int i=0; i<200; i++) {
        f_funcs_host[i]=NULL;
    };
    specialty="LaplaceNFEqField";
    initLaplaceNFSolver();
    cout << "Initiating a LaplaceNFEqField: " << name()<<"."<<endl;
};


// --------------------------------------------------------------
void LaplaceNFEqField::getRHSAddi (int i_field) {
    solveLaplaceNFEq(i_field);
    // cout << "Solving the poisson equation." <<endl;
};


// --------------------------------------------------------------
void LaplaceNFEqField::solveLaplaceNFEq (int i_field) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, f[i_field], f[i_field], k2s_dev, Nx, Ny, Nbx, Nby, 0);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_FORWARD);
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, f[i_field], rhs[i_field], k2s_dev, Nx, Ny, Nbx, Nby, 1);
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_INVERSE);
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, f[i_field], f[i_field], k2s_dev, Nx, Ny, Nbx, Nby, 2);
    applyBounCondPeriGPU(f[i_field]);    
};


//===============================================================
void LaplaceNFEqField::initLaplaceNFSolver() {
  // Creating wavenumber array
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    k2s_host=new real[Nx*Ny];
    cudaMalloc((void **)&k2s_dev, (Nx*Ny)*sizeof(double));
    cudaMalloc((void **)&phi_complex, sizeof(cufftDoubleComplex)*Nx*Ny);
    setk2s();    
    cufftPlan2d(&cufftPlan, Ny, Nx, CUFFT_Z2Z);
}

//==============================================================
void LaplaceNFEqField::setk2s() {
    double kx,ky;
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    double dx=gridSize().x;
    double dy=gridSize().y;
    
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
	    k2s_host[idx]=prefactors[0];
	    for (int k=1; k<=max_power; k++) {
	      k2s_host[idx]=k2s_host[idx]+pow(-kx*kx-ky*ky, k);
	    }
        }
    }
    if (abs(prefactors[0])<0.0000000000000001) {
      k2s_host[0]=1;
    } else {
      k2s_host[0]=prefactors[0];
    };
    cudaMemcpy(k2s_dev,k2s_host,sizeof(double)*Nx*Ny,cudaMemcpyHostToDevice);
};

//===============================================================
void LaplaceNFEqField::setSolver (int max_power_t, double* prefactors_t) {
    max_power=max_power_t;
    for (i=0; i<=max_power; i++) {
      prefactors[i]=prefactors_t[i];
    }
    setk2s();
    cout << "LaplaceNFEqField Solver has been reset!" <<endl;
};

// ==============================================================

#endif
