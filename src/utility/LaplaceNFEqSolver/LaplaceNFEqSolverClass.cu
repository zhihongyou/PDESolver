#ifndef LAPLACENFEQSOLVERCLASS_CU
#define LAPLACENFEQSOLVERCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "LaplaceNFEqSolverClass.h"


using namespace std;


// =============================================================
// Constructors
LaplaceNFEqSolver::LaplaceNFEqSolver () {};


// ------------------------------------------------------------
LaplaceNFEqSolver::LaplaceNFEqSolver (Mesh* mesh_ptr_t) {
    // cout <<"haha"<<endl;
    initLaplaceNFEqSolver(mesh_ptr_t);
    // cout <<"haha"<<endl;
};


// ------------------------------------------------------------
LaplaceNFEqSolver::LaplaceNFEqSolver (Mesh* mesh_ptr_t, string name_t) {
    name=name_t;
    initLaplaceNFEqSolver(mesh_ptr_t);
};


//===============================================================
void LaplaceNFEqSolver::initLaplaceNFEqSolver(Mesh* mesh_ptr_t) {
    int Nx=(*mesh_ptr_t).host.grid_number.x;    
    int Ny=(*mesh_ptr_t).host.grid_number.y;
    if (k2s_host==NULL) {
        k2s_host=new double[Nx*Ny];
        cudaMalloc((void **)&k2s_dev, (Nx*Ny)*sizeof(double));
        cudaMalloc((void **)&phi_complex, sizeof(cufftDoubleComplex)*Nx*Ny);
        cufftPlan2d(&cufftPlan, Ny, Nx, CUFFT_Z2Z);
    }
    mesh_ptr=mesh_ptr_t;
    double prefactors_t[2]={0, 1};
    setLaplaceNFEqSolver(1, prefactors_t);
}


// ============================================================
void LaplaceNFEqSolver::setLaplaceNFEqSolver(int max_power_t, double* prefactors_t) {
    max_power=max_power_t;
    for (int i=0; i<=max_power; i++) {
        prefactors[i]=prefactors_t[i];
    }
    setk2s();
};


//==============================================================
void LaplaceNFEqSolver::setk2s() {
    double kx,ky;
    int Nx=(*mesh_ptr).host.grid_number.x;
    int Ny=(*mesh_ptr).host.grid_number.y;
    double dx=(*mesh_ptr).host.grid_size.x;
    double dy=(*mesh_ptr).host.grid_size.y;
    
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
	      k2s_host[idx]=k2s_host[idx]+prefactors[k]*pow(-kx*kx-ky*ky, k);
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


// --------------------------------------------------------------
void LaplaceNFEqSolver::solveLaplaceNFEq (double* phi, double* f) {
    int Nx=(*mesh_ptr).host.grid_number.x;
    int Ny=(*mesh_ptr).host.grid_number.y;
    int Nbx=(*mesh_ptr).host.grid_number_boun.x;
    int Nby=(*mesh_ptr).host.grid_number_boun.y;

    // Assign f to the real part of phi_complex
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, phi, f, k2s_dev, Nx, Ny, Nbx, Nby, 0);
    // Fourier transform phi_complex
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_FORWARD);
    // Set phi_complex=phi_complex/k2s
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, phi, f, k2s_dev, Nx, Ny, Nbx, Nby, 1);
    // Inverse Fourier transform phi_complex
    cufftExecZ2Z(cufftPlan,phi_complex,phi_complex,CUFFT_INVERSE);
    // Assign the real part of phi_complex to phi
    solveLaplaceNFEqCoreGPU<<<Ny,Nx>>>(phi_complex, phi, f, k2s_dev, Nx, Ny, Nbx, Nby, 2);
    // applyBounCondPeriGPU(phi);    
};


// ==============================================================

#endif
