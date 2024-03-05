#ifndef LIVINGLCPOLARFIELDCLASS_CU
#define LIVINGLCPOLARFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "LivingLCPolarFieldClass.h"
#include "LivingLCPolarFieldClassGPU.cu"

using namespace std;

// =============================================================
// Constructors
// -------------------------------------------------------------
LivingLCPolarField::LivingLCPolarField (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=0;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.location = "both";
    traits_host.expo_data = "on";    
    initPolarField ();
};


// -------------------------------------------------------------
LivingLCPolarField::LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = "sin";
    traits_host.expo_data = "on";
    initPolarField ();
};


// -------------------------------------------------------------
LivingLCPolarField::LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = "periodic";
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = "on";
    initPolarField ();
};


// -------------------------------------------------------------
LivingLCPolarField::LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond = boun_cond_t;
    traits_host.init_cond = init_cond_t;
    traits_host.expo_data = expo_data_t;
    initPolarField ();
};


// -------------------------------------------------------------
void LivingLCPolarField::initPolarField() {
    
};


// -------------------------------------------------------------
void LivingLCPolarField::postProcessing(int i_field) {
    // Get velocity from vorticity field
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx=gridSize().x;
    double dy=gridSize().y;

    getLivingLCPxPyThetaGPU<<<Ny,Nx>>>((*ptr_Pxx).f[i_field], f[i_field], (*ptr_px).f[i_field], (*ptr_py).f[i_field], (*ptr_theta).f[i_field], (*ptr_theta_old).f[i_field], Nx, Ny, Nbx, Nby);
    getLivingLCFlipGPU<<<Ny,Nx>>>((*ptr_theta_old).f[i_field], (*ptr_theta).f[i_field], (*ptr_flip).f[i_field], Nx, Ny, Nbx, Nby);
    (*ptr_px).applyBounCondPeriGPU((*ptr_px).f[i_field]);
    (*ptr_py).applyBounCondPeriGPU((*ptr_py).f[i_field]);
    (*ptr_theta).applyBounCondPeriGPU((*ptr_theta).f[i_field]);
    (*ptr_flip).applyBounCondPeriGPU((*ptr_flip).f[i_field]);
};


// =============================================================

#endif
