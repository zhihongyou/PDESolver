#ifndef LAPLACENFEQFIELDCLASS_CU
#define LAPLACENFEQFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "LaplaceNFEqFieldClass.h"


using namespace std;

// ==============================================================
// Constructors
// -------------------------------------------------------------
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
    traits_host.priority=1;
    specialty="LaplaceNFEqField";
    // The equation is be default a Poisson equation.
    double prefactors_t[2]={0,1};
    LaplNFEqSolver.initLaplaceNFEqSolver(traits_host.mesh_ptr);
    setLaplaceNFEq(1, prefactors_t);
    cout << "Initiating a LaplaceNFEqField: " << name()<<"."<<endl;
};


// -------------------------------------------------------------
void LaplaceNFEqField::setLaplaceNFEq (int max_power_t, double* prefactors_t) {
    max_power=max_power_t;
    for (int i=0; i<=max_power; i++) {
        prefactors[i]=prefactors_t[i];
    }
    LaplNFEqSolver.setLaplaceNFEqSolver(max_power, prefactors);
};


// --------------------------------------------------------------
void LaplaceNFEqField::solveLaplaceNFEq (int i_field) {
    LaplNFEqSolver.solveLaplaceNFEq(f[i_field],f[i_field]);
    // cout << "Solving the poisson equation." <<endl;
};


// ==============================================================

#endif
