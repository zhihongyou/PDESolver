#ifndef INCOMPFLOWCLASS_CU
#define INCOMPFLOWCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "IncompFlowClass.h"

using namespace std;

// =============================================================
// Constructors
IncompFlow::IncompFlow (Mesh* mesh_ptr_t, string name_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=0;
    boun_cond = "periodic";
    init_cond = "sin";
    location = "both";
    expo_data = "on";
    initFields ();
};


// -------------------------------------------------------------
IncompFlow::IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = "sin";
    expo_data = "on";
    initFields ();    
};


// -------------------------------------------------------------
IncompFlow::IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = init_cond_t;
    expo_data = "on";
    initFields ();
};

// -------------------------------------------------------------
IncompFlow::IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = boun_cond_t;
    init_cond = init_cond_t;
    expo_data = expo_data_t;
    initFields ();
};


// -------------------------------------------------------------
void IncompFlow::initFields () {
    setFieldProperties(&omega, name+".omega",priority, init_cond);
    if (priority == 0) {
        setFieldProperties(&vx, name+".vx",-1, init_cond);
        setFieldProperties(&vy, name+".vy",-1, init_cond);    
        setFieldProperties(&phi, name+".phi",-1, "sin");
    } else {
        setFieldProperties(&vx, name+".vx",priority, init_cond);
        setFieldProperties(&vy, name+".vy",priority, init_cond);
        setFieldProperties(&phi, name+".phi",priority, "sin");
    };
    // This can trigger calculating stream function and velocity after
    //   a new omega has been obtained
    omega.specialty="IncompFlowOmega";
    omega.ptr_vx=&vx;
    omega.ptr_vy=&vy;
    omega.ptr_phi=&phi;
    omega.setOmegaField(mesh_ptr, 1);
    // No RHS for vx, vy, phi, as they need special functions to get values
    vx.setRhsTerms({});
    vy.setRhsTerms({});
    phi.setRhsTerms({});
    
    vx.allocField<double>(vx.f_host[0], "cpu");
    vy.allocField<double>(vy.f_host[0], "cpu");
    omega.allocField<double>(omega.f_host[0], "cpu");
    phi.allocField<double>(phi.f_host[0], "cpu");
    setFieldValues();
};


// -------------------------------------------------------------
void IncompFlow::setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t) {
    (*field_ptr).traits_host.mesh_ptr=mesh_ptr;
    (*field_ptr).traits_host.name=field_name;
    (*field_ptr).traits_host.priority=priority_t;
    (*field_ptr).traits_host.boun_cond = boun_cond;
    (*field_ptr).traits_host.init_cond = init_cond_t;
    (*field_ptr).traits_host.expo_data = expo_data;
    (*field_ptr).num_f_funcs=0;
    for (int i=0; i<200; i++) {
        (*field_ptr).f_funcs_host[i]=NULL;
    };
};


// -------------------------------------------------------------
void IncompFlow::setFieldValues () {
    int Nx=vx.gridNumber().x;
    int Ny=vx.gridNumber().y;
    int Nbx=vx.gridNumberBoun().x;
    int Nby=vx.gridNumberBoun().y;
    double dx=vx.gridSize().x;
    double dy=vx.gridSize().y;   
    
    if (vx.traits_host.init_cond=="Gaussian") {
        vx.initFieldGaus(0,10,1);
        vy.initFieldGaus(0,10,1);
    } else if (vx.traits_host.init_cond=="ones") {
        vx.initFieldConst(1);
        vy.initFieldConst(1);
    } else if (vx.traits_host.init_cond=="sin") {
        vx.initFieldSin(0.01,1,0);
        vy.initFieldSin(0.01,1,0);
    };

    FFuncType d1x=f_func_map_all[{"d1x","CentralDifferenceO4Iso2D"}];
    FFuncType d1y=f_func_map_all[{"d1y","CentralDifferenceO4Iso2D"}];
    FFuncArgs f_func_args(Nx, Ny, Nbx, Nby, dx, dy);
    // omega.getVelocity(phi.f_host[0], vx.f_host[0], vy.f_host[0], d1x, d1y, f_func_args);
    vx.applyBounCondPeriCPU(vx.f_host[0]);
    vy.applyBounCondPeriCPU(vy.f_host[0]);
    setOmegaValues();
};


// -------------------------------------------------------------
void IncompFlow::setOmegaValues () {
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


// =============================================================

#endif
