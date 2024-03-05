#ifndef LIVINGLCCLASS_CU
#define LIVINGLCCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "LivingLCClass.h"

using namespace std;

// =============================================================
// Constructors
LivingLC::LivingLC (Mesh* mesh_ptr_t, string name_t) {
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
LivingLC::LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = "sin";
    expo_data = "on";
    initFields ();    
};


// -------------------------------------------------------------
LivingLC::LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = "periodic";
    init_cond = init_cond_t;
    expo_data = "on";
    initFields ();
};

// -------------------------------------------------------------
LivingLC::LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=priority_t;
    boun_cond = boun_cond_t;
    init_cond = init_cond_t;
    expo_data = expo_data_t;
    initFields ();
};


// -------------------------------------------------------------
void LivingLC::initFields () {
    setFieldProperties(&Pxx, name+".Pxx", priority, init_cond);
    setFieldProperties(&Pxy, name+".Pxy", priority, init_cond);
    setFieldProperties(&theta_old, name+".theta_old",-1, init_cond);
    setFieldProperties(&theta, name+".theta",-1, init_cond);
    setFieldProperties(&px, name+".px",-1, "sin");
    setFieldProperties(&py, name+".py",-1, "sin");
    setFieldProperties(&flip, name+".flip",-1, "sin");
    flip.initFieldConst(1, 0);
    
    Pxy.specialty="LivingLCPolar";
    Pxy.ptr_Pxx=&Pxx;
    Pxy.ptr_theta_old=&theta_old;
    Pxy.ptr_theta=&theta;
    Pxy.ptr_px=&px;
    Pxy.ptr_py=&py;
    Pxy.ptr_flip=&flip;
    
    // No RHS for vx, vy, phi, as they need special functions to get values
    theta_old.setRhsTerms({});
    theta.setRhsTerms({});
    px.setRhsTerms({});
    py.setRhsTerms({});
    flip.setRhsTerms({});   
};


// -------------------------------------------------------------
void LivingLC::setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t) {
    (*field_ptr).traits_host.mesh_ptr=mesh_ptr;
    (*field_ptr).traits_host.name=field_name;
    (*field_ptr).traits_host.priority=priority_t;
    (*field_ptr).traits_host.boun_cond = boun_cond;
    (*field_ptr).traits_host.init_cond = init_cond_t;
    (*field_ptr).traits_host.expo_data = expo_data;
    (*field_ptr).num_f_funcs=0;
    
    (*field_ptr).initFieldAddi();
    for (int i=0; i<200; i++) {
        (*field_ptr).f_funcs_host[i]=NULL;
    };
};


// =============================================================

#endif
