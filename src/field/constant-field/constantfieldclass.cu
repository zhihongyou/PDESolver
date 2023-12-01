#ifndef CONSTANTFIELDCLASS_CU
#define CONSTANTFIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include <cufft.h>
#include <cufftXt.h>
#include "constantfieldclass.h"


using namespace std; 

// ===================================================================
// Constructors
ConstantField::ConstantField (Mesh* mesh_ptr_t, string name_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    priority=-1;
    boun_cond = "periodic";
    init_cond = "";
    location = "both";
    expo_data = "off";
    initFields ();
};


// -------------------------------------------------------------------
void ConstantField::initFields () {
    setFieldProperties(&one, name+".one",-1);
    setFieldProperties(&zero, name+".zero",-1);
    setFieldProperties(&pi, name+".pi",-1);
    
    one.setRhsTerms({});
    zero.setRhsTerms({});
    pi.setRhsTerms({});
    
    one.allocField<double>(one.f_host[0], "cpu");
    zero.allocField<double>(zero.f_host[0], "cpu");
    pi.allocField<double>(pi.f_host[0], "cpu");

    setFieldValues();
};


// -------------------------------------------------------------------
void ConstantField::setFieldProperties (Field* field_ptr, string field_name, int priority_t) {
    (*field_ptr).traits_host.mesh_ptr=mesh_ptr;
    (*field_ptr).traits_host.name=field_name;
    (*field_ptr).traits_host.priority=priority_t;
    (*field_ptr).traits_host.boun_cond = boun_cond;
    (*field_ptr).traits_host.expo_data = "off";
    (*field_ptr).num_f_funcs=0;
    for (int i=0; i<200; i++) {
        (*field_ptr).f_funcs_host[i]=NULL;
    };
};


// -------------------------------------------------------------------
void ConstantField::setFieldValues () {
    int Nx=one.gridNumber().x;
    int Ny=one.gridNumber().y;
    int Nbx=one.gridNumberBoun().x;
    int Nby=one.gridNumberBoun().y;
    
    for (int i=0; i<Ny; i++) {
        for (int j=0; j<Nx; j++) {
            int idx=(i+Nby)*(Nx+2*Nbx)+j+Nbx;
            one.f_host[0][idx]=1;
            zero.f_host[0][idx]=0;
            pi.f_host[0][idx]=3.14159265358979323846;
        };
    };
};


// =====================================================================

#endif
