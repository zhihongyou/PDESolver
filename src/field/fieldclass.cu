#ifndef FIELDCLASS_CU
#define FIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "fieldclass.h"
#include "fieldclassGPU.cu"
#include "fieldclass_fieldfunctionsGPU.cu"
#include "fieldclass_boundaryconditionGPU.cu"
#include "fieldclass_initialconditionGPU.cu"
#include "fieldclass_fieldfunctions.cu"
#include "fieldclass_boundarycondition.cu"
#include "fieldclass_initialcondition.cu"
#include "fieldclass_export.cu"


using namespace std;


// ==============================================================
// Constructors
Field::Field (Mesh* mesh_ptr_t, string name_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.rank=1;
    traits_host.priority=0;
    traits_host.boun_cond="periodic";
    traits_host.init_cond="sin";
    traits_host.expo_data="on";
    initFieldAddi();
};


// ==============================================================
// Constructors
Field::Field (Mesh* mesh_ptr_t, string name_t, int priority_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.rank=1;
    traits_host.priority=priority_t;
    traits_host.boun_cond="periodic";
    traits_host.init_cond="sin";
    traits_host.expo_data="on";
    initFieldAddi();
};


// -------------------------------------------------------------
Field::Field (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.rank=0;
    traits_host.priority=priority_t;
    traits_host.boun_cond="periodic";
    traits_host.init_cond=init_cond_t;
    traits_host.expo_data="on";
    initFieldAddi();
};


// -------------------------------------------------------------
Field::Field (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.rank=0;
    traits_host.priority=priority_t;
    traits_host.boun_cond=boun_cond_t;
    traits_host.init_cond=init_cond_t;
    traits_host.expo_data=expo_data_t;
    initFieldAddi();
};


// =============================================================
// Other field class methods
// -------------------------------------------------------------
void Field::initFieldAddi() {
    // Additional setups when initializing fields
    allocField<double>(f_host[0], "cpu");
    if (traits_host.init_cond=="Gaussian") {
        initFieldGaus(0,10,1);
    } else if (traits_host.init_cond=="zeros") {
        initFieldConst(0);
    } else if (traits_host.init_cond=="ones") {
        initFieldConst(1);
    } else if (traits_host.init_cond=="sin") {
        initFieldSin(0.01,2,0);
    };
    num_f_funcs=0;
    for (int i=0; i<200; i++) {
        f_funcs_host[i]=NULL;
    };
};


// -------------------------------------------------------------
void Field::setFieldConstCPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    setFieldConstCPUCore(f_t, f_val, Nx, Ny, Nbx, Nby);
};


// ------------------------------------------------------------
void Field::setFieldConstGPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    setFieldConstGPUCore<<<Ny,Nx>>>(f_t, f_val, Nx, Ny, Nbx, Nby);
};


// -----------------------------------------------------------
void Field::setRhsTerms(vector<rhsTerm> rhs_terms_t) {
    rhs_terms=rhs_terms_t;
    for (int i=0; i<rhs_terms.size(); i++) {        
        for (int j=0; j<rhs_terms[i].f_funcs.size(); j++) {
            rhs_terms[i].f_funcs[j].f_func_args.Nx=gridNumber().x;
            rhs_terms[i].f_funcs[j].f_func_args.Ny=gridNumber().y;
            rhs_terms[i].f_funcs[j].f_func_args.Nbx=gridNumberBoun().x;
            rhs_terms[i].f_funcs[j].f_func_args.Nby=gridNumberBoun().y;
            rhs_terms[i].f_funcs[j].f_func_args.di=gridNumber().x+2*gridNumberBoun().x;
            rhs_terms[i].f_funcs[j].f_func_args.dj=1;
            rhs_terms[i].f_funcs[j].f_func_args.dx=gridSize().x;
            rhs_terms[i].f_funcs[j].f_func_args.dy=gridSize().y;
        };
    };
};


// -------------------------------------------------------------
// string Field::equation() {
//     string eqn;
//     if (priority==0) {
//         eqn="p_t "+name()+"=";
//     } else {
//         eqn=name()+"=";
//     };
//     for (auto rhs_term_i : rhsTerms()) {
//         if (rhs_term_i != rhsTerms().begin()) {
            
//         };
//         for (auto f_func_i : rhs_term_i.f_function) {
//             eqn=eqn+
//                 cout<<"*"<<f_func_i.f_operator<<"("<<(*f_func_i.field_ptr).name() <<")";
//                 evalOperator(f_ptr_i,f_func_i,f_func_ptrs[N_funcs],i_field);
//                 N_funcs+=1;
//             };
//             addRHSTerm(f_ptr_i,i_field,rhs_term_i,f_func_ptrs,N_funcs);
//         };

//         if ((*f_ptr_i).priority()>0 && (*f_ptr_i).bounCond()=="periodic") {
//             if (device=="cpu") {
//                 (*f_ptr_i).applyBounCondPeriCPU(i_field);
//             } else if (device=="gpu"){
//                 (*f_ptr_i).applyBounCondPeriGPU(i_field);
//             };
//         };
// };

// -------------------------------------------------------------
void Field::updateAnyFieldDev (double* f_dev_ptr, double * f_host_ptr) {
    // Copy any field data from CPU to GPU

    allocField<double>(f_dev_ptr, "gpu");
    cudaMemcpy(f_dev_ptr, f_host_ptr, gridNumberAll()*sizeof(double),cudaMemcpyHostToDevice);
};

// -------------------------------------------------------------
void Field::updateAnyFieldHost (double* f_host_ptr, double * f_dev_ptr) {
    // Copy any field data from GPU to CPU
    cudaMemcpy(f_host_ptr, f_dev_ptr, gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
};

// -------------------------------------------------------------
void Field::updateMainFieldDev () {
    // Copy main field data from CPU to GPU
    allocField<double>(f[0], "gpu");
    cudaMemcpy(f[0], f_host[0], gridNumberAll()*sizeof(double),cudaMemcpyHostToDevice);
};

// -------------------------------------------------------------
void Field::updateMainFieldHost () {
    // Copy main field data from GPU to CPU
    allocField<double>(f[0], "cpu");
    cudaMemcpy(f_host[0], f[0], gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
};

// -------------------------------------------------------------
template <typename T>
void Field::allocField (T* &f_t, string location) {
    // Allocate memory for fields
    if (f_t==NULL) {
        if (location=="cpu") {
            f_t=new T[gridNumberAll()];
        } else if (location=="gpu") {
            cudaMalloc(&f_t, gridNumberAll()*sizeof(T));
        };
    };    
};

// -------------------------------------------------------------
double* Field::getFFuncPtr(string f_operator) {
    // Get pointer to a field function
    double* f_func_ptr;
    for (int i=0; i<num_f_funcs; i++) {
        if (f_funcs_rhs[i].f_operator==f_operator) {
            f_func_ptr=f_funcs_host[f_funcs_rhs[i].f_func_idx];
        };
    };
    
    return f_func_ptr;
        
};


// -------------------------------------------------------------
void Field::addFunctoRHS(FFuncDef f_func_i, string device, string func_scheme) {
    // Add each function on the RHS to a list of function to be evaluated
    FFuncItem f_func_item;
    f_func_item.f_operator=f_func_i.f_operator;
    f_func_item.f_func_idx=num_f_funcs;
    f_func_item.f_func_args=f_func_i.f_func_args;    
    
    if (device=="cpu") {
        // Allocate a memory to store the function field
        allocField<double>(f_funcs_host[num_f_funcs], "cpu");        
        if (f_func_map_all[{f_func_i.f_operator,func_scheme}]==0) {
            // Onsite function
            f_func_item.f_func=f_func_map_all[{f_func_i.f_operator,""}];
        } else {
            // Finite difference operator
            f_func_item.f_func=f_func_map_all[{f_func_i.f_operator,func_scheme}];
        };
    } else if (device=="gpu") {
        allocField<double>(f_funcs_host[num_f_funcs], "gpu");
        if (f_func_map_all[{f_func_i.f_operator,func_scheme}]==0) {
            f_func_item.f_func=f_func_map_all_dev[{f_func_i.f_operator,""}];
        } else {
            f_func_item.f_func=f_func_map_all_dev[{f_func_i.f_operator,func_scheme}];
        };
    };
    f_funcs_rhs[num_f_funcs]=f_func_item;
    num_f_funcs+=1;
};


// -------------------------------------------------------------
void Field::getRHSAddi(int i_field) {
    // cout << "Getting additional rhs." <<endl;
};


// =============================================================

#endif
