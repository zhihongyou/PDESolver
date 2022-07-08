#ifndef EVOLVERCLASS_CU
#define EVOLVERCLASS_CU

#include <iostream> 
#include <vector>
#include "evolverclass.h"
#include "evolverclassGPU.cu"
#include "evolverclass_integrator.cu"

using namespace std;


// ======================================================================

// ----------------------------------------------------------------------
// Update system over time. This is the global control.
void Evolver::run() {
    cout <<"Start running simulation ..." <<endl;
    time_start_wall = clock();
    
    // Initiate function calls and pointers.
    initEvolver();
    
    // Export initial configuration.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).expoData() =="on") {
            (*f_ptr_i).export_conf("0",device,1);
        };
    };    

    // Running through time
    for (time_now=time_start; time_now<=time_stop+0.1*time_step; time_now+=time_step) {

        // Move a step forward.
        if (scheme=="EulerForward") {
            EulerForward();
        } else if (scheme=="RK4") {
            RK4();
        };

        // Reach export time.
        if ((int) (time_now/time_export) > (int) ((time_now-time_step)/time_export)) {
            // Export field configuration.
            int te=floor(time_now/time_export);    
            string str_t=to_string(te);
            for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
                if ((*f_ptr_i).expoData() =="on") {
                    (*f_ptr_i).export_conf(str_t,device,1);
                    // (*f_ptr_i).export_conf_any((*f_ptr_i).rhs[0],"dt"+(*f_ptr_i).name(),str_t, device, 1);
                };        
            };
            // Show progress of simulation.
            showProgress();
        };
        
    };
    // Simulation finished.
    time_stop_wall = clock();
    std::cout.flush();
    cout <<"Simulation finished. Total spent:";
    cout<<(time_stop_wall-time_start_wall)/CLOCKS_PER_SEC<<" seconds.";
    cout <<"                                                        ";
    cout<<endl;
};


// ----------------------------------------------------------------------
void Evolver::initEvolver() {
    // Initialize field function map used to identify functions to call
    //   based on the operator name.
    setFFuncMap();
    
    // Initialize fields
    initFields();
    
    // Initialize 
    initRHSs();
};


// ----------------------------------------------------------------------
// Allocate memory for fields 
// Location of memory depends on device of the evolver.
void Evolver::initFields () {            

    // Initiate fields and LHS&RHS
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        // Set numerical scheme for each field.
        (*f_ptr_i).FDMScheme=FDMScheme;
        
        int num_grid=(*f_ptr_i).gridNumberAll();
        if (device=="cpu") {            
            for (int i_f_copy=0; i_f_copy<num_field_copy; i_f_copy++) {
                (*f_ptr_i).f[i_f_copy]=new double[num_grid];
                if ((*f_ptr_i).priority()==0) {
                    (*f_ptr_i).rhs[i_f_copy]=new double[num_grid];
                    (*f_ptr_i).lhs[i_f_copy]=new double[num_grid];
                };
            };
            // Copy data from f_host to device
            for (int idx=0; idx<num_grid; idx++) {
                (*f_ptr_i).f[0][idx]=(*f_ptr_i).f_host[0][idx];
            };
            
        } else if (device=="gpu") {
            (*f_ptr_i).f_temp_host=new double[(*f_ptr_i).gridNumberAll()];
            for (int i_f_copy=0; i_f_copy<num_field_copy; i_f_copy++) {
                cudaMalloc(&(*f_ptr_i).f[i_f_copy], num_grid*sizeof(double));
                if ((*f_ptr_i).priority()==0) {
                    cudaMalloc(&(*f_ptr_i).rhs[i_f_copy], num_grid*sizeof(double));
                    cudaMalloc(&(*f_ptr_i).lhs[i_f_copy], num_grid*sizeof(double));
                };
            };
            cudaMemcpy((*f_ptr_i).f[0], (*f_ptr_i).f_host[0], num_grid*sizeof(double),cudaMemcpyHostToDevice);            
        };
    };
    
};


// ----------------------------------------------------------------------
void Evolver::initRHSs() {
    
    // Loop over fields to get functions appearing in their RHS.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        int num_grid=(*f_ptr_i).gridNumberAll();
        int num_terms=(*f_ptr_i).rhsTerms().size();        
        int num_funcs=0;
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {
            num_funcs+=rhs_term_i.f_funcs.size();
        };
        (*f_ptr_i).rhs_ptrs_host.num_terms=new int[2];
        (*f_ptr_i).rhs_ptrs_host.num_funcs_1term=new int[num_terms];
        (*f_ptr_i).rhs_ptrs_host.prefactors=new double[num_terms];
        (*f_ptr_i).rhs_ptrs_host.schemes=new int[num_funcs];
        (*f_ptr_i).rhs_ptrs_host.f_func_ptrs=new double * [num_funcs];
        
        
        // Loop over terms on the RHS of each field.
        // cout <<"For field "<<(*f_ptr_i).name()<<": "<<endl;
        int i_func=0;
        int i_term=0;
        int num_terms_expl=0;
        // cout <<"RHS for field : "<<(*f_ptr_i).name()<<endl;
        // Add explicit terms
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {            
            // double* f_func_ptrs[10];
            int i_func_1term=0;
            if (rhs_term_i.scheme=="explicit") {
                num_terms_expl+=1;
                // (*f_ptr_i).rhs_ptrs_host.schemes[i_term]=1;
                
                // Evaluate each operator applied on field
                for (auto f_func_i : rhs_term_i.f_funcs) {
                    // Add this term to function list of this field. Each function appears only once.                    
                    // Check if function is already in the list
                    // cout <<"Check term: "<<f_func_i.f_operator<<":"<<(*f_func_i.field_ptr).name()<<endl;
                    int toAdd=1;
                    for (int i=0; i<(*f_func_i.field_ptr).num_f_funcs; i++) {
                        FFuncItem f_func_i1=(*f_func_i.field_ptr).f_funcs_rhs[i];
                        
                        if (f_func_i.f_operator == f_func_i1.f_operator) {
                            toAdd=0;
                        };
                        // cout <<f_func_i.f_operator<<"="<<f_func_i1.f_operator<<", toAdd="<<toAdd<<endl;
                    };
                    if (toAdd==1) {
                        // cout <<"Add term: "<<f_func_i.f_operator<<":"<<(*f_func_i.field_ptr).name()<<endl;
                        (*f_func_i.field_ptr).addFunctoRHS(f_func_i,device,FDMScheme);
                    };
                    
                    // Add field function pointer to rhs_ptrs_host.
                    (*f_ptr_i).rhs_ptrs_host.f_func_ptrs[i_func]=(*f_func_i.field_ptr).getFFuncPtr(f_func_i.f_operator);

                    i_func+=1;
                    i_func_1term+=1;
                };
                (*f_ptr_i).rhs_ptrs_host.num_funcs_1term[i_term]=i_func_1term;
                (*f_ptr_i).rhs_ptrs_host.prefactors[i_term]=rhs_term_i.prefactor;
                i_term+=1;
            };
        };

        int num_terms_impl=0;
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {            
            if (rhs_term_i.scheme=="semiImplicit") {
                int i_func_1term=0;
                num_terms_impl+=1;
                // (*f_ptr_i).rhs_ptrs_host.schemes[i_term]=-1;                
                // Evaluate each operator applied on field
                for (auto f_func_i : rhs_term_i.f_funcs) {                    
                
                    // Add this term to function list of this field. Each function appears only once.                    
                    int toAdd=1;
                    // Check if function is already in the list
                    for (int i=0; i<(*f_func_i.field_ptr).num_f_funcs; i++) {
                        FFuncItem f_func_i1=(*f_func_i.field_ptr).f_funcs_rhs[i];
                        if (f_func_i.f_operator == f_func_i1.f_operator) {
                            toAdd=0;
                        };
                    };
                    if (toAdd==1) {                            
                        // (*f_ptr_i).f_funcs_rhs.push_back(f_func_i);
                        (*f_func_i.field_ptr).addFunctoRHS(f_func_i,device,FDMScheme);
                    };

                    (*f_ptr_i).rhs_ptrs_host.f_func_ptrs[i_func]=(*f_func_i.field_ptr).getFFuncPtr(f_func_i.f_operator);

                    i_func+=1;
                    i_func_1term+=1;
                };
                (*f_ptr_i).rhs_ptrs_host.num_funcs_1term[i_term]=i_func_1term;
                (*f_ptr_i).rhs_ptrs_host.prefactors[i_term]=rhs_term_i.prefactor;
                i_term+=1;
            };
        };
        
        (*f_ptr_i).rhs_ptrs_host.num_terms[0]=num_terms_expl;
        (*f_ptr_i).rhs_ptrs_host.num_terms[1]=num_terms_impl;

        // Copy these values to device.
        if (device=="gpu") {                        
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.num_terms,2*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.num_funcs_1term,num_terms*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.schemes,num_terms*sizeof(int));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.prefactors,num_terms*sizeof(double));
            cudaMalloc(&(*f_ptr_i).rhs_ptrs_dev.f_func_ptrs,num_funcs*sizeof(double*));
            cudaMalloc(&(*f_ptr_i).f_funcs_dev,200*sizeof(double*));
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.num_terms, (*f_ptr_i).rhs_ptrs_host.num_terms, 2*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.num_funcs_1term, (*f_ptr_i).rhs_ptrs_host.num_funcs_1term, num_terms*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.schemes, (*f_ptr_i).rhs_ptrs_host.schemes, num_terms*sizeof(int),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.prefactors, (*f_ptr_i).rhs_ptrs_host.prefactors, num_terms*sizeof(double),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).rhs_ptrs_dev.f_func_ptrs, (*f_ptr_i).rhs_ptrs_host.f_func_ptrs, num_funcs*sizeof(double*),cudaMemcpyHostToDevice);
            cudaMemcpy((*f_ptr_i).f_funcs_dev, (*f_ptr_i).f_funcs_host, 200*sizeof(double*),cudaMemcpyHostToDevice);
        };

    };    
};


// ----------------------------------------------------------------------
void Evolver::getRHS(int i_field) {

    // Evaluate field functions for priority=0
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).priority()<=0) {
            evalFieldFuncs(f_ptr_i,i_field);            
        };
    };    

    // Get new field values for those with priority>0
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {        
        if ((*f_ptr_i).priority()>0) {
            updateRHS(f_ptr_i,i_field);
            // Apply periodic boundary condition.
            if (device=="cpu") {
                (*f_ptr_i).applyBounCondPeriCPU((*f_ptr_i).f[i_field]);
            } else if (device=="gpu"){
                (*f_ptr_i).applyBounCondPeriGPU((*f_ptr_i).f[i_field]);
            };

            if ((*f_ptr_i).specialty=="IncompressibleFlow.omega") {
                (*(*system_ptr).incomFlow_ptr).getVelocity(i_field);
            };
        };
    };

    // Evaluate field functions for priority>0 fields
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {        
        if ((*f_ptr_i).priority()>0) {
            evalFieldFuncs(f_ptr_i,i_field);
        };
    };

    // Update RHS of fields with priority=0
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).priority()==0) {
            updateRHS(f_ptr_i,i_field);
        };
    };

};


// ----------------------------------------------------------------------
void Evolver::evalFieldFuncs(Field* f_ptr_i, int i_field) {

    for (int i=0; i<(*f_ptr_i).num_f_funcs; i++) {
        FFuncItem f_func_i=(*f_ptr_i).f_funcs_rhs[i];
        
        if (device=="cpu") {            
            (*f_ptr_i).getFFuncCPU((*f_ptr_i).f_funcs_host[f_func_i.f_func_idx], i_field, f_func_i.f_func, f_func_i.f_func_args, "new");            
        } else if (device=="gpu") {
            (*f_ptr_i).getFFuncGPU((*f_ptr_i).f_funcs_host[f_func_i.f_func_idx], i_field, f_func_i.f_func, f_func_i.f_func_args, "new");
        };

    }; 
        
};


// -----------------------------------------------------------------------
void Evolver::updateRHS(Field* f_ptr_t, int i_field) {
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;    
    double* rhs_temp;
    double* lhs_temp;
    
    if ((*f_ptr_t).priority() == 0) {
        // Fields with time derivatives
        rhs_temp=(*f_ptr_t).rhs[i_field];
        lhs_temp=(*f_ptr_t).lhs[i_field];
    } else {
        // Fields without time derivative.
        if (device == "cpu") {
            rhs_temp=(*f_ptr_t).f[i_field];
            lhs_temp=(*f_ptr_t).f[i_field];
        } else if (device == "gpu") {
            rhs_temp=(*f_ptr_t).f[i_field];
            lhs_temp=(*f_ptr_t).f[i_field];
        };
    };    

    (*f_ptr_t).allocField<double>(rhs_temp, device);
    (*f_ptr_t).allocField<double>(lhs_temp, device);    
    
    if (device=="cpu") {
        updateRHSCoreCPU((*f_ptr_t).rhs_ptrs_host, rhs_temp, lhs_temp, Nx, Ny, Nbx, Nby);
    } else if (device=="gpu") {        
        updateRHSCoreGPU<<<Ny,Nx>>>((*f_ptr_t).rhs_ptrs_dev, rhs_temp, lhs_temp, Nx, Ny, Nbx, Nby);
    };
    
};

// ----------------------------------------------------------------------
void Evolver::updateRHSCoreCPU(rhsPtrs rhs_ptrs, double* rhs_temp, double* lhs_temp, int Nx, int Ny, int Nbx, int Nby) {
    
    double temp;
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            rhs_temp[idx]=0;
            lhs_temp[idx]=0;
            int i_func=0;
            for (int i_term=0; i_term<rhs_ptrs.num_terms[0]; i_term++) {
                temp=rhs_ptrs.prefactors[i_term];
                for (int i_func1=0; i_func1<rhs_ptrs.num_funcs_1term[i_term]; i_func1++) {
                    temp=temp*rhs_ptrs.f_func_ptrs[i_func][idx];
                    i_func+=1;
                };
                rhs_temp[idx]+=temp;
            };
            for (int i_term=rhs_ptrs.num_terms[0]; i_term<rhs_ptrs.num_terms[0]+rhs_ptrs.num_terms[1]; i_term++) {
                double temp=rhs_ptrs.prefactors[i_term];
                for (int i_func1=0; i_func1<rhs_ptrs.num_funcs_1term[i_term]-1; i_func1++) {
                    temp=temp*rhs_ptrs.f_func_ptrs[i_func][idx];
                    i_func+=1;
                };
                lhs_temp[idx]+=temp;
            };
        };
    };    
};


// ----------------------------------------------------------------------
void Evolver::fieldsUpdate(int i_f_new, int i_f_old, int i_df, double time_step_t) {
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).priority() ==0) {
            if (device=="cpu") {            
                fieldUpdateCPU(f_ptr_i,i_f_new,i_f_old,i_df,time_step_t);
                if ((*f_ptr_i).bounCond()=="periodic") {
                    (*f_ptr_i).applyBounCondPeriCPU((*f_ptr_i).f[i_f_new]);
                };
            } else if (device=="gpu") {
                fieldUpdateGPU(f_ptr_i,i_f_new,i_f_old,i_df,time_step_t);
                if ((*f_ptr_i).bounCond()=="periodic") {
                    (*f_ptr_i).applyBounCondPeriGPU((*f_ptr_i).f[i_f_new]);
                };
            };
            if ((*f_ptr_i).specialty=="IncompressibleFlow.omega") {
                (*(*system_ptr).incomFlow_ptr).getVelocity(i_f_new);
            };
        };        
    };
};

// ----------------------------------------------------------------------
void Evolver::fieldUpdateCPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t) {
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;
    
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            (*f_ptr_t).f[i_f_new][idx]=((*f_ptr_t).f[i_f_old][idx]+(*f_ptr_t).rhs[i_df][idx]*time_step_t)/(1+(*f_ptr_t).lhs[i_df][idx]*time_step_t);
        };
    };
};

// ----------------------------------------------------------------------
void Evolver::fieldUpdateGPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t) {
    
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;
    fieldUpdateGPUCore<<<Ny,Nx>>>((*f_ptr_t).f[i_f_new], (*f_ptr_t).f[i_f_old], (*f_ptr_t).rhs[i_df], (*f_ptr_t).lhs[i_df], time_step_t, Nx, Ny, Nbx, Nby);
};


// ----------------------------------------------------------------------
void Evolver::showProgress() {
  // Print progress.
    double progress=(time_now-time_start)/(time_stop-time_start);
    int barWidth = 50;
    time_now_wall = clock();
    double time_used_wall=double(time_now_wall-time_start_wall)/CLOCKS_PER_SEC;
  
    std::cout << "Progress: ";  
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %";
    if (time_now==0) {
        std::cout <<"\r";
    } else {
        std::cout << ".  " << floor(time_used_wall/progress*(1-progress)) << "s remains.\r";
    }

    std::cout.flush();
}


// ======================================================================

#endif
