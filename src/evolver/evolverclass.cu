#ifndef EVOLVERCLASS_CPP
#define EVOLVERCLASS_CPP

#include <iostream> 
#include <vector>
#include "evolverclass.h"
#include "evolverclassGPU.cu"

using namespace std;


// ======================================================================

// ----------------------------------------------------------------------
// Update system over time. This is the global control.
void Evolver::run() {
    cout <<"Start running simulation ..." <<endl;
    time_start_wall = clock();
    
    // Export initial configuration.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {                
        if (device=="gpu") {
            (*f_ptr_i).updateMainFieldDev();            
        };       
                
        if ((*f_ptr_i).expoData() =="on") {
            (*f_ptr_i).export_conf("0",device,1);
        };
    };

    // Running through time
    for (time_now=time_start; time_now<=time_stop+0.1*time_step; time_now+=time_step) {
        
        // Move a step forward.
        if (scheme=="EulerForward") {
            EulerForward();
        };

        // Reach export time.
        if ((int) (time_now/time_export) > (int) ((time_now-time_step)/time_export)) {
            // Export field configuration.
            int te=floor(time_now/time_export);    
            string str_t=to_string(te);
            for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
                if ((*f_ptr_i).expoData() =="on") {
                    (*f_ptr_i).export_conf(str_t,device,1);                   
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
// Euler forward scheme to evolve over time.
void Evolver::EulerForward() {
    
    getRHSs(0);        
    fieldsUpdate(0,0,0);
    // for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
    //     (*f_ptr_i).export_conf("0",device,1);
    // };
};

// ----------------------------------------------------------------------
void Evolver::fieldsUpdate(int i_f_new, int i_f_old, int i_df) {
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).priority() ==0) {
            if (device=="cpu") {
                fieldUpdateCPU(f_ptr_i,i_f_new,i_f_old,i_df,time_step);                
                if ((*f_ptr_i).bounCond()=="periodic") {
                    (*f_ptr_i).applyBounCondPeriCPU(i_f_new);
                };
            } else if (device=="gpu") {
                fieldUpdateGPU(f_ptr_i,i_f_new,i_f_old,i_df,time_step);
                if ((*f_ptr_i).bounCond()=="periodic") {
                    (*f_ptr_i).applyBounCondPeriGPU(i_f_new);
                };
                
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
            (*f_ptr_t).f_host[i_f_new][idx]=(*f_ptr_t).f_host[i_f_old][idx]+
                (*f_ptr_t).f_rhs[i_df][idx]*time_step_t;
        };
    };
};

// ----------------------------------------------------------------------
void Evolver::fieldUpdateGPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t) {
    
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;
    fieldUpdateGPUCore<<<Ny,Nx>>>((*f_ptr_t).f_dev[i_f_new], (*f_ptr_t).f_dev[i_f_old], (*f_ptr_t).f_rhs[i_df], time_step_t, Nx, Ny, Nbx, Nby);
};


// ----------------------------------------------------------------------
void Evolver::getRHSs(int i_field) {
    
    // Loop over fields.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        
        clearRHS(f_ptr_i,i_field);
        
        // Loop over terms on the RHS of each field.
        // cout <<"For field "<<(*f_ptr_i).name()<<": "<<endl;        
        for (auto rhs_term_i : (*f_ptr_i).rhsTerms()) {            
            // double* f_func_ptrs[10];
            double** f_func_ptrs=new double*[10];
            int N_funcs=0;
            // Evaluate each operator applied on field
            for (auto f_func_i : rhs_term_i.f_function) {
                // cout<<""<<f_func_i.f_operator<<"("<<(*f_func_i.field_ptr).name() <<"),";
                evalOperator(f_ptr_i,f_func_i,f_func_ptrs[N_funcs],i_field);                
                N_funcs+=1;
            };
            
            addRHSTerm(f_ptr_i,i_field,rhs_term_i,f_func_ptrs,N_funcs);
                       
        };

        if ((*f_ptr_i).priority()>0 && (*f_ptr_i).bounCond()=="periodic") {
            if (device=="cpu") {
                (*f_ptr_i).applyBounCondPeriCPU(i_field);
            } else if (device=="gpu"){
                (*f_ptr_i).applyBounCondPeriGPU(i_field);
            };
        };        
    };
};

// ----------------------------------------------------------------------
void Evolver::clearRHS(Field* f_ptr_t, int i_field) {
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;
    if ((*f_ptr_t).priority()==0) {
        if (device=="cpu") {
            if ((*f_ptr_t).f_rhs[i_field] == NULL) {
                (*f_ptr_t).f_rhs[i_field]=new double[(*f_ptr_t).gridNumberAll()];
            };
            (*f_ptr_t).setFieldConstCPU((*f_ptr_t).f_rhs[i_field], 0, Nx, Ny, Nbx, Nby);
        } else if (device=="gpu") {
            if ((*f_ptr_t).f_rhs[i_field] == NULL) {
                cudaMalloc(&(*f_ptr_t).f_rhs[i_field],
                (*f_ptr_t).gridNumberAll()*sizeof(double));
            };
            (*f_ptr_t).setFieldConstGPU((*f_ptr_t).f_rhs[i_field], 0, Nx, Ny, Nbx, Nby);
        };
        
    } else {
        if (device=="cpu") {
            if ((*f_ptr_t).f_host[i_field] == NULL) {
                (*f_ptr_t).f_host[i_field]=new double[(*f_ptr_t).gridNumberAll()];
            };
            (*f_ptr_t).setFieldConstCPU((*f_ptr_t).f_host[i_field], 0, Nx, Ny, Nbx, Nby);
        } else if (device=="gpu") {
            if ((*f_ptr_t).f_dev[i_field] == NULL) {
                cudaMalloc(&(*f_ptr_t).f_dev[i_field],
                (*f_ptr_t).gridNumberAll()*sizeof(double));
            };
            (*f_ptr_t).setFieldConstGPU((*f_ptr_t).f_dev[i_field], 0, Nx, Ny, Nbx, Nby);
        };
        
    };

};


// -----------------------------------------------------------------------
void Evolver::evalOperator(Field* f_ptr_t, field_function f_func_t, double*& f_func_ptr, int i_field) {
    if (f_func_t.f_operator=="1") {
        if (device=="cpu") {
            f_func_ptr=(*f_func_t.field_ptr).f_host[i_field];
        } else {
            f_func_ptr=(*f_func_t.field_ptr).f_dev[i_field];
        }
    };
    if (f_func_t.f_operator=="laplace") {
        if (device=="cpu") {
            f_func_ptr=(*f_func_t.field_ptr).getLaplaceCPU(i_field,"new");
        } else if (device=="gpu") {
            f_func_ptr=(*f_func_t.field_ptr).getLaplaceGPU(i_field,"new");
        }
    };    
};


// -----------------------------------------------------------------------
void Evolver::addRHSTerm(Field* f_ptr_t, int i_field, rhs_term rhs_term_t, double** f_func_ptrs, int N_funcs) {
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;    
    double* f_rhs_temp_ptr;
    
    if ((*f_ptr_t).priority() == 0) {
        if (rhs_term_t.scheme=="explicit") {
            f_rhs_temp_ptr=(*f_ptr_t).f_rhs[i_field];
        } else {
            f_rhs_temp_ptr=(*f_ptr_t).f_lhs[i_field];
        };
    } else {
        if (device == "cpu") {
            f_rhs_temp_ptr=(*f_ptr_t).f_host[i_field];
        } else if (device == "gpu") {
            f_rhs_temp_ptr=(*f_ptr_t).f_dev[i_field];
        };
    };    
    
    if (f_rhs_temp_ptr == NULL) {
        if (device=="cpu") {
            f_rhs_temp_ptr=new double[(*f_ptr_t).gridNumberAll()];
            (*f_ptr_t).setFieldConstCPU(f_rhs_temp_ptr, 0, Nx, Ny, Nbx, Nby);
        } else if (device=="gpu") {
            cudaMalloc(&f_rhs_temp_ptr,
            (*f_ptr_t).gridNumberAll()*sizeof(double));
            (*f_ptr_t).setFieldConstGPU(f_rhs_temp_ptr, 0, Nx, Ny, Nbx, Nby);
        };
    };
    
    if (device=="cpu") {        
        addRHSTermCPU(rhs_term_t, f_rhs_temp_ptr, f_func_ptrs, N_funcs, Nx, Ny, Nbx, Nby);
    } else if (device=="gpu") {
        double** f_func_ptrs_dev;
        cudaMalloc(&f_func_ptrs_dev,10*sizeof(double*));
        cudaMemcpy(f_func_ptrs_dev, f_func_ptrs, 10*sizeof(double*),cudaMemcpyHostToDevice);
        addRHSTermGPU<<<Ny,Nx>>>(rhs_term_t, f_rhs_temp_ptr, f_func_ptrs_dev, N_funcs, Nx, Ny, Nbx, Nby);
    };
    
};

// ----------------------------------------------------------------------
void Evolver::addRHSTermCPU(rhs_term rhs_term_t, double* f_rhs_temp_ptr, double** f_func_ptrs, int N_funcs, int Nx, int Ny, int Nbx, int Nby) {
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double temp=rhs_term_t.prefactor;
            for (int i_func=0; i_func<N_funcs;i_func++) {
                temp=temp*f_func_ptrs[i_func][idx];
            };
            f_rhs_temp_ptr[idx]+=temp;
            
        };
    };    
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
