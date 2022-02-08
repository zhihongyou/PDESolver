#ifndef EVOLVERCLASS_CPP
#define EVOLVERCLASS_CPP

#include <iostream> 
#include <vector>
// #include "../mesh/meshclass.h"
// #include "../field/fieldclass.h"
#include "evolverclass.h"

using namespace std;


// ======================================================================

// ----------------------------------------------------------------------
// Update system over time. This is the global control.
void Evolver::run() {
    cout <<"Start running simulation ..." <<endl;
    time_start_wall = clock();
    
    // Export initial configuration.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).expo_data =="on") {
            (*f_ptr_i).export_conf("0",0);
            // (*f_ptr_i).getLaplace(0,"new");
            // (*f_ptr_i).export_conf_any(
            //             (*f_ptr_i).laplace,
            //             "laplace_"+(*f_ptr_i).name,"0");
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
                if ((*f_ptr_i).expo_data =="on") {
                    (*f_ptr_i).export_conf(str_t,0);
                    // (*f_ptr_i).export_conf_any(
                    //     (*f_ptr_i).laplace,
                    //     "laplace_"+(*f_ptr_i).name,str_t);
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
};

// ----------------------------------------------------------------------
void Evolver::fieldsUpdate(int i_f_new, int i_f_old, int i_df) {
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        if ((*f_ptr_i).priority ==0) {
            fieldUpdate(f_ptr_i,i_f_new,i_f_old,i_df,time_step);
            if ((*f_ptr_i).boun_cond=="periodic") {
                (*f_ptr_i).applyBounCondPeri(i_f_new);
            };
        };        
    };
};

// ----------------------------------------------------------------------
void Evolver::fieldUpdate(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t) {
    int Nx=(*(*f_ptr_t).mesh_ptr).grid_number.x;
    int Ny=(*(*f_ptr_t).mesh_ptr).grid_number.y;
    int Nbx=(*(*f_ptr_t).mesh_ptr).grid_number_boun.x;
    int Nby=(*(*f_ptr_t).mesh_ptr).grid_number_boun.y;
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nby)*(Ny+2*Nby)+i+Nbx;
            (*f_ptr_t).f_cpu[i_f_new][idx]=(*f_ptr_t).f_cpu[i_f_old][idx]+
                (*f_ptr_t).f_rhs[i_df][idx]*time_step_t;
        };
    };    
};


// ----------------------------------------------------------------------
void Evolver::getRHSs(int i_field) {
    // Loop over fields.
    for (auto f_ptr_i : (*system_ptr).field_ptrs ) {
        clearRHS(f_ptr_i,i_field);
        // Loop over terms on the RHS of each field.
        // cout <<"For field "<<(*f_ptr_i).name<<": "<<endl;        
        for (auto rhs_term_i : (*f_ptr_i).rhs_terms) {            
            double* f_func_ptrs[10];
            int N_funcs=0;
            // cout <<"term "<<N_funcs<<" with operator "<<rhs_term_i.rhs_operator<<" and number "<<rhs_term_i.numbers[0]<<" :";
            // Evaluate each operator applied on field
            for (auto f_func_i : rhs_term_i.f_function) {
                // cout<<" "<<f_func_i.f_operator<<"("<<(*f_func_i.field_ptr).name <<"), ";
                evalOperator(f_ptr_i,f_func_i,f_func_ptrs[N_funcs],i_field);                
                N_funcs+=1;
            };
            // cout <<endl;
            // Add one RHS term to field increments
            addRHSTerm(f_ptr_i,i_field,rhs_term_i,f_func_ptrs,N_funcs);
        };
        if ((*f_ptr_i).priority>0 && (*f_ptr_i).boun_cond=="periodic") {
            (*f_ptr_i).applyBounCondPeri(i_field);
        };
    };
};

//  -----------------------------------------------------------------------
void Evolver::clearRHS(Field* f_ptr_t, int i_field) {
    double* f_temp;
    if ((*f_ptr_t).priority==0) {
        f_temp=(*f_ptr_t).f_rhs[i_field];
    } else {
        f_temp=(*f_ptr_t).f_cpu[i_field];
    };
    int Nx=(*(*f_ptr_t).mesh_ptr).grid_number.x;
    int Ny=(*(*f_ptr_t).mesh_ptr).grid_number.y;
    int Nbx=(*(*f_ptr_t).mesh_ptr).grid_number_boun.x;
    int Nby=(*(*f_ptr_t).mesh_ptr).grid_number_boun.y;
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nbx)*(Nx+2*Nbx)+i+Nbx;
            f_temp[idx]=0;
        };                
    };
};

// -----------------------------------------------------------------------
void Evolver::evalOperator(Field* f_ptr_t, field_function f_func_t, double*& f_func_ptr, int i_field) {
    if (f_func_t.f_operator=="1") {
        if (device=="cpu") {
            f_func_ptr=(*f_func_t.field_ptr).f_cpu[i_field];
        } else {
            f_func_ptr=(*f_func_t.field_ptr).f_gpu[i_field];
        }
    };
    if (f_func_t.f_operator=="laplace") {
        if (device=="cpu") {
            f_func_ptr=(*f_func_t.field_ptr).getLaplace(i_field,"new");
        } else {
            f_func_ptr=(*f_func_t.field_ptr).getLaplace(i_field,"new");
        }
    };
    if (f_func_t.f_operator=="bi_laplace") {
        if (device=="cpu") {
            f_func_ptr=(*f_func_t.field_ptr).getBiLaplace(i_field,"new");
        } else {
            f_func_ptr=(*f_func_t.field_ptr).getBiLaplace(i_field,"new");
        }
    };
    
};


// -----------------------------------------------------------------------
void Evolver::addRHSTerm(Field* f_ptr_t, int i_field, rhs_term rhs_term_t, double** f_func_ptrs, int N_funcs) {
    double* f_rhs_temp_ptr;
    if ((*f_ptr_t).priority == 0) {
        if (rhs_term_t.scheme=="explicit") {
            f_rhs_temp_ptr=(*f_ptr_t).f_rhs[i_field];
        } else {
            f_rhs_temp_ptr=(*f_ptr_t).f_lhs[i_field];
        };
    } else {
        if (device == "cpu") {
            f_rhs_temp_ptr=(*f_ptr_t).f_cpu[i_field];
        } else {
            f_rhs_temp_ptr=(*f_ptr_t).f_gpu[i_field];
        };
    };
    if (f_rhs_temp_ptr == NULL) {
        f_rhs_temp_ptr=new double[(*(*f_ptr_t).mesh_ptr).getGridNumberAll()];
    };
    
    int Nx=(*(*f_ptr_t).mesh_ptr).grid_number.x;
    int Ny=(*(*f_ptr_t).mesh_ptr).grid_number.y;
    int Nbx=(*(*f_ptr_t).mesh_ptr).grid_number_boun.x;
    int Nby=(*(*f_ptr_t).mesh_ptr).grid_number_boun.y;
    for (int j=0; j<Ny;j++) {
        for (int i=0; i<Nx; i++) {        
            int idx=(j+Nbx)*(Nx+2*Nbx)+i+Nbx;
            if (rhs_term_t.rhs_operator == "+") {
                for (int i_func=0; i_func<N_funcs;i_func++) {
                    f_rhs_temp_ptr[idx]+=f_func_ptrs[i_func][idx];
                };
                for (auto n_i : rhs_term_t.numbers) {
                    f_rhs_temp_ptr[idx]+=n_i;
                };
            } else if (rhs_term_t.rhs_operator == "*") {
                double temp=1;
                for (int i_func=0; i_func<N_funcs;i_func++) {
                    temp=temp*f_func_ptrs[i_func][idx];
                };
                for (auto n_i : rhs_term_t.numbers) {                    
                    temp=temp*n_i;
                };
                f_rhs_temp_ptr[idx]+=temp;
            };
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
