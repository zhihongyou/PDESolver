#ifndef EVOLVERCLASS_INTEGRATOR_CU
#define EVOLVERCLASS_INTEGRATOR_CU

// --------------------------------------------------------------
// Euler forward scheme to evolve over time.
void Evolver::EulerForward() {        

    getRHS(0);
    // Update fields with zero priority.
    fieldsUpdate(0, 0, 0, time_step);
};


// --------------------------------------------------------------
// Euler forward scheme to evolve over time.
void Evolver::RK4() {    
    getRHS(0);
    fieldsUpdate(1, 0, 0, 0.5*time_step);
    getRHS(1);
    fieldsUpdate(2, 0, 1, 0.5*time_step);
    getRHS(2);
    fieldsUpdate(3, 0, 2, time_step);
    getRHS(3);
    fieldsUpdate(0, 0, 0, 1.0/6.0*time_step);
    fieldsUpdate(0, 0, 1, 1.0/3.0*time_step);
    fieldsUpdate(0, 0, 2, 1.0/3.0*time_step);
    fieldsUpdate(0, 0, 3, 1.0/6.0*time_step);
};


// --------------------------------------------------------------
void Evolver::fieldsUpdate(int i_f_new, int i_f_old, int i_df, double time_step_t) {
    // Update fields with priority=0
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
            // Get velocity
            if ((*f_ptr_i).specialty=="IncompFlowOmega") {
                IncompFlowOmegaField* f_ptr_temp = (IncompFlowOmegaField*) f_ptr_i;
                (*f_ptr_temp).getVelocity(i_f_new);
            };
        };        
    };
};


// --------------------------------------------------------------
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


// --------------------------------------------------------------
void Evolver::fieldUpdateGPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t) {
    
    int Nx=(*f_ptr_t).gridNumber().x;
    int Ny=(*f_ptr_t).gridNumber().y;
    int Nbx=(*f_ptr_t).gridNumberBoun().x;
    int Nby=(*f_ptr_t).gridNumberBoun().y;
    fieldUpdateGPUCore<<<Ny,Nx>>>((*f_ptr_t).f[i_f_new], (*f_ptr_t).f[i_f_old], (*f_ptr_t).rhs[i_df], (*f_ptr_t).lhs[i_df], time_step_t, Nx, Ny, Nbx, Nby);
};


// ==============================================================

#endif
