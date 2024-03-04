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


// ==============================================================

#endif
