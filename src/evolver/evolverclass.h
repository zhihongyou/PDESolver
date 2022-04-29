#ifndef EVOLVERCLASS_H
#define EVOLVERCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <time.h>
#include "../system/systemclass.cpp"


using namespace std;


//................Class .................................................

class Evolver{

public:

    // Which device to use to do computationally expensive parts.
    // "cpu" (default) or "gpu".
    std::string device="cpu";
    
    // Explicit or semi-implicit schemes for potential terms.
    // "explicit" (default) or "semiImplicit".
    std::string scheme="EulerForward";

    System* system_ptr;

    // Number of copies of fields. Determined by the scheme.
    int num_field_copy;

    // Simulation time:
    // Simulation starts at timeStart and stops at timeStop, with an
    //   initial time step timeStep.
    // Fields are exported every timeExport.
    double time_now;
    double time_start=0, time_stop=1, time_step=0.001, time_export=0.1;
    clock_t time_start_wall, time_now_wall, time_stop_wall;
    
    // Define Constructors ()
    // Evolver();
    Evolver(System * system_ptr1, double time_start_t=0, double time_stop_t=1, double time_step_t=0.001, double time_export_t=0.1, std::string scheme_t="EulerForward", string device_t="cpu") {
        system_ptr=system_ptr1;
        time_start=time_start_t;
        time_stop=time_stop_t;
        time_step=time_step_t;
        time_export=time_export_t;
        scheme=scheme_t;
        device=device_t;
        if (scheme=="EulerForward") {
            num_field_copy=1;
        } else if (scheme=="PredictorCorrector") {
            num_field_copy=2;
        };
    };

    
    // ==================================================================
    // Methods
    // void getRhsOneOperate (Field* f_t, rhs_operates operate_t);
    // void getRHSs (Field* f_t, std::vector<rhs_operates> phi_rhs);
    void initEvolver();
    void initRHSs();
    void initFields();
    void allocateRHS(Field* f_ptr_t, int i_field);
    void getRHS(int i_field);
    void evalFieldFuncs(Field* f_ptr_t, int i_field);    
    void updateRHS(Field* f_ptr_t, int i_field);
    void updateRHSCoreCPU(rhsPtrs rhs_ptrs, double* rhs_temp, double* lhs_temp, int Nx, int Ny, int Nbx, int Nby);
    void fieldsUpdate(int i_f_new, int i_f_old, int i_df);
    void fieldUpdateCPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t);
    void fieldUpdateGPU(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t);    
    void EulerForward();
    void run();
    void showProgress();
    
};

// =======================================================================
// Other methods useful for evolver.


#endif
