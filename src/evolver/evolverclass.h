#ifndef EVOLVERCLASS_H
#define EVOLVERCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <time.h>
// #include "../utility/data_structs.h"
// #include "../mesh/meshclass.h"
// #include "../mesh/meshclass.cpp"
// #include "../field/fieldclass.h"
// #include "../field/fieldclass.cpp"
// #include "../system/systemclass.h"
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

    // Simulation time:
    // Simulation starts at timeStart and stops at timeStop, with an
    //   initial time step timeStep.
    // Fields are exported every timeExport.
    double time_now;
    double time_start=0, time_stop=1, time_step=0.001, time_export=0.1;
    clock_t time_start_wall, time_now_wall, time_stop_wall;
    
    // Define Constructors ()
    // Evolver();
    Evolver(System * system_ptr1, double time_start_t=0, double time_stop_t=1, double time_step_t=0.001, double time_export_t=0.1, std::string scheme_t="EulerForward") {
        system_ptr=system_ptr1;
        time_start=time_start_t;
        time_stop=time_stop_t;
        time_step=time_step_t;
        time_export=time_export_t;
        scheme=scheme_t;
    };

    
    // ==================================================================
    // Methods
    // void getRhsOneOperate (Field* f_t, rhs_operates operate_t);
    // void getRHSs (Field* f_t, std::vector<rhs_operates> phi_rhs);
    void clearRHS(Field* f_ptr_t, int i_field);
    void getRHSs(int i_f_copy);
    void getRHSField(Field* field_ptr_t, int i_field);
    void addRHSTerm(Field* f_ptr_t, int i_field, rhs_term rhs_term_t, double** f_func_ptrs, int N_funcs);
    void evalOperator(Field* f_ptr_t, field_function f_func_t, double*& f_func_ptr, int i_field);
    void fieldsUpdate(int i_f_new, int i_f_old, int i_df);
    void fieldUpdate(Field* f_ptr_t, int i_f_new, int i_f_old, int i_df, double time_step_t);
    void EulerForward();
    void run();
    void showProgress();
    
};

// =======================================================================
// Other methods useful for evolver.


#endif
