#ifndef INCOMPRESSIBLEFLOWCLASS_H
#define INCOMPRESSIBLEFLOWCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// =======================================================================
class IncompressibleFlow {    
    
    public:

    string name;
    int priority=0;
    string boun_cond = "periodic";
    string init_cond = "sin";
    int rank=0;
    string location = "both";
    string expo_data = "on";
    string equation = "";
    Mesh* mesh_ptr=NULL;
    
    Field vx;
    Field vy;
    Field omega;                // omega=2*vorticity
    Field phi;                  // stream function
    cufftDoubleComplex* phi_complex;

    string FDMScheme;

    double* poisson_k2_host;
    double* poisson_k2_dev;

    cufftHandle cufftPlan;
    
    // ===================================================================
    // Methods
    
    // ===================================================================
    // Constructor
    IncompressibleFlow (Mesh* mesh_ptr_t, string name_t);
    IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t);
    IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    IncompressibleFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);

    // Fields
    void initFields ();
    void setFieldProperties (Field* field_ptr, string field_name, int priority_t);
    void setFieldValues (string init_cond_t);
    void setVValues();
    void setOmegaValues ();

    // Get velocity field
    void initPoissonSolver();
    void getVelocity (int i_field);
    void getVCoreCPU(double* phi, double* vx, double* vy, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args);
    
        
    // ===================================================================
};



#endif
