#ifndef INCOMPFLOWOMEGASTOKESFRICFIELDCLASS_H
#define INCOMPFLOWOMEGASTOKESFRICFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "../LaplaceNFEqField/LaplaceNFEqFieldClass.cu"


using namespace std;


// ===============================================================
class IncompFlowOmegaStokesFricField :public LaplaceNFEqField {
    // This class is used to characterize the vorticity field,
    //   \omega, of an incompressible flow.
    
public:    
          
    // ===========================================================
    Field* ptr_vx;
    Field* ptr_vy;
    Field* ptr_phi;                  // stream function
    cufftDoubleComplex* phi_complex;
    double* poisson_k2_host;
    double* poisson_k2_dev;
    cufftHandle cufftPlan;
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    IncompFlowOmegaStokesFricField () {};
    IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t);
    IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    IncompFlowOmegaStokesFricField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Get velocity field
    void initPoissonSolver();
    void getVelocity(int i_field);
    void getIncompFlowVCoreCPU(int i_field, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args);    
    
    // ===========================================================
};



#endif
