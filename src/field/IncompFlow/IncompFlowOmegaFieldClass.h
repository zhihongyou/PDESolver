#ifndef INCOMPFLOWOMEGAFIELDCLASS_H
#define INCOMPFLOWOMEGAFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "../../utility/LaplaceNFEqSolver/LaplaceNFEqSolverClass.cu"


using namespace std; 


// ===============================================================
class IncompFlowOmegaField :public Field {
    // This class is used to characterize the vorticity field,
    //   \omega, of an incompressible flow.
    
public:
          
    // ===========================================================
    Field* ptr_vx;
    Field* ptr_vy;
    Field* ptr_phi;                  // stream function
    // A Poisson solver to get stream function from vorticity
    LaplaceNFEqSolver LaplSolverW2Phi;
    // Solver to get vorticity
    LaplaceNFEqSolver LaplSolverW;
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    IncompFlowOmegaField () {};
    IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t);
    IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    IncompFlowOmegaField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Get velocity field
    void initFieldAddi ();
    void getVelocity(int i_field);
    void getIncompFlowVCoreCPU(int i_field, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args);
    void setStokesFricSolver (double Gamma, double eta);
    
    // ===========================================================
};



#endif
