#ifndef INCOMPFLOWOMEGAFIELDCLASS_H
#define INCOMPFLOWOMEGAFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// ===============================================================
class IncompFlowOmegaField :public Field {
    // This class is used to characterize the vorticity field,
    //   \omega, of an incompressible flow.
    // The Eq for \omega can take one of the three forms:
    //   1. The standard NSE: dt\omega=RHS or \omega=RHS.
    //   2. The Stokes Eq: -\eta*\nabla=RHS.
    //   3. The frictional Stokes Eq:
    //          gammma*\omega-eta*Laplace(\omega)=RHS.
    // The default setting uses the 1st form.
    // To enable the 2nd&3rd form, simply set omegaEqType to 2,
    //    it will activate additional steps to calculate
    //    \omega field.
    
public:
          
    // ===========================================================
    Field* ptr_vx;
    Field* ptr_vy;
    Field* ptr_phi;                  // stream function
    // A Poisson solver to get stream function from vorticity
    LaplaceNFEqSolver LaplSolverW2Phi;
    // Solver to get vorticity
    LaplaceNFEqSolver LaplSolverW;
    int omegaEqType=1;
    // Substrate friction
    double gamma=0;
    // Viscosity, only applicable when omegaEqMode > 1
    double eta=1;
    
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
    void initOmegaField ();
    void setOmegaEq (int omegaEqType_t, double gamma_t, double eta_t);
    void getOmegaAddi (int i_field);
    void postProcessing(int i_field);
    void getIncompFlowVCoreCPU(int i_field, FFuncType d1x, FFuncType d1y, FFuncArgs f_func_args);
    
    // ===========================================================
};



#endif
