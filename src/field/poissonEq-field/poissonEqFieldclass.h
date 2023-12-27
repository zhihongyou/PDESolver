#ifndef POISSONEQFIELDCLASS_H
#define POISSONEQFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// ===============================================================
class poissonEqField :public Field{
    // This class is used to characterize a filed governed by the
    //   Poisson equation \nabla^2 \phi = f,
    //   where f is specified by the RHS of field \phi.
    
public:    
       
    cufftDoubleComplex* phi_complex;
    double* poisson_k2_host;
    double* poisson_k2_dev;
    cufftHandle cufftPlan;
    
    // ===========================================================
    // Constructor
    poissonEqField (Mesh* mesh_ptr_t, string name_t);
    poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    poissonEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Methods
    void initFieldAddi ();
    void initPoissonSolver();    
    void solvePoissonEq (int i_field);
    void getRHSAddi (int i_field);
    // ===========================================================
};



#endif
