#ifndef LAPLACENFEQFIELDCLASS_H
#define LAPLACENFEQFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// ===============================================================
class LaplaceNFEqField :public Field {
    // This class is used to characterize a filed governed by the
    //   Poisson-like equation (\nabla^2)^n_laplace \phi = f,
    //   where f is specified by the RHS of field \phi.
    // When n_laplace=1, the equation is the Poisson equation,
    //   and n_laplace=2 gives the bi-harmonic equation.
    
public:    
       
    LaplaceNFEqSolver LaplNFEqSolver;
    int max_power=1;
    double* prefactors=new double[10];    
    
    // ===========================================================
    // Constructor
    LaplaceNFEqField ();
    LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t);
    LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    LaplaceNFEqField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Methods
    void initFieldAddi ();
    void setLaplaceNFEq (int max_power_t, double* prefactors_t);
    void solveLaplaceNFEq (int i_field);
    // ===========================================================
};



#endif
