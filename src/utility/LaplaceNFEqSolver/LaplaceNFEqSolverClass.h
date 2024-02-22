#ifndef LAPLACENFEQSOLVERCLASS_H
#define LAPLACENFEQSOLVERCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../../mesh/meshclass.cpp"
#include "LaplaceNFEqSolverClassGPU.cu"


using namespace std; 


// ===============================================================
class LaplaceNFEqSolver {
    // This class is used to characterize a filed governed by the
    //   Poisson-like equation (\nabla^2)^n_laplace \phi = f,
    //   where f is specified by the RHS of field \phi.
    // When n_laplace=1, the equation is the Poisson equation,
    //   and n_laplace=2 gives the bi-harmonic equation.
    
public:    

    Mesh* mesh_ptr;
    string name="LaplaceNFEqSolver";
    
    cufftDoubleComplex* phi_complex;
    double* k2s_host=NULL;
    double* k2s_dev=NULL;
    double* prefactors=new double[10];
    int max_power=1;
    cufftHandle cufftPlan;
    
    // ===========================================================
    // Constructor
    LaplaceNFEqSolver ();
    LaplaceNFEqSolver (Mesh* mesh_ptr_t);
    LaplaceNFEqSolver (Mesh* mesh_ptr_t, string name_t);
    
    // Methods
    void initLaplaceNFEqSolver(Mesh* mesh_ptr_t);
    void setLaplaceNFEqSolver(int max_power_t, double* prefactors_t);
    void setk2s ();
    void solveLaplaceNFEq (double* phi, double* f);
    // ===========================================================
};



#endif
