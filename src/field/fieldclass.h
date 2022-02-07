#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../mesh/meshclass.h"
#include "../mesh/meshclass.cpp"

using namespace std; 


// =======================================================================
class Field {

    
    public:

    // -------------------------------------------------------------------
    // Field properties
    std::string name = "fieldTemp";
    int rank=0;
    int priority=1;
    std::string location = "both";
    std::string boun_cond = "periodic";
    std::string init_cond = "none";
    std::string expo_data = "on";
    Field * field_ptr;
    
    // Field data
    double * f_cpu;
    double * f_gpu;
    double * f_rhs;
    double * f_copy1;
    double * f_copy2;
    double * f_copy3;
    double * f_copy4;
    double * f_copy5;
    double * f_rhs_copy1;
    double * f_rhs_copy2;
    double * f_rhs_copy3;
    double * f_rhs_copy4;
    double * f_rhs_copy5;

    
    // ===================================================================
    // Methods
    // -------------------------------------------------------------------
    // Field initialization
    // Constant field
    void initFieldConst(Mesh mesh_t, double*& f_t, double f_value);
    // Random field with uniform distribution
    void initFieldRandUnif(Mesh mesh_t, double*& f_t, double f_mean, double f_var);
    // Random field with normal distribution
    void initFieldRandNorm(Mesh mesh_t, double*& f_t, double f_mean, double f_var);
    // Gaussian profile in space f(r)=gaus_amplitude*exp(-(r-r_center)^2/r_decay^2).
    void initFieldGaus(Mesh mesh_t, double*& f_t, double r_center, double r_decay, double gaus_amplitude);

    // -------------------------------------------------------------------
    // Field boundary condition
    void applyBounCond(Mesh mesh_t, double*& f_t);

    // -------------------------------------------------------------------
    // Export field
    void exportField(Mesh mesh_t, double*& f_t);

    
    // ===================================================================
    // Constructor
    Field (Mesh f_mesh, std::string f_name, int f_rank, int f_priority, std::string f_boun_cond, std::string f_init_cond, std::string f_expo_data) {
        name=f_name;
        rank=f_rank;
        priority=f_priority;
        boun_cond=f_boun_cond;
        init_cond=f_init_cond;
        expo_data=f_expo_data;
        initFieldConst(f_mesh,f_cpu,2);
    };
    
    // ===================================================================
};

#endif
