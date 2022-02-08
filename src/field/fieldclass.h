#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../mesh/meshclass.h"
#include "../mesh/meshclass.cpp"

using namespace std; 

struct rhs_term;

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
    std::string f_equation = "";
    Mesh* mesh_ptr=NULL;
    std::vector<rhs_term> rhs_terms;
    
    // Field data
    double* f_cpu[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_gpu[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_rhs[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_lhs[5]={NULL,NULL,NULL,NULL,NULL};
    
    // Operators
    double * d1xf=NULL;
    double * d1yf=NULL;
    double * d1zf=NULL;
    double * d1x1yf=NULL;
    double * d1x1zf=NULL;
    double * d1y1zf=NULL;
    double * d2xf=NULL;
    double * d2yf=NULL;
    double * d2zf=NULL;
    double * d2x2yf=NULL;
    double * d2x2zf=NULL;
    double * d2y2zf=NULL;
    double * one_over_f=NULL;
    double * laplace=NULL;
    double * bi_laplace=NULL;

    
    // ===================================================================
    // Methods
    
    // ===================================================================
    // Constructor
    Field (Mesh* mesh_ptr_t, std::string name_t, int rank_t, int priority_t, std::string boun_cond_t, std::string init_cond_t, std::string expo_data_t);


    // ------------------------------------------------------------------
    // Get rhs
    void getRHS(int i_f_copy);
    // -------------------------------------------------------------------
    // Field initialization
    // Constant field
    void initFieldConst(double f_value);
    // Random field with uniform distribution
    void initFieldRandUnif(double f_mean, double f_var);
    // Random field with normal distribution
    void initFieldRandNorm(double f_mean, double f_var);
    // Gaussian profile in space f(r)=gaus_amplitude*exp(-(r-r_center)^2/r_decay^2).
    void initFieldGaus(double r_center, double r_decay, double gaus_amplitude);
    void initFieldSin(double sin_amplitude, int sin_period, double sin_phase);

    // ===================================================================
    // Differential operators
    double* getLaplace(int i_field,std::string method);
    double* getBiLaplace(int i_field,std::string method);

    // -------------------------------------------------------------------
    // Field boundary condition
    void applyBounCondPeri(int i_field);
    void applyBounCondPeriAny(double* f_t);

    // -------------------------------------------------------------------
    // Export field
    void export_conf(string str_t, int include_boun);
    void export_conf_any(double* f_t, string f_name, string str_t, int include_boun);
    void getEqn();
    
    
    
    // ===================================================================
};

#endif
