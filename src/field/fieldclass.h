#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../mesh/meshclass.cpp"

using namespace std; 


struct rhs_term;
// ======================================================================
// Data structure to store the traits of fields.
struct FieldTraits {
    std::string name = "fieldTemp";
    int rank=0;
    int priority=1;
    std::string location = "both";
    std::string boun_cond = "periodic";
    std::string init_cond = "none";
    std::string expo_data = "on";
    std::string equation = "";
    Mesh* mesh_ptr=NULL;
    std::vector<rhs_term> rhs_terms;    
};

// =======================================================================
class Field {

    
    public:

    FieldTraits traits_host;
    FieldTraits* traits_dev_ptr;
    
    // Field data
    double* f_host[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_dev[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_rhs[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_lhs[5]={NULL,NULL,NULL,NULL,NULL};
    // Use to store temp field to write to file;
    double* f_temp_host=NULL;
    
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
    double * laplace_dev=NULL;
    double * bi_laplace_dev=NULL;
    
    
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
    void setFieldConstCPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby);
    void setFieldConstGPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby);

    void setRhsTerms(vector<rhs_term> rhs_terms_t);
    // ===================================================================
    // Differential operators
    double* getLaplaceCPU(int i_field,std::string method);
    double* getBiLaplaceCPU(int i_field,std::string method);
    double* getLaplaceGPU(int i_field, std::string method);    

    // -------------------------------------------------------------------
    // Field boundary condition    
    void applyBounCondPeriCPU(int i_field);
    void applyBounCondPeriGPU(int i_field);
    void applyBounCondPeriAnyCPU(double* f_t);

    // -------------------------------------------------------------------
    void updateAnyFieldDev (double* f_dev_ptr, double * f_host_ptr);
    void updateAnyFieldHost (double* f_host_ptr, double * f_dev_ptr);
    void updateMainFieldDev();
    void updateMainFieldHost();

    // -------------------------------------------------------------------
    // Export field
    void export_conf(string str_t, int include_boun);
    void export_conf_any(double* f_t, string f_name, string str_t, int include_boun, string location);
    void getEqn();
    
    // ------------------------------------------------------------------
    // Fetch field traits
    string name () {
        return traits_host.name;
    };
    int rank () {
        return traits_host.rank;
    };
    int priority () {
        return traits_host.priority;
    };
    string location () {
        return traits_host.location;
    };
    string bounCond() {
        return traits_host.boun_cond;
    };
    string initCond() {
        return traits_host.init_cond;
    };
    string expoData() {
        return traits_host.expo_data;
    };
    string equation() {
        return traits_host.equation;
    };
    Mesh* meshPtr() {
        return traits_host.mesh_ptr;
    };    
    vector<rhs_term> rhsTerms() {
        return traits_host.rhs_terms;
    };
    int gridNumberAll() {
        return ((*traits_host.mesh_ptr).host.grid_number.x+
        2*(*traits_host.mesh_ptr).host.grid_number_boun.x)*
            ((*traits_host.mesh_ptr).host.grid_number.y+
            2*(*traits_host.mesh_ptr).host.grid_number_boun.y)*
            ((*traits_host.mesh_ptr).host.grid_number.z+
            2*(*traits_host.mesh_ptr).host.grid_number_boun.z);
    };
    int gridNumberBulk() {
        return (*traits_host.mesh_ptr).host.grid_number.x*
            (*traits_host.mesh_ptr).host.grid_number.y*
            (*traits_host.mesh_ptr).host.grid_number.z;
    };
    int spaceDim() {
        return (*traits_host.mesh_ptr).host.space_dim;
    };
    Vector3<int> gridNumber () {
        return (*traits_host.mesh_ptr).host.grid_number;
    };
    Vector3<int> gridNumberBoun () {
        return (*traits_host.mesh_ptr).host.grid_number_boun;
    };    
    Vector3<double> boxSize () {
        return (*traits_host.mesh_ptr).host.box_size;
    };
    Vector3<double> gridSize () {
        return (*traits_host.mesh_ptr).host.grid_size;
    };    
    
    // ===================================================================
};

__global__ void getLaplaceGPUCore(double* f, int Nx, int Ny, int Nbx, int Nby, double dx, double dy);

#endif
