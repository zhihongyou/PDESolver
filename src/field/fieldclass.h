#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../mesh/meshclass.cpp"
#include "../utility/fieldFunctionClass.h"
#include "../utility/finiteDifferenceCentralO2Isotropic2D.h"
#include "../utility/finiteDifferenceCentralO4Isotropic2D.h"


using namespace std; 

class Field;
// -----------------------------------------------------------------------
struct fieldFunction {
    // ******************************************************************
    // Operation struct for a single function call to a field
    // f_perator is a string used to determine function call
    // field_ptr is the field on which the f_operator act.
    // ******************************************************************
    
    string f_operator;
    Field* field_ptr;
    
    fieldFunction (string f_operator1, Field* field_ptr1) {
        f_operator=f_operator1;
        field_ptr=field_ptr1;
    };
};

// -----------------------------------------------------------------------
struct rhsTerm {
    // ******************************************************************
    // Operation struct for collecting a single RHS term of a field.
    // f_funcs provides information on all functions of fields, which are
    //   to be multiplied to form the RHS term.
    // prefactor is the prefactor of this term.
    // scheme determines whether this term is calculated explicitly or
    //   implicitly
    // ******************************************************************
    
    vector<fieldFunction> f_funcs;
    double prefactor;
    string scheme;
    
    rhsTerm(double prefactor1, vector<fieldFunction> f_funcs1, string scheme1) {
        prefactor=prefactor1;
        f_funcs=f_funcs1;
        scheme=scheme1;
    };
};

// ----------------------------------------------------------------------
struct rhsPtrs {
    // ******************************************************************
    // Struct that determines the RHS of a field
    // num_terms is the number of terms on the RHS
    // prefactors are collections of prefactors of each term.
    // num_func_1term is the number of fields in each term
    // f_func_ptrs collect all pointers on the RHS.
    // ******************************************************************
    int* num_terms;
    double* prefactors;
    int* num_funcs_1term;
    double** f_func_ptrs;
    int* schemes;
};


// ======================================================================
// Data structure to store the traits of fields.
struct FieldTraits {
    string name = "fieldTemp";
    int rank=0;
    int priority=1;
    string location = "both";
    string boun_cond = "periodic";
    string init_cond = "none";
    string expo_data = "on";
    string equation = "";
    Mesh* mesh_ptr=NULL;        
};

// =======================================================================
class Field {
    
    
    public:

    FieldTraits traits_host;
    FieldTraits traits_dev_ptr;
    
    vector<rhsTerm> rhs_terms;
    // These will store main pointers to field functions.
    vector<fieldFunction> f_funcs_rhs;
    // These will store device pointers to field functions.
    rhsPtrs rhs_ptrs_host;
    rhsPtrs rhs_ptrs_dev;

    // Array of finite difference schemes
    FiniteDifference ** FDM_ptrs;
    // Index of finite difference scheme
    int FDM_idx;
    
    // Use to store temp field to write to file;
    double* f_host[5]={NULL,NULL,NULL,NULL,NULL};
    double* f_temp_host=NULL;
    
    // Field data: those are to be stored on host or device depending
    //   on the "device" used.
    // In case "GPU" is used, these are host pointers to device data.
    // Host can use f[0][idx], but device cannot.
    // Instead, device can only use f_dp[0][idx], where f_dp is a device
    //   pointer point to device data.
    double* f[5]={NULL,NULL,NULL,NULL,NULL};
    double* rhs[5]={NULL,NULL,NULL,NULL,NULL};
    double* lhs[5]={NULL,NULL,NULL,NULL,NULL};    
    
    // Operators
    double* f_now=NULL;
    double* d1x=NULL;
    double* d1y=NULL;
    double* d1z=NULL;
    double* d1x1y=NULL;
    double* d1x1z=NULL;
    double* d1y1z=NULL;
    double* d2x=NULL;
    double* d2y=NULL;
    double* d2z=NULL;
    double* d2x2y=NULL;
    double* d2x2z=NULL;
    double* d2y2z=NULL;    
    double* laplace=NULL;
    double* bi_laplace=NULL;
    double* one_over_f=NULL;
    double* sinf=NULL;
    double* cosf=NULL;
    double* tanf=NULL;
    double* cotf=NULL;
    
    
    // ===================================================================
    // Methods
    
    // ===================================================================
    // Constructor
    Field (Mesh* mesh_ptr_t, string name_t, int rank_t, int priority_t, string boun_cond_t, string init_cond_t, string expo_data_t);    
    
    // ------------------------------------------------------------------
    // Get rhs
    void getRHS(int i_f_copy);
    void allocField (double* &f_t, string location);
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
    void setRhsTerms(vector<rhsTerm> rhs_terms_t);
    // ===================================================================
    // Functions of fields
    double* getLaplaceCPU(int i_field,string method);
    double* getLaplaceGPU(int i_field, string method);
    double* getBiLaplaceCPU(int i_field,string method);    
    double* getBiLaplaceGPU(int i_field, string method);
    double* getFieldFunctionCPU(int i_field, string method);    
    template<class FDM_class>
    double* getFFuncCPU(double* f_func_ptr, int i_field, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), string method);
    template<class FDM_class>
    double* getFFuncGPU(double* f_func_ptr, int i_field, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), string method);

    // -------------------------------------------------------------------
    // Field boundary condition    
    void applyBounCondPeriCPU(double* f_t);
    void applyBounCondPeriGPU(double* f_t);
    void applyBounCondPeriAnyCPU(double* f_t);

    // -------------------------------------------------------------------
    void updateAnyFieldDev (double* f_dev_ptr, double * f_host_ptr);
    void updateAnyFieldHost (double* f_host_ptr, double * f_dev_ptr);
    void updateMainFieldDev();
    void updateMainFieldHost();

    // -------------------------------------------------------------------
    // Export field
    void export_conf(string str_t, string device, int include_boun);
    void export_conf_any(double* f_t, string f_name, string str_t, string location_t, int include_boun);
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
    vector<rhsTerm> rhsTerms() {
        return rhs_terms;
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
    rhsPtrs fieldFuncPtrs () {
        return rhs_ptrs_host;
    };
    
    // ===================================================================
};

// __global__ void getLaplaceGPUCore(double* f, int Nx, int Ny, int Nbx, int Nby, double dx, double dy);

#endif
