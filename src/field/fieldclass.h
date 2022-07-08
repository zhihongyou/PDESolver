#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../mesh/meshclass.cpp"
#include "fieldFunction.cu"


using namespace std; 


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
    
    // Array of finite difference schemes
    string FDMScheme="CentralDifferenceO2Iso2D";
    
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
    double** f_funcs_host=new double*[200];
    double** f_funcs_dev;
    
    // Functions of fields
    vector<rhsTerm> rhs_terms;
    // These will store main pointers to field functions.
    // f_func_rhs is replaced by f_func_rhs_host
    // vector<FFuncDef> f_funcs_rhs;
    // These will store device pointers to field functions.
    rhsPtrs rhs_ptrs_host;
    rhsPtrs rhs_ptrs_dev;    
    // Matrices that store field data
    
    FFuncType* f_funcs=new FFuncType[200];
    // int num_f_funcs_rhs;
    // Number of functions to be used
    int num_f_funcs=0;
    FFuncItem* f_funcs_rhs=new FFuncItem[200];
    // FFuncItem* f_funcs_rhs_dev=new FFuncItem[200];
    // List that link function name to its index in f_funcs.
    map<pair<string,string>, FFuncItem> f_func_map;

    // Specify specialty of fields that can trigger special events.
    string specialty="";
    
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
    Field () {};
    Field (Mesh* mesh_ptr_t, string name_t);
    Field (Mesh* mesh_ptr_t, string name_t, int priority_t);
    Field (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    Field (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // ------------------------------------------------------------------
    // Get rhs
    void getRHS(int i_f_copy);
    template <typename T>
    void allocField(T* &f_t, string location);
    double* getFFuncPtr(string f_operator);
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
    void setFFuncArgs();
    // ===================================================================    

    double* getFFuncByName(string f_operator, int i_field, FFuncArgs f_func_args);
    double* getFFuncCPU(double* &f_func_ptr, int i_field, FFuncType f_func, FFuncArgs f_func_args, string method);
    double* getFFuncGPU(double* &f_func_ptr, int i_field, FFuncType f_func, FFuncArgs f_func_args, string method);
    void addFunctoRHS(FFuncDef f_func_i, string device, string func_scheme);
    
    template <typename T>
    T* getFFuncCPU1(T* &f_func_ptr, int i_field, T f_func(double*,int), string method);
    
    template <typename T>
    T* getFFuncGPU1(T* &f_func_ptr, int i_field, T f_func(double*,int), string method);    

    template<class FFuncClass>
    double* getFFuncCPU2(double* &f_func_ptr, int i_field, FFuncClass& FFuncScheme, double (FFuncClass::*f_func)(double*,int,int,int,double,double), string method);       
    
    template<class FFuncClass>
    double* getFFuncGPU2(double* &f_func_ptr, int i_field, FFuncClass& FFunccheme, double (FFuncClass::*f_func)(double*,int,int,int,double,double), string method);

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
    void export_f_func(string f_operator, string str_t, string device, int include_boun);
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



#endif
