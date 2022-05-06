#ifndef FIELDCLASS_FIELDFUNCTIONS_GENERIC_CU
#define FIELDCLASS_FIELDFUNCTIONS_GENERIC_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>


using namespace std;


// -----------------------------------------------------------------------
// To call this, use:
// FDMCentralO4Iso2D FDM_test;
// phia.getFFuncCPU<FDMCentralO4Iso2D>(phia.laplace, 0, FDM_test, &FDMCentralO4Iso2D::laplace, "new");
template<class FDM_class>
double* Field::getFFuncCPU1(double* f_func_ptr, int i_field, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), string method="new") {
    int get_new=1;
    if (method=="old" && f_func_ptr != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(f_func_ptr, "cpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        int di=1;
        int dj=Nx+2*Nbx;
        double dx=gridSize().x;
        double dy=gridSize().y;

        for (int j=0; j<Ny;j++) {
            for (int i=0; i<Nx; i++) {            
                int idx=(j+Nby)*dj+i+Nbx;
                f_func_ptr[idx]=(FDM_scheme.*f_func)(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return f_func_ptr;
};


// -----------------------------------------------------------------------------------------
template <typename T>
T* Field::getFFuncCPU(T* &f_func_ptr, int i_field, T f_func(double*,int), string method="new") {
    // To call this function, use:
    // phia.getFFuncCPU(phia.one_over_f, 0, FieldFunction::oneOverF, "new");
    int get_new=1;
    if (method=="old" && f_func_ptr != NULL) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<T>(f_func_ptr, "cpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        int dj=Nx+2*Nbx;

        for (int j=0; j<Ny;j++) {
            for (int i=0; i<Nx; i++) {            
                int idx=(j+Nby)*dj+i+Nbx;
                f_func_ptr[idx]=f_func(f[i_field],idx);
            };
        };
    };
    return f_func_ptr;
};


// -----------------------------------------------------------------------------------------
template <typename T>
T* Field::getFFuncGPU(T* &f_func_ptr, int i_field, T f_func(double*,int), string method="new") {
    // To call this function, use:
    // phia.getFFuncGPU(phia.one_over_f, 0, FieldFunction::oneOverF, "new");
    int get_new=1;
    if (method=="old" && f_func_ptr != NULL) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(f_func_ptr, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;        
        getFFuncGPUCore<T><<<Ny,Nx>>>(f_func_ptr,f[i_field],f_func,Nx,Ny,Nbx,Nby);
    };
    return f_func_ptr;
};


// -----------------------------------------------------------------------
// This function is NOT working!!!!
template<class FDM_class>
double* Field::getFFuncGPU1(double* f_func_ptr, int i_field, FDM_class& FDM_scheme, double (FDM_class::*f_func)(double*,int,int,int,double,double), string method="new") {
    int get_new=1;
    if (method=="old" && f_func_ptr != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(f_func_ptr, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        
        // getFFuncGPUCore<FDM_class><<<Ny,Nx>>>(f_func_ptr,f[i_field],FDM_scheme,&FDM_class::*f_func,Nx,Ny,Nbx,Nby,dx,dy);
        // getFFuncGPUCore<*FDM_class><<<Ny,Nx>>>(f_func_ptr,f[i_field],FDM_scheme,&f_func,Nx,Ny,Nbx,Nby,dx,dy);
    };

    return f_func_ptr;
};


// =================================================================

#endif
