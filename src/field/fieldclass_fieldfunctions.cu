#ifndef FIELDCLASS_FIELDFUNCTIONS_CU
#define FIELDCLASS_FIELDFUNCTIONS_CU

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
double* Field::getFNowCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && f_now != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(f_now, "cpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        int dj=Nx+2*Nbx;

        for (int j=0; j<Ny;j++) {
            for (int i=0; i<Nx; i++) {            
                int idx=(j+Nby)*dj+i+Nbx;
                f_now[idx]=f[i_field][idx];
            };
        };
    };
    return f_now;
};


// -----------------------------------------------------------------------
double* Field::getFNowGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && f_now != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(f_now, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        fieldCopy<<<Ny,Nx>>>(f_now, f[i_field], Nx, Ny, Nbx, Nby);
    };
    return f_now;
};


// -----------------------------------------------------------------------
double* Field::getD1xCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1x != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(d1x, "cpu");
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
                d1x[idx]=FDM_ptrs[FDM_idx]->d1x(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return d1x;
};


// -----------------------------------------------------------------------
double* Field::getD1xGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1x != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(d1x, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getD1xGPUCore<<<Ny,Nx>>>(d1x,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return d1x;
};


// -----------------------------------------------------------------------
double* Field::getD1yCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1y != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(d1y, "cpu");
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
                d1y[idx]=FDM_ptrs[FDM_idx]->d1y(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return d1y;
};


// -----------------------------------------------------------------------
double* Field::getD1yGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1y != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(d1y, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getD1yGPUCore<<<Ny,Nx>>>(d1y,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return d1y;
};


// -----------------------------------------------------------------------
double* Field::getD2xCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d2x != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(d2x, "cpu");
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
                d2x[idx]=FDM_ptrs[FDM_idx]->d2x(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return d2x;
};


// -----------------------------------------------------------------------
double* Field::getD2xGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d2x != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(d2x, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getD2xGPUCore<<<Ny,Nx>>>(d2x,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return d2x;
};


// -----------------------------------------------------------------------
double* Field::getD2yCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d2y != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(d2y, "cpu");
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
                d2y[idx]=FDM_ptrs[FDM_idx]->d2y(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return d2y;
};


// -----------------------------------------------------------------------
double* Field::getD2yGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d2y != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(d2y, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getD2yGPUCore<<<Ny,Nx>>>(d2y,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return d2y;
};


// -----------------------------------------------------------------------
double* Field::getD1x1yCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1x1y != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(d1x1y, "cpu");
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
                d1x1y[idx]=FDM_ptrs[FDM_idx]->d1x1y(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return d1x1y;
};


// -----------------------------------------------------------------------
double* Field::getD1x1yGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && d1x1y != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(d1x1y, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getD1x1yGPUCore<<<Ny,Nx>>>(d1x1y,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return d1x1y;
};


// -----------------------------------------------------------------------
double* Field::getLaplaceCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && laplace != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(laplace, "cpu");
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
                laplace[idx]=FDM_ptrs[FDM_idx]->laplace(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return laplace;
};


// -----------------------------------------------------------------------
double* Field::getLaplaceGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && laplace != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(laplace, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getLaplaceGPUCore<<<Ny,Nx>>>(laplace,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return laplace;
};


// -----------------------------------------------------------------------
double* Field::getBiLaplaceCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && bi_laplace != nullptr) {
        get_new=0;
    };

    if (get_new==1) {
        allocField<double>(bi_laplace, "cpu");
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
                bi_laplace[idx] = FDM_ptrs[FDM_idx]->bi_laplace(f[i_field],idx,di,dj,dx,dy);
            };
        };
    };
    return bi_laplace;
};

// -----------------------------------------------------------------------
double* Field::getBiLaplaceGPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && bi_laplace != nullptr) {
        get_new=0;
    };
    
    if (get_new==1) {
        allocField<double>(bi_laplace, "gpu");
        int Nx=gridNumber().x;
        int Ny=gridNumber().y;
        int Nbx=gridNumberBoun().x;
        int Nby=gridNumberBoun().y;
        double dx=gridSize().x;
        double dy=gridSize().y;
        getBiLaplaceGPUCore<<<Ny,Nx>>>(bi_laplace,f[i_field],FDM_ptrs,FDM_idx,Nx,Ny,Nbx,Nby,dx,dy);
    };
    return bi_laplace;
};

// =================================================================

#endif
