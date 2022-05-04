#ifndef FIELDCLASS_CU
#define FIELDCLASS_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "fieldclass.h"
#include "fieldclassGPU.cu"


using namespace std;


// =======================================================================
// Constructor
Field::Field (Mesh* mesh_ptr_t, string name_t, int rank_t, int priority_t, string boun_cond_t, string init_cond_t, string expo_data_t) {
    traits_host.mesh_ptr=mesh_ptr_t;
    traits_host.name=name_t;
    traits_host.rank=rank_t;
    traits_host.priority=priority_t;
    traits_host.boun_cond=boun_cond_t;
    traits_host.init_cond=init_cond_t;
    traits_host.expo_data=expo_data_t;
    // Initiate field on host, which will then be copied to f.
    if (traits_host.init_cond=="Gaussian") {
        initFieldGaus(0,10,1);
    } else if (traits_host.init_cond=="ones") {
        initFieldConst(1);
    } else if (traits_host.init_cond=="sin") {
        initFieldSin(0.01,4,0);        
    };
};

// -----------------------------------------------------------------------
void Field::initFieldConst(double f_value) {
    if (f_host[0] == NULL) {
        f_host[0] = new double[gridNumberAll()];
    };
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            f_host[0][idx]=f_value;
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};

// -----------------------------------------------------------------------
void Field::initFieldGaus(double r_center, double r_decay, double gaus_amplitude) {
    if (f_host[0] == NULL) {
        f_host[0] = new double[gridNumberAll()];
    };
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx2=gridSize().x*gridSize().x;
    double dy2=gridSize().y*gridSize().y;
    double rd2=r_decay*r_decay;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=dx2*(i-0.5*Nx)*(i-0.5*Nx)+dy2*(j-0.5*Ny)*(j-0.5*Ny);
            f_host[0][idx]=gaus_amplitude*(exp(-r2/rd2));
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};

// -----------------------------------------------------------------------
void Field::initFieldSin(double sin_amplitude=1, int sin_period=1, double sin_phase=0) {
    if (f_host[0] == NULL) {
        f_host[0] = new double[gridNumberAll()];
    };
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double r2m=0.25*Nx*Nx+0.25*Ny*Ny;    
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=(i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny);
            f_host[0][idx]=sin_amplitude*sin(2*M_PI*sin_period*i/Nx)*sin(2*M_PI*sin_period*j/Ny);
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};


// ----------------------------------------------------------------------
void Field::setFieldConstCPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    setFieldConstCPUCore(f_t, f_val, Nx, Ny, Nbx, Nby);
};

// ----------------------------------------------------------------------
void Field::setFieldConstGPU(double* f_t, double f_val, int Nx, int Ny, int Nbx, int Nby) {
    setFieldConstGPUCore<<<Ny,Nx>>>(f_t, f_val, Nx, Ny, Nbx, Nby);
};


// ----------------------------------------------------------------------
void Field::setRhsTerms(vector<rhsTerm> rhs_terms_t) {
    rhs_terms=rhs_terms_t;
};

// =======================================================================
// Differential operators
// -----------------------------------------------------------------------
double* Field::getLaplaceCPU(int i_field, string method="new") {
    int get_new=1;
    if (method=="old" && laplace != NULL) {
        get_new=0;
    };

    if (get_new==1) {
        if (laplace == NULL) {
            laplace=new double[gridNumberAll()];
        };
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
        if (laplace == NULL) {
            cudaMalloc(&laplace, gridNumberAll()*sizeof(double));
        };
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
        if (bi_laplace == NULL) {
            bi_laplace=new double[gridNumberAll()];
        };
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
        if (bi_laplace == NULL) {
            cudaMalloc(&bi_laplace, gridNumberAll()*sizeof(double));
        };
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


//=======================================================================
void Field::applyBounCondPeriCPU(double* f_t) {
    applyBounCondPeriAnyCPU(f_t);
};

//=======================================================================
void Field::applyBounCondPeriGPU(double* f_t) {
    applyBounCondPeriAnyGPU<<<gridNumber().y,gridNumber().x>>>
        (f_t,
        gridNumber().x,gridNumber().y,
        gridNumberBoun().x,gridNumberBoun().y);
};

//=======================================================================
void Field::applyBounCondPeriAnyCPU(double* f_t) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    int dj=Nx+2*Nbx;
    int idx,idx1;
    for (int j=Nby; j<Ny+Nby; j++) {
        for (int i=0; i<Nbx; i++) {
            idx=j*dj+i;
            idx1=idx+Nx;
            f_t[idx]=f_t[idx1];
            f_t[idx1+Nbx]=f_t[idx+Nbx];
        };
    };
    for (int j=0; j<Nby; j++) {
        for (int i=0; i<Nx+2*Nbx; i++) {            
            idx=j*dj+i;
            idx1=idx+dj*Ny;
            f_t[idx]=f_t[idx1];
            f_t[idx1+dj*Nby]=f_t[idx+dj*Nby];
        };
    };

};


// ----------------------------------------------------------------------
void Field::export_conf(string str_t, string device, int include_boun=0) {
    if (device=="cpu") {        
        export_conf_any(f[0],name(),str_t, "cpu", include_boun);
    } else if (device=="gpu") {
        export_conf_any(f[0],name(),str_t, "gpu", include_boun);
    };
}

// ----------------------------------------------------------------------
void Field::export_conf_any(double* f_t, string f_name, string str_t, string location_t="cpu" , int include_boun=0) {
    ofstream conf_file;
    int PrecData=8;
    string conf_file_name="data/"+f_name+"_"+ str_t + ".dat";
    conf_file.open(conf_file_name.c_str() );
    
    int idx;
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    int* idx0=new int [4];

    
    if (location_t=="gpu") {
        if (f_temp_host == NULL) {
            f_temp_host=new double[gridNumberAll()];
        };        
        cudaMemcpy(f_temp_host, f_t, gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
        // cout <<f_temp_host[250]<<endl;
        // updateAnyFieldHost(f_temp_host,f_t);
    };        
    
    if (include_boun==0) {
        idx0[0]=0;
        idx0[1]=0;
        idx0[2]=0;
        idx0[3]=0;
    } else {
        idx0[0]=-Nbx;
        idx0[1]=Nbx;
        idx0[2]=-Nby;
        idx0[3]=Nby;
    };

    for (int j=idx0[2]; j<Ny+idx0[3]; j++) {
        for (int i=idx0[0]; i<Nx+idx0[1]; i++) {        
            idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
            // conf_file<<fixed <<setprecision(PrecData) <<f_t[idx]<<endl;
            if (location_t=="cpu") {
                conf_file<<setiosflags(ios::scientific) <<setprecision(PrecData) <<f_t[idx]<<endl;
            } else {
                conf_file<<setiosflags(ios::scientific) <<setprecision(PrecData) <<f_temp_host[idx]<<endl;
            };
        }
    }
    conf_file.close();
}

// ----------------------------------------------------------------------
// string Field::equation() {
//     string eqn;
//     if (priority==0) {
//         eqn="p_t "+name()+"=";
//     } else {
//         eqn=name()+"=";
//     };
//     for (auto rhs_term_i : rhsTerms()) {
//         if (rhs_term_i != rhsTerms().begin()) {
            
//         };
//         for (auto f_func_i : rhs_term_i.f_function) {
//             eqn=eqn+
//                 cout<<"*"<<f_func_i.f_operator<<"("<<(*f_func_i.field_ptr).name() <<")";
//                 evalOperator(f_ptr_i,f_func_i,f_func_ptrs[N_funcs],i_field);
//                 N_funcs+=1;
//             };
//             addRHSTerm(f_ptr_i,i_field,rhs_term_i,f_func_ptrs,N_funcs);
//         };

//         if ((*f_ptr_i).priority()>0 && (*f_ptr_i).bounCond()=="periodic") {
//             if (device=="cpu") {
//                 (*f_ptr_i).applyBounCondPeriCPU(i_field);
//             } else if (device=="gpu"){
//                 (*f_ptr_i).applyBounCondPeriGPU(i_field);
//             };
//         };
// };

// -------------------------------------------------------------------
// Copy any field data from CPU to GPU
void Field::updateAnyFieldDev (double* f_dev_ptr, double * f_host_ptr) {
    if (f_dev_ptr == NULL) {
        cudaMalloc(&f_dev_ptr, gridNumberAll()*sizeof(double));
    };
    cudaMemcpy(f_dev_ptr, f_host_ptr, gridNumberAll()*sizeof(double),cudaMemcpyHostToDevice);
};

// -------------------------------------------------------------------
// Copy any field data from GPU to CPU
void Field::updateAnyFieldHost (double* f_host_ptr, double * f_dev_ptr) {    
    cudaMemcpy(f_host_ptr, f_dev_ptr, gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
};

// -------------------------------------------------------------------
// Copy main field data from CPU to GPU
void Field::updateMainFieldDev () {
    if (f[0] == NULL) {
        cudaMalloc(&f[0], gridNumberAll()*sizeof(double));
    };
    cudaMemcpy(f[0], f_host[0], gridNumberAll()*sizeof(double),cudaMemcpyHostToDevice);
};

// -------------------------------------------------------------------
// Copy main field data from GPU to CPU
void Field::updateMainFieldHost () {
    if (f_host[0] == NULL) {
        f_host[0]=new double[gridNumberAll()];
    };
    cudaMemcpy(f_host[0], f[0], gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
};

// ------------------------------------------------------------------
void Field::allocField (double* f_t, string location) {
    if (f_t==NULL) {
        if (location=="cpu") {
            f_t=new double[gridNumberAll()];
        } else {
            cudaMalloc(&f_t, gridNumberAll()*sizeof(double));
        };
    };    
};

// =================================================================

#endif
