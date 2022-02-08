#ifndef FIELDCLASS_CPP
#define FIELDCLASS_CPP

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>
#include "fieldclass.h"

using namespace std;

// =======================================================================
// Constructor
Field::Field (Mesh* mesh_ptr_t, std::string name_t, int rank_t, int priority_t, std::string boun_cond_t, std::string init_cond_t, std::string expo_data_t) {
    mesh_ptr=mesh_ptr_t;
    name=name_t;
    rank=rank_t;
    priority=priority_t;
    boun_cond=boun_cond_t;
    init_cond=init_cond_t;
    expo_data=expo_data_t;
    if (init_cond=="Gaussian") {
        initFieldGaus(1,1,1);
    } else if (init_cond=="ones") {
        initFieldConst(1);
    } else if (init_cond=="sin") {
        initFieldSin(1,1,0);        
    };
    applyBounCondPeri(0);
    if (priority==0) {
        f_rhs[0] = new double[(*mesh_ptr).getGridNumberAll()];
    };
};

// -----------------------------------------------------------------------
void Field::initFieldConst(double f_value) {
    if (f_cpu[0] == NULL) {
        f_cpu[0] = new double[(*mesh_ptr).getGridNumberAll()];
    };
    int Nx=(*mesh_ptr).grid_number.x;
    int Ny=(*mesh_ptr).grid_number.y;
    int Nbx=(*mesh_ptr).grid_number_boun.x;
    int Nby=(*mesh_ptr).grid_number_boun.y;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            f_cpu[0][idx]=f_value;
        };
    };
};

// -----------------------------------------------------------------------
void Field::initFieldGaus(double r_center, double r_decay, double gaus_amplitude) {
    if (f_cpu[0] == NULL) {
        f_cpu[0] = new double[(*mesh_ptr).getGridNumberAll()];
    };
    int Nx=(*mesh_ptr).grid_number.x;
    int Ny=(*mesh_ptr).grid_number.y;
    int Nbx=(*mesh_ptr).grid_number_boun.x;
    int Nby=(*mesh_ptr).grid_number_boun.y;
    double r2m=0.25*Nx*Nx+0.25*Ny*Ny;    
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=(i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny);
            f_cpu[0][idx]=gaus_amplitude*(exp(-20*r2/r2m)-0.5);
        };
    };
};

// -----------------------------------------------------------------------
void Field::initFieldSin(double sin_amplitude=1, int sin_period=1, double sin_phase=0) {
    if (f_cpu[0] == NULL) {
        f_cpu[0] = new double[(*mesh_ptr).getGridNumberAll()];
    };
    int Nx=(*mesh_ptr).grid_number.x;
    int Ny=(*mesh_ptr).grid_number.y;
    int Nbx=(*mesh_ptr).grid_number_boun.x;
    int Nby=(*mesh_ptr).grid_number_boun.y;
    double r2m=0.25*Nx*Nx+0.25*Ny*Ny;    
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=(i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny);
            f_cpu[0][idx]=sin_amplitude*sin(2*M_PI*sin_period*i/Nx)*sin(2*M_PI*sin_period*j/Ny);
        };
    };
};

// =======================================================================
// Differential operators
// -----------------------------------------------------------------------
double* Field::getLaplace(int i_field,std::string method="new") {
    int get_new=1;
    if (method=="old" && laplace != nullptr) {
        get_new=0;
    };

    if (get_new==1) {
        if (laplace == NULL) {
            laplace=new double[(*mesh_ptr).getGridNumberAll()];
        };
        int Nx=(*mesh_ptr).grid_number.x;
        int Ny=(*mesh_ptr).grid_number.y;
        int Nbx=(*mesh_ptr).grid_number_boun.x;
        int Nby=(*mesh_ptr).grid_number_boun.y;
        int di=1;
        int dj=Nx+2*Nbx;
        double dx=(*mesh_ptr).grid_size.x;
        double dy=(*mesh_ptr).grid_size.y;

        for (int j=0; j<Ny;j++) {
            for (int i=0; i<Nx; i++) {            
                int idx=(j+Nby)*dj+i+Nbx;
                laplace[idx]=1.0/(dx*dy)*( -21.0/5.0*f_cpu[i_field][idx]
			  +13.0/15.0*( f_cpu[i_field][idx+di] + f_cpu[i_field][idx-di] + f_cpu[i_field][idx+dj] + f_cpu[i_field][idx-dj] )
			  +4.0/15.0*( f_cpu[i_field][idx+di+dj] + f_cpu[i_field][idx-di+dj] + f_cpu[i_field][idx+di-dj] + f_cpu[i_field][idx-di-dj] )
			  -1.0/60.0*( f_cpu[i_field][idx+2*di] + f_cpu[i_field][idx-2*di] + f_cpu[i_field][idx+2*dj] + f_cpu[i_field][idx-2*dj] )
			  -1.0/30.0*( f_cpu[i_field][idx+di+2*dj] + f_cpu[i_field][idx-di+2*dj] + f_cpu[i_field][idx+di-2*dj] + f_cpu[i_field][idx-di-2*dj]
				     + f_cpu[i_field][idx+2*di+dj] + f_cpu[i_field][idx-2*di+dj] + f_cpu[i_field][idx+2*di-dj] + f_cpu[i_field][idx-2*di-dj] )
			  );
            };
        };
    };
    return laplace;
};

// -----------------------------------------------------------------------
double* Field::getBiLaplace(int i_field,std::string method="new") {
    int get_new=1;
    if (method=="old" && laplace != nullptr) {
        get_new=0;
    };

    if (get_new==1) {
        if (bi_laplace == NULL) {
            bi_laplace=new double[(*mesh_ptr).getGridNumberAll()];
        };
        int Nx=(*mesh_ptr).grid_number.x;
        int Ny=(*mesh_ptr).grid_number.y;
        int Nbx=(*mesh_ptr).grid_number_boun.x;
        int Nby=(*mesh_ptr).grid_number_boun.y;
        int di=1;
        int dj=Nx+2*Nbx;
        double dx=(*mesh_ptr).grid_size.x;
        double dy=(*mesh_ptr).grid_size.y;

        for (int j=0; j<Ny;j++) {
            for (int i=0; i<Nx; i++) {            
                int idx=(j+Nby)*dj+i+Nbx;
                bi_laplace[idx]=1.0/(dx*dx*dy*dy)*( 779.0/45.0*f_cpu[i_field][idx]
			  -191.0/45.0*( f_cpu[i_field][idx+di] + f_cpu[i_field][idx-di] + f_cpu[i_field][idx+dj] + f_cpu[i_field][idx-dj] )
			  -187.0/90.0*( f_cpu[i_field][idx+di+dj] + f_cpu[i_field][idx-di+dj] + f_cpu[i_field][idx+di-dj] + f_cpu[i_field][idx-di-dj] )
			  +7.0/30.0*( f_cpu[i_field][idx+2*di] + f_cpu[i_field][idx-2*di] + f_cpu[i_field][idx+2*dj] + f_cpu[i_field][idx-2*dj] )
			  +47.0/45.0*( f_cpu[i_field][idx+di+2*dj] + f_cpu[i_field][idx-di+2*dj] + f_cpu[i_field][idx+di-2*dj] + f_cpu[i_field][idx-di-2*dj]
				       + f_cpu[i_field][idx+2*di+dj] + f_cpu[i_field][idx-2*di+dj] + f_cpu[i_field][idx+2*di-dj] + f_cpu[i_field][idx-2*di-dj] )
			  -29.0/180.0*( f_cpu[i_field][idx+2*di+2*dj] + f_cpu[i_field][idx-2*di+2*dj] + f_cpu[i_field][idx+2*di-2*dj] + f_cpu[i_field][idx-2*di-2*dj] )
			  +1.0/45.0*( f_cpu[i_field][idx+3*di] + f_cpu[i_field][idx-3*di] + f_cpu[i_field][idx+3*dj] + f_cpu[i_field][idx-3*dj] )
			  -17.0/180.0*( f_cpu[i_field][idx+di+3*dj] + f_cpu[i_field][idx-di+3*dj] + f_cpu[i_field][idx+di-3*dj] + f_cpu[i_field][idx-di-3*dj]
				       + f_cpu[i_field][idx+3*di+dj] + f_cpu[i_field][idx-3*di+dj] + f_cpu[i_field][idx+3*di-dj] + f_cpu[i_field][idx-3*di-dj] )
			  );
            };
        };
    };
    return bi_laplace;
};


//=======================================================================
void Field::applyBounCondPeri(int i_field) {
    applyBounCondPeriAny(f_cpu[i_field]);
};


//=======================================================================
void Field::applyBounCondPeriAny(double* f_t) {
    int Nx=(*mesh_ptr).grid_number.x;
    int Ny=(*mesh_ptr).grid_number.y;
    int Nbx=(*mesh_ptr).grid_number_boun.x;
    int Nby=(*mesh_ptr).grid_number_boun.y;
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
void Field::export_conf(string str_t, int include_boun=0) {
    export_conf_any(f_cpu[0],name,str_t,include_boun);
}

// ----------------------------------------------------------------------
void Field::export_conf_any(double* f_t, string f_name, string str_t, int include_boun=0) {
    ofstream conf_file;
    int PrecData=8;
    std::string conf_file_name="data/"+f_name+"_"+ str_t + ".dat";
    conf_file.open(conf_file_name.c_str() );

    int idx;
    int Nx=(*mesh_ptr).grid_number.x;
    int Ny=(*mesh_ptr).grid_number.y;
    int Nbx=(*mesh_ptr).grid_number_boun.x;
    int Nby=(*mesh_ptr).grid_number_boun.y;
    int* idx0=new int [4];
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
            conf_file<<setiosflags(ios::scientific) <<setprecision(PrecData) <<f_t[idx]<<endl;
        }
    }
    conf_file.close();
}

// ----------------------------------------------------------------------
void Field::getEqn() {
    // Loop over fields.
    // for (auto rhs_term_i : rhs_terms) {            
    //     double* f_func_ptrs[10];
    //         // cout <<"term "<<N_funcs<<" with operator "<<rhs_term_i.rhs_operator<<" and number "<<rhs_term_i.numbers[0]<<" :";
    //         // Evaluate each operator applied on field
    //     for (auto f_func_i : rhs_term_i.f_function) {
    //             // cout<<" "<<f_func_i.f_operator<<"("<<(*f_func_i.field_ptr).name <<"), ";
    //             evalOperator(f_ptr_i,f_func_i,f_func_ptrs[N_funcs],i_field);
    //         };
    //         // cout <<endl;
    //         // Add one RHS term to field increments
    //         addRHSTerm(f_ptr_i,i_field,rhs_term_i,f_func_ptrs,N_funcs);
    //     };
    // };
};


// -----------------------------------------------------------------------
// Operation struct for a single function call to a field
struct field_function {
    std::string f_operator;
    Field* field_ptr;
    // Field* field_func_ptr;
    field_function (std::string f_operator1,Field* field_ptr1) {
        f_operator=f_operator1;
        field_ptr=field_ptr1;
        // field_func_ptr=field_func_ptr1;
    };
};


// -----------------------------------------------------------------------
// Operation struct for collecting RHS of field.
struct rhs_term {
    std::string rhs_operator;
    std::vector<field_function> f_function;
    std::vector<double> numbers;    
    std::string scheme;
    
    rhs_term(std::string rhs_operator1, std::vector<field_function> f_function1, std::vector<double> numbers1, std::string scheme1) {
        rhs_operator=rhs_operator1;
        f_function=f_function1;
        numbers=numbers1;
        scheme=scheme1;
    };
};


#endif
