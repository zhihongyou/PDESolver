#ifndef FIELDCLASS_INITIALCONDITION_CU
#define FIELDCLASS_INITIALCONDITION_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>

using namespace std;
using std::default_random_engine;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
int seed = time(0);
default_random_engine rng(seed);


// --------------------------------------------------------------
void Field::initFieldImport(string str_t="0", int include_boun=0) {            
    int idx;
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
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
    string conf_file_name=dire_expo()+name()+"_"+ str_t + ".dat";
    FILE * conf_file;
    conf_file = fopen(const_cast<char*>(conf_file_name.c_str()), "r");
    if (conf_file == NULL ) {
        cout << "Unable to open file " <<conf_file_name <<endl;
        exit(1);
    };

    for (int j=idx0[2]; j<Ny+idx0[3]; j++) {
        for (int i=idx0[0]; i<Nx+idx0[1]; i++) {
            idx=(Nx+2*Nbx)*(j+Nby)+i+Nbx;
            fscanf(conf_file, "%lf", &f_host[0][idx]);            
        }
    }
    fclose(conf_file);
    applyBounCondPeriCPU(f_host[0]);
};


// --------------------------------------------------------------
void Field::initFieldConst(double f_value, double f_dev=0.0000000000000001) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    uniform_real_distribution<double> randUR;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            f_host[0][idx]=f_value+f_dev*(randUR(rng)-0.5);
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};

// --------------------------------------------------------------
void Field::initFieldGaus(double r_center, double r_decay, double gaus_amplitude) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double dx2=gridSize().x*gridSize().x;
    double dy2=gridSize().y*gridSize().y;
    double rd2=r_decay*r_decay;
    uniform_real_distribution<double> randUR;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=dx2*(i-0.5*Nx)*(i-0.5*Nx)+dy2*(j-0.5*Ny)*(j-0.5*Ny);
            f_host[0][idx]=gaus_amplitude*(exp(-r2/rd2))+0.00000001*(randUR(rng)-0.5);
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};

// -----------------------------------------------------------------------
void Field::initFieldSin(double sin_amplitude=1, int sin_period=1, double sin_phase=0) {
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    double r2m=0.25*Nx*Nx+0.25*Ny*Ny;
    uniform_real_distribution<double> randUR;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            double r2=(i-0.5*Nx)*(i-0.5*Nx)+(j-0.5*Ny)*(j-0.5*Ny);
            f_host[0][idx]=sin_amplitude*sin(4*M_PI*sin_period*i/Nx)*sin(2*M_PI*sin_period*j/Ny)+0.0001*(randUR(rng)-0.5);           
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};


// =================================================================

#endif
