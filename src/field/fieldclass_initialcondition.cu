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


// -----------------------------------------------------------------------
void Field::initFieldConst(double f_value) {    
    int Nx=gridNumber().x;
    int Ny=gridNumber().y;
    int Nbx=gridNumberBoun().x;
    int Nby=gridNumberBoun().y;
    uniform_real_distribution<double> randUR;
    for (int j=0; j<Ny; j++) {
        for (int i=0; i<Nx; i++) {
            int idx=(j+Nby)*(Nx+2*Nbx)+i+Nbx;
            f_host[0][idx]=f_value+0.0001*(randUR(rng)-0.5);
        };
    };
    applyBounCondPeriCPU(f_host[0]);
};

// -----------------------------------------------------------------------
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
            f_host[0][idx]=gaus_amplitude*(exp(-r2/rd2))+0.0001*(randUR(rng)-0.5);
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
