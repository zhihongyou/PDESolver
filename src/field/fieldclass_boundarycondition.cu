#ifndef FIELDCLASS_BOUNDARYCONDITION_CU
#define FIELDCLASS_BOUNDARYCONDITION_CU

#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h>
#include <cmath>
#include <sstream>
#include <iomanip>


using namespace std;


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


// =================================================================

#endif
