#include "systemclass.h"
#include <iostream> 
#include <vector>

using namespace std; 


System:: System(double dx,double dy, double dz){
    // define unit steps that previously defined in the header file. This is replaced by the updated version in grid_access. 
    dx=Lx/(Nx+0.0);
    dy=Ly/(Ny+0.0);
    dz=Lz/(Nz+0.0);
};


