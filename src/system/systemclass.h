#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>

using namespace std; 


//................Class .................................................

class System {
    // define parameters for a mesh. Usually, we only consider a rectanguar mesh since most of the problem can be comformally mapped onto it.
    
    // protected so that one can access on different files.
    // private:
    
    
    // one can call outside of the class.
    public: 

    // Dimension of the system.
    int Dimension=1;
    // Size of the system in different dimensions.
    double Lx=2, Ly=1, Lz=1;
    // Number of grids in different dimensions.
    int Nx=64, Ny=1, Nz=1;
    // Number of boundary grids in different dimensions.
    int Nbx=2, Nby=2, Nbz=2;
    // Spatial resolution in different directions.
    double Dx, Dy, Dz;

    
    // Constructor, sets the default values of parameters
    System (int dimension=1, double lx=1, double ly=1, double lz=1, int nx=64, int ny=1, int nz=1, int nbx=2, int nby=2, int nbz=2) {
        Dimension=dimension;
        Lx=lx;
        Ly=ly;
        Lz=lz;
        Nx=nx;
        Ny=ny;
        Nz=nz;
        Nbx=nbx;
        Nby=nby;
        Nbz=nbz;
        Dx=Lx/Nx;
        Dy=Ly/Ny;
        Dz=Lz/Nz;
    };


//................Methods ................................................

    
};



#endif
