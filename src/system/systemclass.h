#ifndef SYSTEMCLASS
#define SYSTEMCLASS
#include <iostream> 
#include <vector>
#include <string>

using namespace std; 


//................Class .................................................

class System{
    // define parameters for a mesh. Usually, we only consider a rectanguar mesh since most of the problem can be comformally mapped onto it.
    
    // protected so that one can access on different files.
    protected:
    
    // Size of the system in different dimensions.
    real Lx, Ly, Lz;
    // Number of grids in different dimensions.
    int Nx, Ny, Nz;
    // Number of boundary grids in different dimensions.
    int Nbx, Nby, Nbz;
    // Spatial resolution in different directions.
    real dx, dy, dz;
    

    public: // one can call outside of the class. 
    
    // Define Constructors (length x ,length y, density)
    System(real Lx=1, real Ly=1, real Lz=1, int Nx=60, int Ny=1, int Nz=1);


//................Methods ................................................

    
};

#endif
