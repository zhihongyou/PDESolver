#ifndef MESHCLASS_H
#define MESHCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../utility/vector3class.h"

using namespace std; 

//................Class .................................................

class Mesh {
    // private:
    
    // -------------------------------------------------------------------
    // one can call outside of the class.
    public: 

    // Dimension of the system.
    int space_dim=1;
    // Size of the system in different dimensions.
    Vector3<double> box_size = Vector3<double>(1,1,1);
    // Number of grids in different dimensions.
    Vector3<int> grid_number = Vector3<int>(64,1,1);
    // Number of boundary grids in different dimensions.
    Vector3<int> grid_number_boun = Vector3<int>(2,0,0);
    // Spatial resolution in different directions.
    Vector3<double> grid_size = Vector3<double>(1,1,1);
    
    
    // -------------------------------------------------------------------
    // Constructor, sets the default values of parameters
    Mesh (int dimension=1) {
        space_dim=dimension;
        switch(space_dim) {
            case 2:
                grid_number.y=64;
                grid_number_boun.y=2;
            case 3:
                grid_number.y=64;
                grid_number_boun.y=2;
                grid_number.z=64;
                grid_number_boun.z=2;
        };
        grid_size.x=box_size.x/grid_number.x;
        grid_size.y=box_size.y/grid_number.y;
        grid_size.z=box_size.z/grid_number.z;
    };


//................Methods ................................................
    // -------------------------------------------------------------------
    void setGridSize () {
        grid_size.x=box_size.x/grid_number.x;
        grid_size.y=box_size.y/grid_number.y;
        grid_size.z=box_size.z/grid_number.z;
    };

    // -------------------------------------------------------------------
    void setSpaceDim (int dimension) {
        space_dim=dimension;
        switch(space_dim) {
            case 1:
                grid_number.y=1;
                grid_number.z=1;
            case 2:                
                grid_number.z=1;
        };
        setGridSize();
    };
    
    // -------------------------------------------------------------------
    void setBoxSize (double lx, double ly, double lz) {
        box_size.x=lx;
        switch(space_dim) {
            case 2:
                box_size.y=ly;
            case 3:
                box_size.y=ly;
                box_size.z=lz;
        };
        setGridSize();        
    };

    // -------------------------------------------------------------------
    void setGridNumber (int nx, int ny, int nz) {
        grid_number.x=nx;
        switch(space_dim) {
            case 2:
                grid_number.y=ny;
            case 3:
                grid_number.y=ny;
                grid_number.z=nz;
        };
        setGridSize();        
    };

    // -------------------------------------------------------------------
    void setGridNumberBoun (int nx, int ny, int nz) {
        grid_number_boun.x=nx;
        switch(space_dim) {
            case 2:
                grid_number_boun.y=ny;
            case 3:
                grid_number_boun.y=ny;
                grid_number_boun.z=nz;
        };
    };

    // -------------------------------------------------------------------
    int getGridNumberAll () {
        return (grid_number.x+2*grid_number_boun.x)*(grid_number.y+2*grid_number_boun.y)*(grid_number.z+2*grid_number_boun.z);
    };

    // -------------------------------------------------------------------
    int getGridNumberAllBulk () {
        return grid_number.x*grid_number.y*grid_number.z;
    };

    // ==================================================================
};



#endif
