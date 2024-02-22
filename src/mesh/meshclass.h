#ifndef MESHCLASS_H
#define MESHCLASS_H

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif

#include "cuda_runtime.h"
#include "cuda.h"
#include "device_launch_parameters.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <cstdio>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <random> // The header for the generators.
#include <iomanip>
#include <iostream> 
#include <vector>
#include <string>
#include "../utility/Vector3/vector3class.h"

#define real double
#define Pi 3.1415926535897932384626433832795
#define Zero 0
typedef double2 Complex;

using namespace std;


// ======================================================================
// Data wrapper of Mesh class. This allows storage of
struct MeshData {
    int space_dim=1;
    // Size of the system in different dimensions.
    Vector3<double> box_size = Vector3<double>(64,1,1);
    // Number of grids in different dimensions.
    Vector3<int> grid_number = Vector3<int>(64,1,1);
    // Number of boundary grids in different dimensions.
    Vector3<int> grid_number_boun = Vector3<int>(3,0,0);
    // Spatial resolution in different directions.
    Vector3<double> grid_size = Vector3<double>(1,1,1);
};


// ======================================================================
// Mesh class.
class Mesh {
    
// -------------------------------------------------------------------
public: 

    // Mesh data on CPU;
    MeshData host;
    // Mesh data on GPU;
    MeshData* dev_ptr;    
    
    // -------------------------------------------------------------------
    // Constructor, sets the default values of parameters
    Mesh (int dimension=1) {
        host.space_dim=dimension;
        switch(host.space_dim) {
            case 2:
                host.box_size.y=64;
                host.grid_number.y=64;                
                host.grid_number_boun.y=3;
                break;
            case 3:
                host.box_size.y=64;
                host.box_size.z=64;
                host.grid_number.y=64;
                host.grid_number_boun.y=3;
                host.grid_number.z=64;
                host.grid_number_boun.z=3;
        };
        host.grid_size.x=host.box_size.x/host.grid_number.x;
        host.grid_size.y=host.box_size.y/host.grid_number.y;
        host.grid_size.z=host.box_size.z/host.grid_number.z;
        // Create mesh data on GPU, and copy from CPU to GPU.
        cudaMalloc(&dev_ptr, sizeof(struct MeshData));
        cudaMemcpy(dev_ptr, &host, sizeof(struct MeshData),cudaMemcpyHostToDevice);
    };

    
// ======================================================================
// Methods
    
    // -------------------------------------------------------------------
    // Copy mesh data from GPU to CPU
    void updateMeshDev () {
        cudaMemcpy(dev_ptr, &host, sizeof(struct MeshData),cudaMemcpyHostToDevice);
    };

    // -------------------------------------------------------------------
    // Copy mesh data from GPU to CPU
    void updateMeshHost () {
        cudaMemcpy(&host, dev_ptr, sizeof(struct MeshData),cudaMemcpyHostToDevice);
    };
        
    // -------------------------------------------------------------------
    void setGridSize () {
        host.grid_size.x=host.box_size.x/host.grid_number.x;
        host.grid_size.y=host.box_size.y/host.grid_number.y;
        host.grid_size.z=host.box_size.z/host.grid_number.z;
        updateMeshDev();
    };

    // -------------------------------------------------------------------
    void setSpaceDim (int dimension) {
        host.space_dim=dimension;
        switch(host.space_dim) {
            case 1:
                host.grid_number.y=1;
                host.grid_number.z=1;
                break;
            case 2:                
                host.grid_number.z=1;
        };
        setGridSize();
    };
    
    // -------------------------------------------------------------------
    void setBoxSize (double lx, double ly, double lz) {
        host.box_size.x=lx;
        switch(host.space_dim) {
            case 2:
                host.box_size.y=ly;
                break;
            case 3:
                host.box_size.y=ly;
                host.box_size.z=lz;
        };
        setGridSize();        
    };

    // -------------------------------------------------------------------
    void setGridNumber (int nx, int ny, int nz) {
        host.grid_number.x=nx;
        switch(host.space_dim) {
            case 2:
                host.grid_number.y=ny;
                break;
            case 3:
                host.grid_number.y=ny;
                host.grid_number.z=nz;
        };
        setGridSize();        
    };

    // -------------------------------------------------------------------
    void setGridNumberBoun (int nx, int ny, int nz) {
        host.grid_number_boun.x=nx;
        switch(host.space_dim) {
            case 2:
                host.grid_number_boun.y=ny;
                break;
            case 3:
                host.grid_number_boun.y=ny;
                host.grid_number_boun.z=nz;
        };
        updateMeshDev();
    };

    // -------------------------------------------------------------------
    int GridNumberAll () {
        return (host.grid_number.x+2*host.grid_number_boun.x)*(host.grid_number.y+2*host.grid_number_boun.y)*(host.grid_number.z+2*host.grid_number_boun.z);
    };

    // -------------------------------------------------------------------
    int GridNumberAllBulk () {
        return host.grid_number.x*host.grid_number.y*host.grid_number.z;
    };

    // -------------------------------------------------------------------

    // ==================================================================
    
};



#endif
