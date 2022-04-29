#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "cuda_runtime.h"
#include "cuda.h"
#include "device_launch_parameters.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
// #include "src/mesh/meshclass.cuh"
// #include "src/field/fieldclass.cu"
// #include "src/system/systemclass.cpp"
#include "src/evolver/evolverclass.h"
#include "src/evolver/evolverclass.cu"
// #include "src/utility/finiteDifference.h"
// #include "src/utility/finiteDifferenceCentralO2Isotropic.h"

using namespace std;

__global__ void testGPU (rhsPtrs rhs_ptrs_dev, int a) {
    rhs_ptrs_dev.num_terms[0]+=a;
    // rhs_ptrs_dev.
    // f_dev[0][10]=1000;
};

__global__ void testGPU1 (double* prefactors) {
    // rhs_ptrs_dev.num_terms[0]+=0;
    prefactors[0]+=11;
    // rhs_ptrs_dev.num_funcs_1term[0]+=12;
};


// ======================================================================
int main() {
    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.    
    Field mu(&mesh, "mu", 0, 1, "periodic", "Gaussian", "on");
    Field phi(&mesh, "phi", 0, 0, "periodic", "Gaussian", "on");

    // Set field equations.
    mu.setRhsTerms({
        {2.5,{{"laplace",&phi}},"explicit"},
        {-1.5,{{"1",&phi},{"laplace",&phi},{"1",&phi}},"explicit"}
    });
    phi.setRhsTerms({
        {2.5,{{"laplace",&phi}},"explicit"}
    });
    
    // Add fields to the system.
    // mySys.field_ptrs.push_back(&mu);
    mySys.field_ptrs.push_back(&phi);
    // Print system information.
    // mySys.printSysInfo();
    
    // Creating an evolver:
    string device="cpu";
    Evolver evolver(&mySys,0,1000,0.001,100,"EulerForward",device);    
    evolver.run();

    // ------------------------------------------------------------
    // Initiate function calls and pointers.
    // evolver.initEvolver();
    // cout<<"Operator is: " <<phi.f_funcs_rhs[0].f_operator<<endl;
    // std::cout<< phi.rhs_ptrs_host.num_terms[0] <<endl;
    // std::cout<< phi.rhs_ptrs_host.num_funcs_1term[0] <<endl; 
    // std::cout<< phi.rhs_ptrs_host.f_func_ptrs[0] <<endl; 
    // std::cout<< phi.laplace <<endl;
    
    // cout <<"No. terms: " << phi.rhs_ptrs_host.num_terms[0] <<endl;
    // testGPU<<<1,1>>>(phi.rhs_ptrs_dev,1);
    // cudaMemcpy(phi.rhs_ptrs_host.num_terms, phi.rhs_ptrs_dev.num_terms, 2*sizeof(int),cudaMemcpyDeviceToHost);
    // cout <<"No. terms now: " << phi.rhs_ptrs_host.num_terms[0] <<endl;

    // cout <<"No. terms: " << phi.rhs_ptrs_host.num_terms[0] <<endl;
    // cout <<"Function address: " << phi.rhs_ptrs_host.f_func_ptrs[0] <<endl;
    // testGPU<<<1,1>>>(phi.rhs_ptrs_dev,1);
    // cudaMemcpy(phi.rhs_ptrs_host.num_terms, phi.rhs_ptrs_dev.num_terms, 2*sizeof(int),cudaMemcpyDeviceToHost);
    // cout <<"No. terms now: " << phi.rhs_ptrs_host.num_terms[0] <<endl;
    // cout <<phi.laplace<<endl;
        
    
    evolver.EulerForward();
    // ------------------------------------------------------------    
    // phi.export_conf("1",device,1);
    // phi.export_conf_any(phi.laplace[0],"laplace","0",device,1);
    // phi.export_conf_any(phi.rhs[0],"rhs","0",device,1);
    
    // Export initial configuration.
   
    // phi.export_conf("0",device,1);
    // cout <<"Export 0 successfully"<<endl;
    
    // evolver.EulerForward();

        
    // cudaMalloc(phi.laplace[0],phi.gridNumberAll()*sizeof(double));
    

    // cudaMemcpy(phi.f_host[0], phi.laplace[0], phi.gridNumberAll()*sizeof(double),cudaMemcpyDeviceToHost);
    // cudaMemcpy(phi.rhs_ptrs_host.prefactors, phi.rhs_ptrs_dev.prefactors, phi.rhs_ptrs_host.num_terms[0]*sizeof(double),cudaMemcpyDeviceToHost);
    // cout << phi.f_host[0][250]<<", "<<phi.rhs_ptrs_host.prefactors[0]<<endl;

    // cout <<phi.f_rhs[0][250]<<endl;
    
    // testGPU1<<<1,1>>>(phi.rhs_ptrs_dev.prefactors,phi.f_rhs);
    // cout<<phi.rhs_ptrs_host.prefactors[0]<<endl;
    // cudaMemcpy(phi.rhs_ptrs_host.prefactors, phi.rhs_ptrs_dev.prefactors, phi.rhs_ptrs_host.num_terms[0]*sizeof(double),cudaMemcpyDeviceToHost);
    // cout<<phi.rhs_ptrs_host.prefactors[0]<<endl;
    
    // evolver.EulerForward();
        

    // cout<<phi.f_funcs_rhs[0].f_operator<<endl;
    // cout<<&phi<<", "<<phi.f_funcs_rhs[0].field_ptr<<endl;
    
    // Run simulations.
    
    
    return 0;
};

