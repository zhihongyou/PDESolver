#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "src/system/systemclass.h"
#include "src/system/systemclass.cpp"
// #include "src/utility/dataType.h"

using namespace std;

// =======================================================================
void testFunction(double* u, int N) {
    for (int i=0; i<N; i++) {
        u[i]=i;
    };
};

// Field* newField(System mySys, std::string f_name, int f_rank, int f_priority, std::string f_boun_cond, std::string f_init_cond, std::string f_export) {
//     Field field_temp(mySys, f_name, f_rank, f_priority, f_boun_cond, f_init_cond, f_export);
//     return &field_temp;
//     // Field phi(mySys, "phi", 0, 0, "periodic", "ones", "on");
// };


// =======================================================================
int main() {    

    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field phi(mesh, "phi", 0, 0, "periodic", "ones", "on");
    // Add the new field to the system.
    mySys.field_ptrs.push_back(&phi);
    Field rho(mesh, "rho", 0, 0, "periodic", "ones", "on");    
    mySys.field_ptrs.push_back(&rho);

    // Printing system's information.
    mySys.printSysInfo();
    
    
    
    return 0;
};


