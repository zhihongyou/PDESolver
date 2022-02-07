#include "systemclass.h"
#include <iostream> 
#include <vector>

using namespace std; 


// ======================================================================
void System::printSysInfo () {
    cout << "System's info:" <<endl;
    // Printing grid info.
    if (mesh_ptr != nullptr) {
        cout <<"This is a "<<(*mesh_ptr).space_dim<<" dimensional system.\n";
        switch ((*mesh_ptr).space_dim) {
            case 1:
                cout <<"The box size is: "<<(*mesh_ptr).box_size.x<<".\n";
                cout <<"There are "<<(*mesh_ptr).grid_number.x<<" grid points, ";
                cout <<"with "<<(*mesh_ptr).grid_number_boun.x<<" layers of boundary grids.\n";
            case 2:
                cout <<"The box size is "<<(*mesh_ptr).box_size.x;
                cout <<"*"<<(*mesh_ptr).box_size.y<<".\n";
                cout <<"There are "<<(*mesh_ptr).grid_number.x<<"*";
                cout <<(*mesh_ptr).grid_number.y<<" grid points, ";
                cout <<"and "<<(*mesh_ptr).grid_number_boun.x<<"*";
                cout <<(*mesh_ptr).grid_number_boun.y<<" layers of boundary grids.\n";
            case 3:
                cout <<"The box size is "<<(*mesh_ptr).box_size.x;
                cout <<"*"<<(*mesh_ptr).box_size.y;
                cout <<"*"<<(*mesh_ptr).box_size.z<<".\n";
                cout <<"There are "<<(*mesh_ptr).grid_number.x<<"*";
                cout <<(*mesh_ptr).grid_number.y<<"*";
                cout <<(*mesh_ptr).grid_number.z<<" grid points, ";
                cout <<"and "<<(*mesh_ptr).grid_number_boun.x<<"*";
                cout <<(*mesh_ptr).grid_number_boun.y<<"*";
                cout <<(*mesh_ptr).grid_number_boun.z<<" layers of boundary grids.\n";
        };
        // cout <<(*mesh_ptr).
    } else {
        cout <<"The system contains NO mesh!" <<endl;
    };
    
    // Printing field info.
    if (field_ptrs.size()==0) {
        cout << "There are no field in this system."<<endl;
    } else {
        cout << "There are " <<field_ptrs.size()<<" fields: ";
        for (auto iter = field_ptrs.begin(); iter != field_ptrs.end(); ) {        
            std::cout << (*iter)->name;
            if ( ++iter == field_ptrs.end() ) {
                cout<<"." <<endl;
            } else {
                cout <<",";
            };
        };
    };
    
};


// System:: System(double dx,double dy, double dz){
//     // define unit steps that previously defined in the header file. This is replaced by the updated version in grid_access. 
//     dx=Lx/Nx;
//     dy=Ly/Ny;
//     dz=Lz/Nz;
// };


