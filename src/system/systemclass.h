#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../mesh/meshclass.h"
#include "../mesh/meshclass.cpp"
#include "../field/fieldclass.h"
#include "../field/fieldclass.cpp"
// #include "../utility/vector3class.h"

using namespace std; 

//................Class .................................................

class System {
    // private:
    
    // -------------------------------------------------------------------
    // one can call outside of the class.
    public:
    // Name of system.
    std::string name="system";    
    // Pointer to mesh.
    Mesh* mesh_ptr;
    // Pointer to all fields.
    std::list<Field*> field_ptrs;
    

//................Methods ................................................
    // Constructor, sets the default values of parameters
    // System ();
    // Print system information.
    void printSysInfo ();
    
};



#endif
