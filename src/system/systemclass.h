#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>
// #include "../field/fieldclass.h"
#include "../field/fieldclass.cu"


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
    Mesh* mesh_ptr=NULL;
    // Pointer to all fields.
    std::vector<Field*> field_ptrs;
    

//................Methods ................................................
    // Constructor, sets the default values of parameters
    // System ();
    // Print system information.
    void printSysInfo ();
    
};



#endif
