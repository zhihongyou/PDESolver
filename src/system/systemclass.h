#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../field/fieldclass.cu"
#include "../field/incompressible-flow/incompressibleflowclass.cu"


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

    IncompressibleFlow* incomFlow_ptr;

    System () {
        setFFuncMap ();
    };
    

//................Methods ................................................
    // Constructor, sets the default values of parameters
    // System ();
    // Print system information.
    void printSysInfo ();
    void addIncompressibleFlow (IncompressibleFlow* incomFlow_ptr_t);
};



#endif
