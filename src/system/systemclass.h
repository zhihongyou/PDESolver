#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../field/fieldclass.cu"
#include "../field/incompressible-flow/incompressibleflowclass.cu"
#include "../field/constant-field/constantfieldclass.cu"

using namespace std; 

//................Class .................................................

class System {
    // private:
    
    // ---------------------------------------------------------------
    // one can call outside of the class.
    public:
    // Name of system.
    std::string name="system";
    std::string dire_expo="./";
    // Pointer to mesh.
    Mesh* mesh_ptr=NULL;
    // Pointer to all fields.
    std::vector<Field*> field_ptrs;
    
    IncompressibleFlow* incomFlow_ptr;
    ConstantField* constF_ptr;
    
    System () {
        setFFuncMap ();
    };

    System (string dire_expo_t) {
        dire_expo=dire_expo_t;
        setFFuncMap ();
    };
    

//................Methods ................................................
    // Constructor, sets the default values of parameters
    // System ();
    // Print system information.
    void printSysInfo ();
    void addIncompressibleFlow (IncompressibleFlow* incomFlow_ptr_t);
    void addConstantField (ConstantField* constF_ptr_t);
    void addField (Field* f_ptr_t);
    
};



#endif
