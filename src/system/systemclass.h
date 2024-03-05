#ifndef SYSTEMCLASS_H
#define SYSTEMCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../field/fieldClassHeaders.h"
// #include "../field/fieldclass.cu"
// #include "../field/IncompFlow/IncompFlowClass.cu"
// #include "../field/constant-field/constantfieldclass.cu"
// #include "../field/LaplaceNFEqField/LaplaceNFEqFieldClass.cu"

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
    
    IncompFlow* incompFlow_ptr;
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
    void addField (Field* f_ptr_t);
    void addIncompFlow (IncompFlow* incompFlow_ptr_t);
    void addConstantField (ConstantField* constF_ptr_t);
    void addLivingLC (LivingLC* LivingLC_ptr_t);    
};



#endif
