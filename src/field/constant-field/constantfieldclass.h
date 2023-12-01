#ifndef CONSTANTFIELDCLASS_H
#define CONSTANTFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// =======================================================================
class ConstantField {    
    
    public:

    string name;
    int priority=0;
    string boun_cond = "periodic";
    string init_cond = "";
    int rank=0;
    string location = "both";
    string expo_data = "on";
    string equation = "";
    Mesh* mesh_ptr=NULL;
    
    Field one;
    Field zero;
    Field pi;
    
    // ===================================================================
    // Methods
    
    // ===================================================================
    // Constructor
    ConstantField (Mesh* mesh_ptr_t, string name_t);
    void initFields ();
    void setFieldValues ();
    void setFieldProperties (Field* field_ptr, string field_name, int priority_t);
    // ===================================================================
};



#endif
