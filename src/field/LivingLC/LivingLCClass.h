#ifndef LIVINGLCCLASS_H
#define LIVINGLCCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "LivingLCPolarFieldClass.cu"


using namespace std; 


// ===============================================================
class LivingLC {
    
    public:

    string name;
    int priority=0;
    string boun_cond = "periodic";
    string init_cond = "sin";
    int rank=0;
    string location = "both";
    string expo_data = "on";
    string equation = "";
    Mesh* mesh_ptr=NULL;
    string FDMScheme;

    Field cplus;
    Field cminus;
    Field Pxx;
    LivingLCPolarField Pxy;
    Field theta;
    Field theta_old;
    Field px;
    Field py;
    Field flip;
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    LivingLC (Mesh* mesh_ptr_t, string name_t);
    LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t);
    LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    LivingLC (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);

    // Fields
    void initFields ();
    void setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t);
        
    // ===========================================================
};



#endif
