#ifndef INCOMPFLOWCLASS_H
#define INCOMPFLOWCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"
#include "IncompFlowOmegaFieldClass.cu"


using namespace std; 


// ===============================================================
class IncompFlow {
    
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
    
    IncompFlowOmegaField omega;    // omega=2*vorticity
    Field vx;
    Field vy;    
    Field phi;                  // stream function    
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    IncompFlow (Mesh* mesh_ptr_t, string name_t);
    IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t);
    IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    IncompFlow (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);

    // Fields
    void initFields ();

    void setFieldProperties (Field* field_ptr, string field_name, int priority_t, string init_cond_t);
    void setFieldValues ();
    void setVValues();
    void setOmegaValues ();
        
    // ===========================================================
};



#endif
