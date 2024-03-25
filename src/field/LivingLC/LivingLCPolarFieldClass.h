#ifndef LIVINGLCPOLARFIELDCLASS_H
#define LIVINGLCPOLARFIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include <map>
#include "../fieldclass.cu"


using namespace std; 


// ===============================================================
class LivingLCPolarField :public Field {
    
public:
          
    // ===========================================================    
    Field* ptr_cplus;
    Field* ptr_cminus;
    Field* ptr_Pxx;
    Field* ptr_theta;
    Field* ptr_theta_old;    
    Field* ptr_px;
    Field* ptr_py;
    Field* ptr_flip;        
    
    // ===========================================================
    // Methods
    
    // ===========================================================
    // Constructor
    LivingLCPolarField () {};
    LivingLCPolarField (Mesh* mesh_ptr_t, string name_t);
    LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t);
    LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t);
    LivingLCPolarField (Mesh* mesh_ptr_t, string name_t, int priority_t, string init_cond_t, string boun_cond_t, string expo_data_t);
    
    // Get velocity field
    void initPolarField();
    void postProcessing(int i_field);
    
    // ===========================================================
};



#endif
