#ifndef FIELDCLASS_H
#define FIELDCLASS_H
#include <iostream> 
#include <vector>
#include <string>
#include "../system/systemclass.h"

using namespace std; 


// The field on the mesh grid
class Field {
    // define parameters 
  // private: // has nothing to do with the parent class private.

    protected:

    
    public:

    char name="fieldTemp";
    int rank=0;
    int priority=1;
    char location="both";
    char bounCond="none";
    char export="on";

    double * fcpu;
    double * fgpu;
    double * dtf;
    double * fcopy1;
    double * fcopy2;
    double * fcopy3;
    double * fcopy4;
    double * dtfcopy1;
    double * dtfcopy2;
    double * dtfcopy3;
    double * dtfcopy4;

};

#endif
