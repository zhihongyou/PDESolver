#ifndef EVOLVERCLASS
#define EVOLVERCLASS
#include <iostream> 
#include <vector>
#include <string>

using namespace std; 


//................Class .................................................

class Evolver{

    
    // protected so that one can access on different files.
    protected:

    // Which device to use to do computationally expensive parts.
    // "cpu" (default) or "gpu".
    char device;
    
    // Explicit or semi-implicit schemes for potential terms.
    // "explicit" (default) or "semiImplicit".
    char method;

    // Simulation time:
    // Simulation starts at timeStart and stops at timeStop, with an
    //   initial time step timeStep.
    // Fields are exported every timeExport.
    real timeStart, timeStop, timeStep, timeExport;
    

    public: // one can call outside of the class. 
    
    // Define Constructors ()
    Evolver(real timeStart=0, real timeStop=1, real timeStep=0.001, real timeExport=0.1);


//................Methods ................................................

    
};

#endif
