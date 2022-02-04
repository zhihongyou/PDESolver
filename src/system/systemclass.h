#ifndef RAW_GRID
#define RAW_GRID
#include <iostream> 
#include <vector>
#include <string>


using namespace std; 



//................Class ...................................................................

class Raw_Mesh{
    // define parameters for a mesh. Usually, we only consider a rectanguar mesh since most of the problem can be comformally mapped onto it. 
    protected:
    // protected so that one can access on different files. 
      double dx_raw, dy_raw;  // this is for the raw grid. One need to use the dx dy (projected from grid value assignment)
      double unit_density;


    public: // one can call outside of the class. 
    
    // Define Constructors (length x ,length y, density)
   Raw_Mesh(double raw_x_size=60,double raw_y_size=60, double unit_density=60);
   // this function define the range in 2D, with given linear density (same for x and y)




//................Methods ...................................................................
// getting mesh data output 

    
};

#endif