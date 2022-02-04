#include "raw_grid.h"
#include <iostream> 
#include <vector>

using namespace std; 

Raw_Mesh:: Raw_Mesh(double raw_x_size,double raw_y_size, double unit_density){
      // define unit steps that previously defined in the header file. This is replaced by the updated version in grid_access. 
        dx_raw=raw_x_size/unit_density;
        dy_raw=raw_y_size/unit_density;
        unit_density=unit_density;


};


