#include "Field_v2.h"
#include <iostream> 
#include <vector>
#include <fstream>
#include <random>
#include <math.h> 

using namespace std; 
//----------------- function  inital field value Gaussian function------------------------------------------- 
// for a given mesh dx,dy, calculate its field value at each mesh point
std::vector<std::vector<double> > 
Fields:: gaussian_field(
                        double mean,double sigma){
        double current_x,current_y,estimate_field; // define variables 
        std::vector<std::vector<double> >  gaussian_out(size_x , std::vector<double> (size_y,0) );

          // looping both x and y to create mesh_grid for x,y and Gaussian values. 
          for (int i =0; i<size_x;i++){ // row 
          current_x= physical_x0+i*physical_dx; // define current estimate value for x
          
              for (int j=0; j<size_y;j++){ //column 
                  current_y=physical_y0+j*physical_dy;// define current estimate value for x
                  projected_physical_grid_x[i][j]=current_x; // set gridmesh values for both x y. 
                  projected_physical_grid_y[i][j]=current_y;   // set gridmesh values for both x y.            

                  //define the Gaussian function that estimate the field value 
                  estimate_field=exp(-(current_x*current_x+current_y*current_y));
                  gaussian_out[i][j]=estimate_field;
                };
            };
        return gaussian_out;

};

//------------------finite difference laplacian-------------------------------------------- 
std::vector<std::vector<double> >  
Fields:: calculate_field_laplacian(
                            std::vector<std::vector<double> > target_field){
// The method calculate the laplacian for the mesh_grid in the class object. One must run the function_map before they run this. 
        for (int i =0; i<size_x;i++){ // row 
            for (int j=0; j<size_y;j++){ //column 
                // For now, this is only supported for equal spacing
                // setting the boundary grad equal to the one inside  (or zero, which is currently using)
                // The boundary is replaced with the near elements (or zero), so the grad matrix has the same dimension as the mesh matrix. 
                
                // Condition to see if we are on the mesh boundary.
               field_laplacian[i][j]=(target_field[i+shifts+1][j+shifts]+target_field[i+shifts-1][j+shifts]+target_field[i+shifts][j+shifts+1]+target_field[i+shifts][j+shifts-1]-4*target_field[i+shifts][j+shifts])/(physical_dx*physical_dy);
              // cout<< physical_dx*physical_dy<<endl;
              };
          };

      return field_laplacian;
};

    // get the difference field column direction 
    std::vector<std::vector<double> >  Fields:: get_difference_boundary_grid_x(){
        return Physical_Grid :: get_grid(difference_function_x);
    };
    // get the difference field row direction 
    std::vector<std::vector<double> > Fields::  get_difference_boundary_grid_y(){
        return Physical_Grid :: get_grid(difference_function_y);
};







