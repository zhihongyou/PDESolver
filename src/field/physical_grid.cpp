#include "physical_grid.h"
#include <iostream> 
#include <vector>
#include <random>
#include <math.h> 



using namespace std; 



//------------------fixed boundary value------------------------------------------- 

std::vector<std::vector<double> > 
  Physical_Grid:: boundary_binding(
                                    std::vector<std::vector<double> > target_field,// 2D field that require binding 
                                    double fixed_boundary_value){ // fixed boundary value 

      // define boundary value 
            // looping columns to asign values at the right and left boundaries. 
            for (int i=0 ; i<(size_x);i++){
               for (int j=0 ; j<(size_y);j++){
                 field_with_bc[i+shifts][j+shifts]=target_field[i][j]; //alsign inner region 
               };
            };

      // for boundary region 
            // column  
            for (int i=0 ; i<(shifts);i++){ // counting the boundary 
               for (int j=0 ; j<net_size_y;j++){  //alsign boundary values 
                // left locations 
                 field_with_bc[i][j]= fixed_boundary_value; 
                 // right locations
                 field_with_bc[i+shifts+size_x][j]= fixed_boundary_value;
               };
            };
         

            // row  
            for (int i=0 ; i<(shifts);i++){ // counting the boundary 
               for (int j=0 ; j<net_size_x;j++){
                // top locations 
                 field_with_bc[j][i]= fixed_boundary_value;
                 // bottom locations
                 field_with_bc[j][i+shifts+size_x]= fixed_boundary_value;
               };
            };


            // calculate the difference between the boundary and the points on the boundary. This will be used to calculate gradient 
                 // column 
                for (int k=0 ; k<correlator;k++){
                  for (int j=0 ; j<net_size_y;j++){
                  difference_function_x[j][k]=fixed_boundary_value-field_with_bc[shifts][j];
                  difference_function_x[j][correlator+k]=fixed_boundary_value-field_with_bc[net_size_x-shifts-1][j];
                  };
                };


            // calculate the difference between the boundary and the points on the boundary. 
                 // row 
                for (int k=0 ; k<correlator;k++){
                  for (int j=0 ; j<net_size_y;j++){
                  difference_function_y[j][k]=fixed_boundary_value-field_with_bc[j][shifts];
                  difference_function_y[j][correlator+k]=fixed_boundary_value-field_with_bc[j][net_size_x-shifts-1];
                  };
                };
      //Updated current grid

    return field_with_bc;
};


std::vector<std::vector<double> >
Physical_Grid:: boundary_unbinding(
                                   std::vector<std::vector<double> > target_field){
std::vector<std::vector<double> >  grid_reduction(size_x , std::vector<double> (size_y,0) ); // define the reduced field
            for (int i=0 ; i<(size_x);i++){
               for (int j=0 ; j<(size_y);j++){
                 grid_reduction[i][j]=target_field[i+shifts][j+shifts];
               };
            };
      return grid_reduction;
};


  std::vector<std::vector<double> > 
  Physical_Grid:: get_grid(
                          std::vector<std::vector<double> > &current_grid,bool display){
      if (display=true){
        for (int i =0; i<current_grid.size();i++){ // row 
          cout << "-------------" <<endl;
          for (int j=0; j<current_grid[0].size();j++){ //column 
                  cout <<  current_grid[i][j]<< '\t';
          };
          cout<<endl;
        };
      };

  return current_grid;
  };

