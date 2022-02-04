#include "Field_v2.h"
#include <iostream> 
#include <vector>
#include <fstream>


using namespace std; 


//------------------update current mesh_grid to the supplied field-------------------------------------------- 

void  Fields:: update_field_with_bc(std::vector<std::vector<double> > updated_field ){
    field_with_bc=updated_field;
};

//------------------finite difference laplacian-------------------------------------------- 
std::vector<std::vector<double>>  
Fields::  newton_stepping(
                          std::vector<std::vector<double>> gradient_terms,double dt,double D){
// initial variable 
std::vector<std::vector<double> >  next_time_field_grid(size_x , std::vector<double> (size_y,0) ); 
// the new field is updated by the gradient term supplied
      for (int i =0; i<size_x;i++){ // row 
          for (int j=0; j<size_y;j++){ //column 
            next_time_field_grid[i][j]= field_with_bc[i+shifts][j+shifts]+dt*D*(gradient_terms[i][j]-field_with_bc[i+shifts][j+shifts]);
          };

      };
    return next_time_field_grid;

};
