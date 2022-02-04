#ifndef FIELD_H
#define FIELD_H
#include <iostream> 
#include <vector>
#include <string>
#include "physical_grid.h"
#include "physical_grid.cpp"
using namespace std; 

// The field on the mesh grid
class Fields: public Physical_Grid {
    // define parameters 
  // private: // has nothing to do with the parent class private.

    protected:

    // define the gradient function  
    std::vector<std::vector<double> > field_laplacian;
    std::vector<std::vector<double> > next_time_field;

    public:
    // Initialize finite difference (i.e. gradient)
    Fields(         double raw_x_size,
                   double raw_y_size, 
                   int unit_density,
                   double x0,double y0,
                   double x1,double y1,
                   int boundary_width=2,
                   int internal_boundaries=0)
    :Physical_Grid( raw_x_size,
                   raw_y_size, 
                   unit_density,
                   x0, y0,
                   x1, y1,
                   boundary_width,
                   internal_boundaries){


          // Define the vector: std::vector<std::vector<double> > V(ROWS, std::vector <double> (COLUMNS, DEFAULT_VALUE))  
          // Define the laplacian vector: 
          std::vector<std::vector<double> >  grad_field_out(size_x , std::vector<double> (size_y,0) );
          // point the public into protected to call
          field_laplacian=grad_field_out;
          // total order shift when calculate the laplacian. 
          shifts=boundary_width;
          // Each variable on the protected, one must defined and assign their value here. 
          std::vector<std::vector<double> > field_out_next(size_x , std::vector<double> (size_y,0) );         
          next_time_field=field_out_next;
  
    };

//................Methods ...................................................................

    //------------------calculate laplacian-------------------------------------------- 
        // define laplacian with equal spacing
    std::vector<std::vector<double> > calculate_field_laplacian(std::vector<std::vector<double> > target_field);


        // Map the grid onto gaussian function, with given function x-y range. 

    std::vector<std::vector<double> >  gaussian_field(double mean=0,double sigma=1);

    //------------------call original field grid------------------------------------------- 
     std::vector<std::vector<double> > get_original_field(bool display=false);


    //------------------get the mesh grid------------------------------------------- 
     std::vector<std::vector<double> > get_meshgrid_x(bool display=false);    
     std::vector<std::vector<double> > get_meshgrid_y(bool display=false);    

    //------------------call laplacian field grid------------------------------------------- 
     std::vector<std::vector<double> > get_laplacian_field(bool display=false);

    //------------------call laplacian field grid------------------------------------------- 
    std::vector<std::vector<double> >   get_difference_boundary_grid_x();
    std::vector<std::vector<double> >   get_difference_boundary_grid_y();

    //------------------export_current field ------------------------------------------ 
    void out_field(int id=0,bool display=true,string path="output/field.csv");
      
    //------------------export_current field grad ------------------------------------------ 
    void out_field_grad(int id=0,bool display=true,string path="output/field_grad.csv");
    void  out_meshgird(int id=0,bool display=true,string path_x="output/meshgrid_x.csv",string path_y="output/meshgrid_y.csv");

    void initial_file_loc(string path);
    //------------------update current field ------------------------------------------ 
    void  update_field_with_bc(std::vector<std::vector<double> >  updated_field );



    //------------------Stepping Function------------------------------------------ 
    std::vector<std::vector<double>> newton_stepping(std::vector<std::vector<double>> gradient_terms,double dt=0.1,double D=1); // input the time interval 



};
#endif