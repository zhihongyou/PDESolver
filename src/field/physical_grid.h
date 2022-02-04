#ifndef PHYSICAL_GRID
#define PHYSICAL_GRID
#include <iostream> 
#include <vector>
#include <string>
#include "raw_grid.h"
#include "raw_grid.cpp"

using namespace std; 

// The grid is the inheritance of the Mesh decoder

class Physical_Grid: public Raw_Mesh {
    // define parameters 
  // private: // has nothing to do with the parent class private.

    protected:
    int size_y,size_x, net_size_x, net_size_y,correlator;
    double physical_dx,physical_dy; 
    double physical_x0,physical_y0,physical_x1,physical_y1;

    // this defines the order of the gradient approximation. 
    int shifts;
   // std::vector<std::vector<double> > mesh_gird; // initializing mesh grid 2D. 
    std::vector<std::vector<double> > projected_physical_grid_x; // initializing mesh with x,y cood values. 
    std::vector<std::vector<double> > projected_physical_grid_y; // initializing mesh with x,y cood values. 


    std::vector<std::vector<double> >  field_with_bc;
  // define two point difference function 
    std::vector<std::vector<double> >  difference_function_x;
    std::vector<std::vector<double> >  difference_function_y;
          

    public:
    //inharent from Mesh_decoder (taking more parameters than the mesh)
    Physical_Grid(double raw_x_size,  // how many unit direction point in x (raw matrix scale)
                  double raw_y_size,  // how many unit direction point in y(raw matrix scale)
                  double unit_density, // density for unit direction point (raw matrix scale)
                  double x0,double y0, // first corner point (physical scale)
                  double x1,double y1, // second corner point (physical scale)
                  int boundary_width=2, // the boundary width (raw matrix scale)
                  int correlat=4, // difference matrix size for each end (raw matrix scale)
                  int internal_boundaries=0)
    :Raw_Mesh(raw_x_size, 
              raw_y_size,
              unit_density){

          //------------------Define new vars--------------------------------------------
            size_y= (int)(raw_y_size*unit_density); //column 
            size_x=(int)(raw_x_size*unit_density); //row 
            net_size_y =size_y+2*boundary_width; // total matrix size includes boundary points 
            net_size_x =size_x+2*boundary_width; // total matrix size includes boundary points 
            physical_x0=x0;
            physical_y0=y0;
            physical_x1=x1;
            physical_y1=y1;
            physical_dx=(x1-x0)/(size_x); // physical interval 
            physical_dy=(y1-y0)/(size_y);
            shifts=boundary_width;
            correlator=correlat;

             
            std::vector<std::vector<double> >  Grid_map_x(size_x , std::vector<double> (size_y,0) );
            std::vector<std::vector<double> >  Grid_map_y(size_x , std::vector<double> (size_y,0) );

            projected_physical_grid_x=Grid_map_x;// Those are place holders that will be replaced in the later stage. // point the public into private to call
            projected_physical_grid_y=Grid_map_y;// Those are place holders that will be replaced in the later stage. // point the public into private to call
           
           // internal boundaries 
          // if (internal_boundaries!=0){
            
          // };
          // DBC how far they enter into the boundary, is defined by the ( b_left,int b_right,int b_top,int b_bot)
            std::vector<std::vector<double> >  target_field_with_bc(net_size_x , std::vector<double> (net_size_y,0) );
          // define two point difference function 
            std::vector<std::vector<double> >  difffunction_x(net_size_x , std::vector<double> (2*correlator,0) );
            std::vector<std::vector<double> >  diff_function_y(net_size_y , std::vector<double> (2*correlator,0) );
            difference_function_x=difffunction_x;
            difference_function_y=diff_function_y;
            field_with_bc=target_field_with_bc;
      

    
    };

//................Methods ...................................................................

        //------------------Grid addition------------------------------------------- 
        // bind a given matrix with the boundary with given parameters
    std::vector<std::vector<double> > boundary_binding(std::vector<std::vector<double> > target_field,double fixed_boundary_value=1.3);

        //------------------Grid reduction------------------------------------------- 
    // remove the boundary 
    std::vector<std::vector<double> > boundary_unbinding(std::vector<std::vector<double> > target_field);


    //------------------value projection on grid-------------------------------------------- 
              
    std::vector<std::vector<double> > get_grid(std::vector<std::vector<double> > &current_grid,bool display=false);


};
#endif

