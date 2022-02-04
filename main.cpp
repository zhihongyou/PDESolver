#include <iostream> 
#include <vector>
#include <string>
#include "Field_v2.h"
#include "field_v2.cpp"
#include "field_access.cpp"
#include "stepping.cpp"

int main() {
    // Setting up the raw grid. 
    Fields mesh(5,5,20,-5,-5,5,5,1,1);  //(a,b,c): a,b means amount of data points span. c means the data point density at the point.

          // initial relevent parameters 
          std::vector<std::vector<double> >  current_laplacian,current_field;
          std::vector<std::vector<double> >  projected_field_with_boundary,projected_field;
          std::vector<std::vector<double> >  next_time_projected_field;

          //define time steps with initial parameters 
          double dt,T;
          int n, steps;
          string field_path, grad_field_path ;
          n=0; // initialize current stepping.
          steps=5; // total stepping asked. 
          dt=0.3; // time interval for each step. 
          T=0; // initialize total time. 
          // Project the field onto the raw mesh. 
          current_field=mesh.gaussian_field(); // One must run this before calculating the laplacian 
         // mesh.get_grid(current_field);
          // modify the field with added boundary condition. It will return the updated field with boundary condition applied 
          projected_field_with_boundary=mesh.boundary_binding(current_field,0); // initial boundary condition
          //mesh.get_grid(projected_field_with_boundary);

          field_path="output/field.csv";
          grad_field_path="output/field_grad.csv";
          mesh.initial_file_loc(field_path);
          mesh.initial_file_loc(grad_field_path);

          //Export the initial field to the file.
          mesh.out_field();
          while (n<steps) // condition on the step
          { 
            // Calculate the field gradient, and return the gradient term. 
            //One can supply any field since it does not have implicit usages for the projected_field defined above. It will calculate the laplacian for the supplied field. 
            current_laplacian=mesh.calculate_field_laplacian(projected_field_with_boundary);
            //mesh.get_original_field(true);
            // stepping toward next time step, and return the projected field. 
            //Similar to gridient term, there is no implcit uses for the projected field. It will step the supplied field. 
            next_time_projected_field=mesh.newton_stepping(current_laplacian,dt,1);
            // apply the boundary condition again. 
            projected_field=mesh.boundary_binding(next_time_projected_field,0);
            // update the current field into the object. It will updated stored field with supplied field. 
            mesh.update_field_with_bc(projected_field);


            // export current fields and field_grad for this time step. 
            mesh.out_field_grad();
            mesh.out_field();

            // stepping 
            n++;
            // Update the total time. 
            T=T+dt; 
          };
          // output the mesh grid for x and y. 
          mesh.out_meshgird();
    return 0;
};
