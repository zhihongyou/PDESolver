#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "src/evolver/evolverclass.cpp"


using namespace std;


// ======================================================================
int main() {    

    // Generating a new system.
    System mySys;
    // Generating a new mesh.
    Mesh mesh(2);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field phi(&mesh, "phi", 0, 0, "periodic", "Gaussian", "on");
    // Set field equations.
    phi.rhs_terms={{"*",{{"laplace",&phi}},{1},"explicit"}};
    // Add fields to the system.
    mySys.field_ptrs.push_back(&phi);
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:
    Evolver evolver(&mySys,0,100,0.001,1,"EulerForward");
    // Run simulations.
    evolver.run();
    
    
    return 0;
};


// Field phi(&mesh, "phi", 0, 0, "periodic", "Gaussian", "on");
// phi.rhs_terms={{"*",{{"laplace",&phi}},{1},"explicit"}};

    // double M=1;
    // double gamma=1;
    // Field mu(&mesh, "mu", 0, 1, "periodic", "sin", "on");    
    // Field rho(&mesh, "rho", 0, 0, "periodic", "sin", "on");
    // mu.rhs_terms={{"*",{{"laplace",&rho}},{-gamma},"explicit"},
    //               {"*",{{"1",&rho}},{-1},"explicit"},
    //               {"*",{{"1",&rho},{"1",&rho},{"1",&rho}},{1},"explicit"}};    
    // rho.rhs_terms={{"*",{{"laplace",&mu}},{M},"explicit"}};


// M=4;
// gamma=5;
// double kab=-0.3;
// double kba=0.3;
// mua.rhs_terms={{"*",{{"laplace",&phia}},{-gamma},"explicit"},
//                {"*",{{"1",&phia}},{-1},"explicit"},
//                {"*",{{"1",&phia},{"1",&phia},{"1",&phia}},{1},"explicit"}};
// phia.rhs_terms={{"*",{{"laplace",&mua}},{M},"explicit"},
//                 {"*",{{"laplace",&phib}},{M*kab},"explicit"}};
// mub.rhs_terms={{"*",{{"laplace",&phib}},{-gamma},"explicit"},
//                {"*",{{"1",&phib}},{-1},"explicit"},
//                {"*",{{"1",&phib},{"1",&phib},{"1",&phib}},{1},"explicit"}};
// phib.rhs_terms={{"*",{{"laplace",&mub}},{M},"explicit"},
//                 {"*",{{"laplace",&phia}},{M*kba},"explicit"}};


    // Field mua(&mesh, "mua", 0, 1, "periodic", "sin", "off");    
//     Field mub(&mesh, "mub", 0, 1, "periodic", "sin", "off");
//     Field phia(&mesh, "phia", 0, 0, "periodic", "Gaussian", "on");    
//     Field phib(&mesh, "phib", 0, 0, "periodic", "Gaussian", "on");
//     // Write field equations
//     double M=4;
//     double gamma=5;
//     double kab=-1;
//     double kba=1;
//     mua.rhs_terms={
//         {"*",{{"laplace",&phia}},{-gamma},"explicit"},
//                   {"*",{{"1",&phia}},{-1},"explicit"},
//                   {"*",{{"1",&phia},{"1",&phia},{"1",&phia}},{1},"explicit"}};
//     phia.rhs_terms={{"*",{{"laplace",&mua}},{M},"explicit"},
//                     {"*",{{"laplace",&phib}},{M*kab},"explicit"}};
//     mub.rhs_terms={{"*",{{"1",&phib}},{1},"explicit"}};
//     phib.rhs_terms={{"*",{{"laplace",&mub}},{M},"explicit"},
//                     {"*",{{"laplace",&phia}},{M*kba},"explicit"}};
//     // Add fields to system.
//     mySys.field_ptrs.push_back(&mua);
//     mySys.field_ptrs.push_back(&mub);
//     mySys.field_ptrs.push_back(&phia);
//     mySys.field_ptrs.push_back(&phib);
// Evolver evolver(&mySys,0,1000,0.002,1,"EulerForward");
