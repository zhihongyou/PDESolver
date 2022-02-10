#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
// #include "src/evolver/evolverclass.cpp"
// #include "src/utility/finiteDifference.h"
#include "src/utility/finiteDifferenceCentralO2Isotropic.h"

using namespace std;


// ======================================================================
int main() {
    // FiniteDifferenceScheme* scheme;
    // scheme=&central_difference_O4_I;
    int N=10000;
    double* a;    
    a=new double[N];
    for (int i=0; i<N; i++) {            
        a[i]=1;
    };
    std::string method="O4I";
    
    clock_t t0=clock();
    for (int j=0; j<100000; j++) {
        for (int i=0; i<N; i++) {
            // if (method=="O4I") {
            // a[i]=FDM<&central_O4I>::add(a,i);
            // } else if (method == "O4") {
            a[i]=FDM<&central_O2I>::add(a,i);
            // } else if (method == "O2") {
            //     a[i]=FDM<&central_difference_O4>::add(a,i);
            // } else if (method == "O2") {
            //     a[i]=FDM<&central_difference_O4>::add(a,i);
            // } else if (method == "O2") {
            //     a[i]=FDM<&central_difference_O4>::add(a,i);
            // } else if (method == "O2") {
            //     a[i]=FDM<&central_difference_O4>::add(a,i);
            // } else if (method == "O2") {
            //     a[i]=FDM<&central_difference_O4>::add(a,i);
            // } else {
            //     a[i]=add(a,i);
            // };
            // a[i]=add(a,i);
        };
    };
    clock_t ts = clock();
    double tused=double(ts-t0)/CLOCKS_PER_SEC;
    cout <<a[1]<<endl;
    cout <<tused<<endl;

    // Generating a new system.
    // System mySys;
    // Generating a new mesh.
    // Mesh mesh(2);
    // Add mesh to system.
    // mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    // Field phi(&mesh, "phi", 0, 0, "periodic", "Gaussian", "on");
    // // Set field equations.
    // phi.rhs_terms={{"*",{{"laplace",&phi}},{1},"explicit"}};
    // // Add fields to the system.
    // mySys.field_ptrs.push_back(&phi);
    // // Print system information.
    // mySys.printSysInfo();
    
    // cout<< (*mySys.field_ptrs[0]).name <<endl;
    // cout<< (*mySys.field_ptrs[0]).f_cpu[0][200] <<endl;
    
    // phi.getLalpace(0);
    // phi.expo_conf_any(
    
    // Creating an evolver:
    // Evolver evolver(&mySys,0,100,0.001,1,"EulerForward");
    // Run simulations.
    // evolver.run();
    
    
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
    // Field mub(&mesh, "mub", 0, 1, "periodic", "sin", "off");
    // Field phia(&mesh, "phia", 0, 0, "periodic", "Gaussian", "on");    
    // Field phib(&mesh, "phib", 0, 0, "periodic", "Gaussian", "on");
    // // Write field equations
    // double M=4;
    // double gamma=5;
    // double kab=-1;
    // double kba=1;
    // mua.rhs_terms={
    //     {"*",{{"laplace",&phia}},{-gamma},"explicit"},
    //               {"*",{{"1",&phia}},{-1},"explicit"},
    //               {"*",{{"1",&phia},{"1",&phia},{"1",&phia}},{1},"explicit"}};
    // phia.rhs_terms={{"*",{{"laplace",&mua}},{M},"explicit"},
    //                 {"*",{{"laplace",&phib}},{M*kab},"explicit"}};
    // mub.rhs_terms={{"*",{{"1",&phib}},{1},"explicit"}};
    // phib.rhs_terms={{"*",{{"laplace",&mub}},{M},"explicit"},
    //                 {"*",{{"laplace",&phia}},{M*kba},"explicit"}};
    // // Add fields to system.
    // mySys.field_ptrs.push_back(&mua);
    // mySys.field_ptrs.push_back(&mub);
    // mySys.field_ptrs.push_back(&phia);
    // mySys.field_ptrs.push_back(&phib);
    // Evolver evolver(&mySys,0,1000,0.002,1,"EulerForward");
