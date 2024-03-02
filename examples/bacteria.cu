#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <iostream> 
#include <vector>
#include <list>
#include <cstdio>
#include <string>
#include "cuda.h"
#include "curand.h"
#include "curand_kernel.h"
#include "cuda_runtime_api.h"
#include <cmath>
#include <ctime>
#include <cufft.h>
#include <cufftXt.h>
#include "userDefinedFunction.h"
#include "../src/evolver/evolverclass.cu"


using namespace std;


// ======================================================================
int main() {
    
    // Simulation parameters
    string direExpo="data/";
    string device="gpu";
    string FDMScheme="CentralDifferenceO4Iso2D";
    string timeScheme="RK4";
    double dt=0.001;
    double T=100;
    double dtExpo=5;
    int    NGrid=128;
    double L=2*NGrid;
    

    // Generating a new system.
    System mySys(direExpo);
    // Generating a new mesh.
    Mesh mesh(2);
    mesh.setGridNumber(NGrid,NGrid,1);
    mesh.setBoxSize(L,L,1);
    // Add mesh to system.
    mySys.mesh_ptr=&mesh;
    
    // Creating fields.
    Field S2(&mesh, "S2",1);
    Field Hxx(&mesh,"Hxx",1);
    Field Hxy(&mesh,"Hxy",1);
    Field q2oq1(&mesh,"q2oq1",1);
    Field QxQy2(&mesh,"QXQY2",1);
    Field theta(&mesh,"theta",1);    
    Field n1(&mesh,"n1",1);
    Field n2(&mesh,"n2",1);
    Field Qxx(&mesh, "Qxx",0,"sin");
    Field Qxy(&mesh, "Qxy",0,"sin");
    Field cp(&mesh,"cp",0);
    Field cm(&mesh,"cm",0);
    Field Exx(&mesh,"Exx",1);
    Field Exy(&mesh,"Exy",1);
    Field A(&mesh,"A",1);
    Field Omega(&mesh,"Omega",1);
    Field trQW(&mesh,"trQW",1);
    Field sigmaA(&mesh,"sigmaA",1);
    Field sigmaS1(&mesh,"sigmaS1",1);
    Field sigmaS2(&mesh,"sigmaS2",1);
    Field trQH(&mesh,"trQH",1);
      
    Field phi0(&mesh,"phi0",1,"zero");
   
    
  
//LaplaceNFEqField omega1(&mesh,"omega1",1);


    // Field Psi(&mesh, "Psi",0,"sin");    
    Field dxvy(&mesh,"dxvy",1);
    Field dyvx(&mesh,"dyvx",1);
    IncompFlow incompFlow(&mesh, "incompFlow",1);
    
    
    cp.initFieldConst(0.11);
    cm.initFieldConst(0.1);
    // Set equations
    double K = 15;
    //double Er = 3.75;
    double Gamma = 1;
    double eta = 0.5;
    double xian = 0;
    double h = 20;
    double xi = 0.9;  
    double a = 0.4;
    double b = 0.8;
    //double phi = 0;
    double l = 5;
    double V0 = 15;
    double tau  = 50;
    double Dc = 200;
    
    double AA =  187.5; 
    // 187.5;
    double zta = 12*eta/h/2;
incompFlow.setOmegaEq(2,zta,eta);
    dxvy.setRhsTerms({
        {1,{{"d1x",&incompFlow.vy}}}
    });
    dyvx.setRhsTerms({
        {1,{{"d1y",&incompFlow.vx}}}
    });
    
    S2.setRhsTerms({
        {4,{{&Qxx},{&Qxx}}},
        {4,{{&Qxy},{&Qxy}}}
    });

    Hxx.setRhsTerms({
        {a,{{&Qxx}}},{-b/2,{{&S2},{&Qxx}}},
        {K,{{"laplace",&Qxx}}}
      });
    Hxy.setRhsTerms({
        {a,{{&Qxy}}},{-b/2,{{&S2},{&Qxy}}},
        {K,{{"laplace",&Qxy}}}
      });
    cp.setRhsTerms({
        {-V0,{{"d1x",&n1},{&cp}}},
        {-V0,{{"d1y",&n2},{&cp}}},
        {-V0,{{"d1x",&cp},{&n1}}},
        {-V0,{{"d1y",&cp},{&n2}}},
        {-1,{{"d1x",&cp},{&incompFlow.vx}}},
        {-1,{{"d1y",&cp},{&incompFlow.vy}}},
        {Dc,{{"laplace",&cp}}},
        {-1/tau,{{&cp}}},
        {1/tau,{{&cm}}},
    });
    cm.setRhsTerms({
        {-V0,{{"d1x",&n1},{&cm}}},
        {-V0,{{"d1y",&n2},{&cm}}},
        {-V0,{{"d1x",&cm},{&n1}}},
        {-V0,{{"d1y",&cm},{&n2}}},
        {-1,{{"d1x",&cm},{&incompFlow.vx}}},
        {-1,{{"d1y",&cm},{&incompFlow.vy}}},
        {Dc,{{"laplace",&cm}}},
        {-1/tau,{{&cm}}},
        {1/tau,{{&cp}}},
    });

    Qxx.setRhsTerms({
        // {-1,{{"d1x",&Qxx},{&incompFlow.vx}}},
        // {-1,{{"d1y",&Qxx},{&incompFlow.vy}}},
        {Gamma,{{&Hxx}}}

        // {2*xi,{{&Qxx},{"d1x",&incompFlow.vx}}},
        // {xi,{{"d1x",&incompFlow.vx}}},

        // {xi,{{&Qxy},{"d1x",&incompFlow.vy}}},
        // {xi,{{&Qxy},{"d1y",&incompFlow.vx}}},
        // {-2*xi,{{&Qxx},{&trQW}}},
        // {-xi,{{&trQW}}},
    });
    Qxy.setRhsTerms({
        // {-1,{{"d1x",&Qxy},{&incompFlow.vx}}},
        // {-1,{{"d1y",&Qxy},{&incompFlow.vy}}},
        {Gamma,{{&Hxy}}}

        
        // {xi,{{&A}}},
        // {-2,{{&Qxx},{&Omega}}},
    
        // {-2*xi,{{&Qxy},{&trQW}}},
    });

    trQW.setRhsTerms({
        {1,{{&Qxy},{"d1x",&incompFlow.vy}}},
        {1,{{&Qxy},{"d1y",&incompFlow.vx}}},
        {2,{{&Qxx},{"d1x",&incompFlow.vx}}}
    });

    A.setRhsTerms({
        {0.5,{{&dyvx}}},
        {0.5,{{&dxvy}}},
        
    });
    Omega.setRhsTerms({
        {0.5,{{&dxvy}}},
        {-0.5,{{&dyvx}}},
    });
    incompFlow.omega.setRhsTerms({
        {1,{{"laplace",&sigmaA}}},
        {1,{{"d2x",&sigmaS2}}},
        {-1,{{"d2y",&sigmaS2}}},
        {-2,{{"d1x1y",&sigmaS1}}}


    });
    sigmaA.setRhsTerms({
        {2,{{&Qxy},{&Hxx}}},
        {-2,{{&Qxx},{&Hxy}}},
    });
    sigmaS1.setRhsTerms({
         {-K,{{"d1x",&Qxx},{"d1x",&Qxx}}},{K,{{"d1y",&Qxy},{"d1y",&Qxy}}},
        {K,{{"d1y",&Qxx},{"d1y",&Qxx}}},{-K,{{"d1x",&Qxy},{"d1x",&Qxy}}},
       
        {-2*xi,{{&Hxy},{&Qxy}}},
        {-2*xi,{{&Hxx},{&Qxx}}},
        {-xi,{{&Hxx}}},
        {2*xi,{{&Qxx},{&trQH}}},
        {xi,{{&trQH}}},
        {AA,{{&cp},{&Qxx}}},
        {AA,{{&cm},{&Qxx}}}


        
        // {-2*K,{{"d1x",&Qxx},{"d1x",&Qxx}}},
        // {-2*K,{{"d1x",&Qxy},{"d1x",&Qxy}}},

    });

    sigmaS2.setRhsTerms({
         {-2*K,{{"d1x",&Qxy},{"d1y",&Qxy}}},{-2*K,{{"d1x",&Qxx},{"d1y",&Qxx}}},
        {-xi,{{&Hxy}}},
        {2*xi,{{&Qxy},{&trQH}}},
        {AA,{{&cp},{&Qxy}}},
        {AA,{{&cm},{&Qxx}}},

    });

    trQH.setRhsTerms({
        {2,{{&Hxx},{&Qxx}}},
        {2,{{&Hxy},{&Qxy}}},
    });




    q2oq1.setRhsTerms({
        {1,{{"1/f",&Qxx},{&Qxy}}}
    });
    QxQy2.setRhsTerms({
       {3,{{"sign",&Qxx}}},{1,{{"sign",&Qxy}}}
    });
    theta.setRhsTerms({
        {0.5,{{"atan",&q2oq1}}},
	{0.5,{{"toatan2",&QxQy2}}}      
    });

    n1.setRhsTerms({
        {1,{{"cos",&theta}}}
    });
    n2.setRhsTerms({
        {1,{{"sin",&theta}}}
    });    
    
    // mySys.addField(&trQH);
    // mySys.addField(&sigmaS1);
    // mySys.addField(&sigmaS2);
    // mySys.addField(&sigmaA);
    // //mySys.addField(&omega1);
    // mySys.addIncompFlow(&incompFlow);
    // mySys.addField(&dxvy);
    // mySys.addField(&dyvx);
    // mySys.addField(&Omega);
    // mySys.addField(&A);
    // mySys.addField(&trQW);
    mySys.addField(&S2);
    mySys.addField(&Hxx);
    mySys.addField(&Hxy);    
    mySys.addField(&q2oq1);
    mySys.addField(&QxQy2);
    mySys.addField(&theta);
    mySys.addField(&n1);
    mySys.addField(&n2);
    mySys.addField(&Qxx);
    mySys.addField(&Qxy);
    // mySys.addField(&cp);
    // mySys.addField(&cm);

    
    // Print system information.
    mySys.printSysInfo();
    
    // Creating an evolver:    
    Evolver evolver(&mySys,0,T,dt,dtExpo,device,timeScheme,FDMScheme);

    // Running simulations
    evolver.run();
    
    
    // -----------------------------------------------------------
    // Testing
    // evolver.initEvolver();
    // evolver.getRHS(0);
    
    // -----------------------------------------------------------
    
    return 0;
};
