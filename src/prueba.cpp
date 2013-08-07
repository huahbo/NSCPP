
#include <iostream>
#include "pressureCorrection.h"
#include "Mesh.h"
#include "Interpol.hpp"
#include "Momentum.hpp"
#include "RungeKutta.h"
#include "MathInterpolation.h"
#include "Cells.h"


using namespace blitz;

int main(){
    MeshTool::MeshBlock Grid(10.0,10.0,400,7);
    Grid.cellDifference();
    cout << Grid.Dx << endl;    
    cout << Grid.Dy << endl;    
    cout << Grid.Nx << endl;    
    cout << Grid.Ny << endl;    
    cout << Grid.cNx << endl;    
    cout << Grid.cNy<< endl;    

    Array<real, 1> deltas(2);
    
    deltas(0) = Grid.Dx; 
    deltas(1) = Grid.Dy; 

    Array<real,2> U(Grid.Nx , Grid.cNy), 
                  V(Grid.cNx , Grid.Ny),
                  Pressure(Grid.cNx, Grid.cNy),
                  vGhostB(Grid.cNx, 6), 
                  uGhostB(Grid.cNy, 6); 

    Array<real,2> uBoundD(2,Grid.Nx),
                  vBoundD(2,Grid.Ny);


    U(Range::all(),Range::all()) = 1.0;

    for(int j = 0; j < Grid.cNy; j++){
        U(0,j) = 2.0;
    }

    MathInterpolation Us(U);

    StaggCellX Ucell(U, uBoundD, Grid);

    Ucell.Solve();

    std::cout << U << std::endl;
    std::cout << Us.I() << std::endl;
/*


    Array<real,2> BoundaryX(Grid.Nx,Grid.Ny),
                  BoundaryY(Grid.Nx,Grid.Ny);


   // Interpolation<real,1> UIJ_to_M(BoundaryX,uBoundD);
   // Interpolation<real,1> VIJ_to_M(BoundaryY,uBoundD);

    
    //Interpolation<real,1> UIJ_to_M(BoundaryX);
    //Interpolation<real,2> VIJ_to_M(BoundaryY);

    Interpolation<real,1> Ucell(U,uBoundD);
    Interpolation<real,2> Vcell(V,vBoundD);

    Momentum<real,1> MomentumX(U,V,Pressure,Ucell,Vcell,deltas);
    Momentum<real,2> MomentumY(U,V,Pressure,Ucell,Vcell,deltas);

    RungeKutta TimeInteger;


//     Ucell.IJ();
    


    
    for(int  i = 0 ; i < Grid.Nx ; i++){
        for(int  j = 0 ; j <  Grid.cNy; j++){
            U(i, j) = 0.0;
        }
   }

    for(int  i = 0 ; i < Grid.cNx ; i++){
        for(int  j = 0 ; j <  Grid.Ny; j++){
            V(i, j) = 0.0;
        }
   }
*/
    
    /**
    Ghost Boundaries (Explicit Method for Pressure Correction)
    **/
/*
    for(int  i = 0 ; i < Grid.cNy ; i++){
        for(int  j = 0 ; j < 6 ; j++){
            uGhostB(i, j) = 1.0;
        }
    }

    for(int  i = 0 ; i < Grid.cNx ; i++){
        for(int  j = 0 ; j < 6 ; j++){
           vGhostB(i, j) = 1.0;
        }
    }
*/

    /**
    Dirichlet Boundaries for Momentum Axis X and Y (Implicit Method) 
    **/
/*
    for(int  i = 0 ; i < Grid.Nx ; i++){
       uBoundD(0,i) = 0.0;
       uBoundD(1,i) = 1.0;
    }
    for(int  j = 0 ; j < Grid.Ny ; j++){
       vBoundD(0,j) = 0.0;
       vBoundD(1,j) = 0.0;
    }


    






    real dt = 0.01;

    PressureCorrect PressureEq(U, V, Pressure, uGhostB, vGhostB, Grid);
    PressureEq.Make();


    MomentumX.Solve();
    MomentumY.Solve();

    //U = TimeInteger.Solve(U,MomentumX.Divergence,0.5*dt);    
    //V = TimeInteger.Solve(V,MomentumY.Divergence,0.5*dt);    


    PressureEq.Solve(0.5*dt);
   // U = U - PressureEq.CorrectU(0.5*dt);
   // V = V - PressureEq.CorrectV(0.5*dt);
    

    std::cout << MomentumX.Divergence << std::endl;
  //  std::cout << BoundaryX << std::endl;
*/



return 0;
}


