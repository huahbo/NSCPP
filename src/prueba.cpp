
#include <iostream>
#include "pressureCorrection.h"
#include "Mesh.h"


using namespace blitz;

int main(){

    MeshTool::MeshBlock Grid(10.0,10.0,100,100);
    Grid.cellDifference();
    cout << Grid.Dx << endl;    
    cout << Grid.Dy << endl;    
    cout << Grid.Nx << endl;    
    cout << Grid.Ny << endl;    
    cout << Grid.cNx << endl;    
    cout << Grid.cNy<< endl;    
  

    Array<real, 2> U(Grid.Nx , Grid.cNy), 
                   V(Grid.cNx , Grid.Ny),
                   Pressure(Grid.cNx, Grid.cNy),
                   uBoundary(Grid.Nx, 6), 
                   vBoundary(6, Grid.Ny); 


    for(int  i = 0 ; i < Grid.Nx ; i++){
        for(int  j = 0 ; j <  Grid.cNy; j++){
            U(i, j) = 1.0;
        }
   }

    
    for(int  i = 0 ; i < Grid.cNx ; i++){
        for(int  j = 0 ; j <  Grid.Ny; j++){
            V(i, j) = 1.0;
        }
   }

    
    for(int  i = 0 ; i < Grid.Nx ; i++){
        for(int  j = 0 ; j < 6 ; j++){
            uBoundary(i, j) = -1.0;
        }
    }

    for(int  i = 0 ; i < Grid.cNy ; i++){
        for(int  j = 0 ; j < 6 ; j++){
            uBoundary(j, i) = -1.0;
        }
    }


    PressureCorrect PressureEq(U, V, Pressure, uBoundary, vBoundary, Grid);
    PressureEq.Make();
    PressureEq.Solve(0.05);

return 0;
}

