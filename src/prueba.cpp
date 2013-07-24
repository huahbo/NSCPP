
#include <iostream>
#include "pressureCorrection.h"
#include "Mesh.h"


using namespace blitz;

int main(){

    MeshTool::MeshBlock Grid(10.0,10.0,10,10);
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


     PressureCorrect PressureEq(U, V, Pressure, uBoundary, vBoundary, Grid);

     PressureEq.Div();
     PressureEq.Update();

return 0;
}

