#include "../include/pressureCorrection.h"
#include "iostream"

void PressureCorrect::Div()
{
    real inverseConst = 1.0 / 24.0;

  
  
    /**
     Divergence direction Axis X
    **/
     
     for(int i = 1; i < Grid.Nx - 2; i++){ //Without X - Boundaries
         for(int j = 0; j < Grid.cNy; j++){
  
             D(i,j) = (-1.0*uS(i+2,j) + 27.0*uS(i+1,j) - 27.0*uS(i,j) + uS(i-1,j))*
                      inverseConst*Grid.inverseDx; 
  
         }
     }
  
   
     /**
      Divergence direction Axis Y + Divergence direction X
     **/
   
      for(int i = 0; i < Grid.cNx; i++){
          for(int j = 1; j < Grid.Ny -2; j++){ //Without  Y - Boundaries
              
              D(i,j) += (-1.0*vS(i,j+2) + 27.0*vS(i,j+1) - 27.0*vS(i,j) + vS(i,j-1))*
                      inverseConst*Grid.inverseDy;
          }
      }
     
   
    /**
     Boundary Axis X
    **/ 
  
  
     for(int j = 0; j < Grid.cNy; j++){
  
         D(0,j) += (-1.0*uS(2,j) + 27.0*uS(1,j) - 27.0*uS(0,j) + uBoundary(0,j))
                 *inverseConst*Grid.inverseDx;
  
         D(Grid.cNx-1,j) += (-1.0*uBoundary(3,j) + 27.0*uS(Grid.Nx-1,j) - 27.0*uS(Grid.Nx-2,j) +  uS(Grid.Nx-3,j))
                     *inverseConst*Grid.inverseDx; 
  
     }
  
    /**
     Boundary Axis Y
    **/ 
  
     for(int i = 0; i < Grid.cNx; i++){
         
         D(i,0) += (1.0*vS(i,2) + 27.0*vS(i,1) - 27.0*vS(i,0) + vBoundary(i,0))
                 *inverseConst*Grid.inverseDy;
  
         D(i,Grid.cNy-1) += (-1.0*vBoundary(i,3) + 27.0*vS(i,Grid.Ny-1) - 27.0*vS(i,Grid.Ny-2)+ vS(i,Grid.Ny-3))
                      *inverseConst*Grid.inverseDy;
  
     }

     std::cout << uS;
     std::cout << vS;
     std::cout << D;

}
