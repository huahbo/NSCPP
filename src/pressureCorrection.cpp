#include "pressureCorrection.h"
#include "MathCore.h"
#include <iostream>

void PressureCorrect::Div(real stage)
{
    real inverseConst = 1.0 / 24.0;

  
  
    /**
     Divergence direction Axis X
    **/
     
     for(int i = 1; i < Grid.Nx - 2; i++){ //Without X - Boundaries
         for(int j = 0; j < Grid.cNy; j++){
  
             D(i,j) = (-1.0*uS(i+2,j) + 27.0*uS(i+1,j) - 27.0*uS(i,j) + uS(i-1,j))*
                      inverseConst*Grid.inverseDx; 


             int offset = i + j*Grid.cNx;
             RHS[offset] = D(i,j) / stage;  
  
         }
     }
  
   
     /**
      Divergence direction Axis Y + Divergence direction X
     **/
   
      for(int i = 0; i < Grid.cNx; i++){
          for(int j = 1; j < Grid.Ny -2; j++){ //Without  Y - Boundaries
              
              D(i,j) += (-1.0*vS(i,j+2) + 27.0*vS(i,j+1) - 27.0*vS(i,j) + vS(i,j-1))*
                      inverseConst*Grid.inverseDy;

             int offset = i + j*Grid.cNx;
             RHS[offset] = D(i,j) / stage;   
          }
      }
     
   
    /**
     Boundary Axis X
    **/ 
  
  
     for(int j = 0; j < Grid.cNy; j++){
  
         D(0,j) += (-1.0*uS(2,j) + 27.0*uS(1,j) - 27.0*uS(0,j) + uBoundary(0,j))
                 *inverseConst*Grid.inverseDx;

         int offset = j*Grid.cNx;
         RHS[offset]  = D(0,j) / stage;
  
         D(Grid.cNx-1,j) += (-1.0*uBoundary(3,j) + 27.0*uS(Grid.Nx-1,j) - 27.0*uS(Grid.Nx-2,j) +  uS(Grid.Nx-3,j))
                     *inverseConst*Grid.inverseDx; 

         offset = Grid.cNx-1 + j*Grid.cNx;
         RHS[offset] = D(Grid.cNx-1,j) / stage;

  
     }
  
    /**
     Boundary Axis Y
    **/ 
  
     for(int i = 0; i < Grid.cNx; i++){
         
         D(i,0) +=  (-1.0*vS(i,2) + 27.0*vS(i,1) - 27.0*vS(i,0) + vBoundary(0,i))
                 *inverseConst*Grid.inverseDy;

          int offset = i;
          RHS[offset] = D(i,0) / stage;

  
         D(i,Grid.cNy-1) +=(-1.0*vBoundary(3,i) + 27.0*vS(i,Grid.Ny-1) - 27.0*vS(i,Grid.Ny-2)+ vS(i,Grid.Ny-3))
                      *inverseConst*Grid.inverseDy;

          offset = i + (Grid.cNy-1)*Grid.cNx;
          RHS[offset] = D(i,Grid.cNy-1) / stage;
  
     }

    
}






void PressureCorrect::Make()
{

    FillPoisson();


    
    for(int i = 0; i < elements; i++){
        id[i] = row2[i] + Grid.cNx*Grid.cNy*col2[i];
    }
         
    insertionSort(id, col2, row2, val2, elements);

    add_ptr(ptr_col, Grid.cNx*Grid.cNy , col2, elements);  

}

void PressureCorrect::Solve(real stage)
{
    Div(stage);
    for(int i = 0; i < Grid.cNx*Grid.cNy; i++){
        x_estimate[i] = (double)1;
    }
    
    pmgmres_ilu_cr(Grid.cNx*Grid.cNy, elements, ptr_col,row2,val2,x_estimate,
                   RHS, 10, 2*Grid.cNx, tol, tol);

/*
    for(int i = 0; i < Grid.cNx*Grid.cNy; i++){
       //std::cout << x_estimate[i] << std::endl;
       std::cout << RHS[i] << std::endl;
    }
*/
    for(int i = 0; i < Grid.cNx ; i++){ //Without X - Boundaries
        for(int j = 0; j < Grid.cNy; j++){
            int offset = i + j*Grid.cNx;
            pressure(i,j) += x_estimate[offset];
        }
    }
    
}
 

void PressureCorrect::FillPoisson() 
{

    real Banded = (-1460.0/576.0)*(
                                   Grid.inverseDx*Grid.inverseDx 
                                 + Grid.inverseDy*Grid.inverseDy
                                  );

    real C1 = (783.0/576.0)*Grid.inverseDx*Grid.inverseDx;
    real C2 = (-54.0/576.0)*Grid.inverseDx*Grid.inverseDx;
    real C3 = (1.0/576.0)*Grid.inverseDx*Grid.inverseDx;

    
    real C4 = (783.0/576.0)*Grid.inverseDy*Grid.inverseDy;
    real C5 = (-54.0/576.0)*Grid.inverseDy*Grid.inverseDy;
    real C6 = (1.0/576.0)*Grid.inverseDy*Grid.inverseDy;

    int memplace = 0;
    int Nx = Grid.cNx;
    int Ny = Grid.cNy;

    MATRIX_double PP(
                     Grid.cNx*Grid.cNy,
                     Grid.cNx*Grid.cNy,
                     0.0
                    );


    MATRIX_double TT(
                     Grid.cNx*Grid.cNy,
                     Grid.cNx*Grid.cNy,
                     0.0
                    );



    for(int i = 0; i < Nx; i++){
       for(int j =  0; j < Ny; j++){
          int k = i + j*Nx;
          
          PP(k,k) = Banded; //main banded

          val2[memplace] = Banded;
          row2[memplace] = k;
          col2[memplace] = k;
          memplace++;

       }
    }   



    for(int i = 3; i < Nx - 3; i++){
       for(int j =  0; j < Ny; j++){
          int k = i + j*Nx;
          
           
          PP(k,k+1) = C1; //RS-C1

          val2[memplace] = C1;
          row2[memplace] = k;
          col2[memplace] = k+1;
          memplace++;

          PP(k,k-1) = C1; //LS-C1

          val2[memplace] = C1;
          row2[memplace] = k;
          col2[memplace] = k-1;
          memplace++;

          PP(k,k+2) = C2; //RS-C2
	
          val2[memplace] = C2;
          row2[memplace] = k;
          col2[memplace] = k+2;
          memplace++;

          PP(k,k-2) = C2; //LS-C2

          val2[memplace] = C2;
          row2[memplace] = k;
          col2[memplace] = k-2;
          memplace++;


          PP(k,k+3) = C3; //RS-C3

          val2[memplace] = C3;
          row2[memplace] = k;
          col2[memplace] = k+3;
          memplace++;

          PP(k,k-3) = C3; //LS-C3

          val2[memplace] = C3;
          row2[memplace] = k;
          col2[memplace] = k-3;
          memplace++;

       }
    }

    /**

    Boundary Values ¡¡¡ Von Neumman !!!

    **/


    // AXIS - X (Interior Matriz)
    for(int j =  0; j < Ny; j++){

       ///////////////////////////////////////////////////////////
       int k = j*Nx;
       PP(k,k+1) = C1*2.0;

       val2[memplace] = C1*2.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;


       PP(k,k+2) = C2*2.0;

        val2[memplace] = C2*2.0;
        row2[memplace] = k;
        col2[memplace] = k+2;
        memplace++;


       PP(k,k+3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;

       //////////////////////////////////////////////////////////// 
       k = 1 + j*Nx;
       PP(k,k+1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = C2*2.0;

       val2[memplace] = C2*2.0;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k+3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;
       
       ////////////////////////////////////////////////////////////
       k = 2 + j*Nx;
       PP(k,k+1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = C2;

       val2[memplace] = C2;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k-2) = C2;

       val2[memplace] = C2;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;


       PP(k,k+3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;
    
       ////////////////////////////////////////////////////////////
       k = (Nx-1) + j*Nx;
       PP(k,k-1) = C1*2.0;

       val2[memplace] = C1*2.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;


       PP(k,k-2) = C2*2.0;

       val2[memplace] = C2*2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k-3;
       memplace++;
       
       //////////////////////////////////////////////////////////// 
       k = (Nx -2)  + j*Nx;
       PP(k,k+1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k-2) = C2*2.0;

       val2[memplace] = C2*2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k-3;
       memplace++;
       
       /////////////////////////////////////////////////////////////
       k = (Nx - 3) + j*Nx;
       PP(k,k+1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = C1;

       val2[memplace] = C1;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = C2;

       val2[memplace] = C2;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k-2) = C2;

       val2[memplace] = C2;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = C3*2.0;

       val2[memplace] = C3*2.0;
       row2[memplace] = k;
       col2[memplace] = k-3;
       memplace++;

    }



    // Y Axis Values
    for(int i = 0; i < Nx; i++)
    {
       for(int j =  3; j < Ny - 3; j++)
       {

          int k = i + j*Nx;
          PP(k,k+1*Nx) = C4; 

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k-1*Nx) = C4; //US-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

          PP(k,k+2*Nx) = C5; //DS-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;

          PP(k,k-2*Nx) = C5; //US-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k+3*Nx) = C6; //DS-C3

          val2[memplace] = C6;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;

          PP(k,k-3*Nx) = C6; //US-C3

          val2[memplace] = C6;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;

       }
   }         



    //  Y Boundary Values 
    for(int i = 0; i < Nx; i++){

      {  
          int k = i;
          PP(k,k+1*Nx) = C4*2.0; //DS-C1

          val2[memplace] = C4*2.0;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = C5*2.0; //DS-C2

          val2[memplace] = C5*2.0;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;


          PP(k,k+3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;
      }



      {  
          int k = i + Nx;
          PP(k,k+1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = C5*2.0; //DS-C2

          val2[memplace] = C5*2.0;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;


          PP(k,k+3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;
      }


      {  
          int k = i + 2*Nx;
          PP(k,k+1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = C5; //DS-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;


          PP(k,k+3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;
      }
 
      /////////////////////////////////////////////////////////////

      {
          int k = i + 1*Nx;
          PP(k,k-1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

      }  


      {
          int k = i + 2*Nx;
          PP(k,k-1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

          PP(k,k-2*Nx) = C5; //DS-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;
      }  

      
   }         


    //  Y Boundary Values 
    for(int i = 0; i < Nx; i++){

        {

          int k = i + (Ny-3)*Nx;
          PP(k,k-1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;


          PP(k,k-2*Nx) = C5; //DS-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k-3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;
      }


        {

          int k = i + (Ny-2)*Nx;
          PP(k,k-1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;


          PP(k,k-2*Nx) = C5*2.0; //DS-C2

          val2[memplace] = C5*2.0;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k-3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;
      }


        {

          int k = i + (Ny-1)*Nx;
          PP(k,k-1*Nx) = C4*2.0; //DS-C1

          val2[memplace] = C4*2.0;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;


          PP(k,k-2*Nx) = C5*2.0; //DS-C2

          val2[memplace] = C5*2.0;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k-3*Nx) = C6*2.0; //DS-C3

          val2[memplace] = C6*2.0;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;
      }


      {
          int k = i + (Ny-2)*Nx;
          PP(k,k+1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;
      }  


      {
          int k = i + (Ny-3)*Nx;
          PP(k,k+1*Nx) = C4; //DS-C1

          val2[memplace] = C4;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = C5; //DS-C2

          val2[memplace] = C5;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;
      }  
   }         


   for(int i = 0; i < elements; i++){
      TT(row2[i] , col2[i]) = val2[i]; 
   }

//    std::cout << PP << std::endl << std::endl;
//    std::cout << TT << std::endl;
    std::cout << elements << std::endl;
    std::cout << memplace << std::endl;

}
