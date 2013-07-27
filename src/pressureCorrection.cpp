#include "pressureCorrection.h"
#include "MathCore.h"
#include <iostream>

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


             int offset = i + j*Grid.cNx;
             b(offset) = D(i,j);  
  
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
             b(offset) = D(i,j);   
          }
      }
     
   
    /**
     Boundary Axis X
    **/ 
  
  
     for(int j = 0; j < Grid.cNy; j++){
  
         D(0,j) += (-1.0*uS(2,j) + 27.0*uS(1,j) - 27.0*uS(0,j) + uBoundary(0,j))
                 *inverseConst*Grid.inverseDx;

         int offset = j*Grid.cNx;
         b(offset) = D(0,j);
  
         D(Grid.cNx-1,j) += (-1.0*uBoundary(3,j) + 27.0*uS(Grid.Nx-1,j) - 27.0*uS(Grid.Nx-2,j) +  uS(Grid.Nx-3,j))
                     *inverseConst*Grid.inverseDx; 

         offset = Grid.cNx-1 + j*Grid.cNx;
         b(offset) = D(Grid.cNx-1,j);

  
     }
  
    /**
     Boundary Axis Y
    **/ 
  
     for(int i = 0; i < Grid.cNx; i++){
         
         D(i,0) += (1.0*vS(i,2) + 27.0*vS(i,1) - 27.0*vS(i,0) + vBoundary(i,0))
                 *inverseConst*Grid.inverseDy;

          int offset = i;
          b(offset) = D(i,0);

  
         D(i,Grid.cNy-1) +=(-1.0*vBoundary(i,3) + 27.0*vS(i,Grid.Ny-1) - 27.0*vS(i,Grid.Ny-2)+ vS(i,Grid.Ny-3))
                      *inverseConst*Grid.inverseDy;

          offset = i + (Grid.cNy-1)*Grid.cNx;
          b(offset) = D(i,Grid.cNy-1);
  
     }

    

  //   std::cout << uS;
    // std::cout << vS;
     //  std::cout << D;

}






void PressureCorrect::Make()
{

    FillPoisson();
    
    for(int i = 0; i < elements; i++){
        id[i] = row2[i] + Grid.cNx*Grid.cNy*col2[i];
    }
         
    insertionSort(id, col2, row2, val2, elements);

    add_ptr(ptr_col, Grid.cNx*Grid.cNy , col2, elements);  

            
    for(int i = 0; i < elements; i++){
        std::cout << id[i] << "  " << col2[i] << " " << row2[i] << std::endl;
    }


    for(int i = 0; i < Grid.cNx*Grid.cNy; i++){

        std::cout << ptr_col[i] << std::endl;
    }
   

    Coord_Mat_double C(
                       Grid.cNx*Grid.cNy,
                       Grid.cNx*Grid.cNy,
                       elements,
                       val2,
                       row2,
                       col2
                      );

    A = C;

    //Matrix de Hessenberg
    MATRIX_double Hr(restart+1,restart, 0.0); 

    H = Hr;
    
}

void PressureCorrect::Solve()
{
    Div();
    


/*
    CompRow_ILUPreconditioner_double Mr(A);


    for(int i= 0; i < x.size(); i++){
         x(i) = (double)i;
    } 


     double *val = new real[81];
     int *row = new int[81];
     int *col = new int[82];

    int checa = 0;

    for(int i= 0; i < Grid.cNx*Grid.cNy; i++){

            val[i] = 2.0;
            row[i] = i;
            col[i] = i;
    }

    col[Grid.cNx*Grid.cNy] = Grid.cNx*Grid.cNy;

    CompCol_Mat_double PP(
                       81,
                       81,
                       Grid.cNx*Grid.cNy,
                       val,
                       row,
                       col
                      );

    VECTOR_double XX(Grid.cNx*Grid.cNy), YY(Grid.cNx*Grid.cNy);

    CompCol_ILUPreconditioner_double MM(PP);
    CompCol_ILUPreconditioner_double M(A);


    for(int i= 0; i < x.size(); i++){
        XX(i) =(double)i* (double)i + (double)i*2.5;
        std::cout << col[i] << "  " << row[i] << std::endl; 
    } 

    YY = XX;
    
    MATRIX_double HH(restart+1,restart, 0.0); 
     GMRES(PP,YY,XX,MM,HH,restart,maxit,tol);*/
    
   /* 
    std::cout << H << std::endl;
    std::cout << XX << std::endl;
    std::cout << YY << std::endl;
    std::cout << b.size()  << std::endl;
    std::cout << b  << std::endl;*/

/*    std::cout << YY << std::endl;
    std::cout << maxit << std::endl;
*/
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
