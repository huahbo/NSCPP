/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*             ********   ***                                 SparseLib++    */
/*          *******  **  ***       ***      ***               v. 1.5c        */
/*           *****      ***     ******** ********                            */
/*            *****    ***     ******** ********              R. Pozo        */
/*       **  *******  ***   **   ***      ***                 K. Remington   */
/*        ********   ********                                 A. Lumsdaine   */
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
/*                                                                           */
/*                                                                           */
/*                     SparseLib++ : Sparse Matrix Library                   */
/*                                                                           */
/*               National Institute of Standards and Technology              */
/*                        University of Notre Dame                           */
/*              Authors: R. Pozo, K. Remington, A. Lumsdaine                 */
/*                                                                           */
/*                                 NOTICE                                    */
/*                                                                           */
/* Permission to use, copy, modify, and distribute this software and         */
/* its documentation for any purpose and without fee is hereby granted       */
/* provided that the above notice appear in all copies and supporting        */
/* documentation.                                                            */
/*                                                                           */
/* Neither the Institutions (National Institute of Standards and Technology, */
/* University of Notre Dame) nor the Authors make any representations about  */
/* the suitability of this software for any purpose.  This software is       */
/* provided ``as is'' without expressed or implied warranty.                 */
/*                                                                           */
/*
*/
/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "gmres.h"
#include <iostream>
#include <stdlib.h>

#include "coord_double.h"
#include "ilupre_double.h"

#include "iotext_double.h"


#include MATRIX_H //Para la matriz de Heissenberg
#define SIZE 6
using namespace std;

int main(int argc, char * argv[])
{

    VECTOR_double x(SIZE), b(SIZE);

    double tol = 1.e-6;
    int maxit = 1   ,
        restart = 32;
    int result;

    
   int Nx = 9,
       Ny = 9; 

   int elements = 13*Nx*Ny - 12*(Nx+Ny);
   
 //   double val2[elements],
  //        row2[elements],
    //      column2[elements];

   double *val2;
   int *row2, *col2;

   val2 = (double*)malloc(sizeof(double)*(elements));
   row2 = (int*)malloc(sizeof(int)*(elements));
   col2 = (int*)malloc(sizeof(int)*(elements));

    

    double val[] = {2.0,-1.0,-1.0,2.0,-1.0,-1.0,2.0,-1.0,-1.0,2.0,-1.0,-1.0,2.0,
                   -1.0,-1.0,2.0};
    int row[] = {0,1,0,1,2,1,2,3,2,3,4,3,4,5,4,5};
 int column[] = {0,0,1,1,1,2,2,2,3,3,3,4,4,4,5,5};



    for(int i = 0; i < SIZE ; i++)
    {
        b(i) = (double)i;
    }
    
    Coord_Mat_double A(6,6,16,val,row,column);


    CompRow_ILUPreconditioner_double Mr(A);

    MATRIX_double H(restart+1,restart, 0.0); //Matrix de Hessenberg 
    MATRIX_double PP(Nx*Ny,Nx*Ny, 0.0); //Matrix de Hessenberg 
    MATRIX_double TT(Nx*Ny,Nx*Ny, 0.0); //Matrix de Hessenberg 

    GMRES(A,x,b,Mr,H,restart,maxit,tol);
  

    int memplace = 0;

    for(int i = 0; i < Nx; i++){
       for(int j =  0; j < Ny; j++){
          int k = i + j*Nx;
          
          PP(k,k) = 7.0; //main banded

          val2[memplace] = 7.0;
          row2[memplace] = k;
          col2[memplace] = k;
          memplace++;

       }
    }   



    for(int i = 3; i < Nx - 3; i++){
       for(int j =  0; j < Ny; j++){
          int k = i + j*Nx;
          
           
          PP(k,k+1) = 1.0; //RS-C1

          val2[memplace] = 1.0;
          row2[memplace] = k;
          col2[memplace] = k+1;
          memplace++;

          PP(k,k-1) = 1.0; //LS-C1

          val2[memplace] = 1.0;
          row2[memplace] = k;
          col2[memplace] = k-1;
          memplace++;

          PP(k,k+2) = 2.0; //RS-C2
	
          val2[memplace] = 2.0;
          row2[memplace] = k;
          col2[memplace] = k+2;
          memplace++;

          PP(k,k-2) = 2.0; //LS-C2

          val2[memplace] = 2.0;
          row2[memplace] = k;
          col2[memplace] = k-2;
          memplace++;


          PP(k,k+3) = 3.0; //RS-C3

          val2[memplace] = 3.0;
          row2[memplace] = k;
          col2[memplace] = k+3;
          memplace++;

          PP(k,k-3) = 3.0; //LS-C3

          val2[memplace] = 3.0;
          row2[memplace] = k;
          col2[memplace] = k-3;
          memplace++;

       }
    }

    for(int j =  0; j < Ny; j++){
       int k = j*Nx;
       PP(k,k+1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;


       PP(k,k+2) = 2.0;

        val2[memplace] = 2.0;
        row2[memplace] = k;
        col2[memplace] = k+2;
        memplace++;


       PP(k,k+3) = 3.0;

       val2[memplace] = 3.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;

        
       k = 1 + j*Nx;
       PP(k,k+1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k+3) = 3.0;

       val2[memplace] = 3.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;
       
       k = 2 + j*Nx;
       PP(k,k+1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k-2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;


       PP(k,k+3) = 3.0;

       val2[memplace] = 3.0;
       row2[memplace] = k;
       col2[memplace] = k+3;
       memplace++;

       k = (Nx-1) + j*Nx;
       PP(k,k-1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;


       PP(k,k-2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = 3.0;

       val2[memplace] = 3.0;
       row2[memplace] = k;
       col2[memplace] = k-3;
       memplace++;
        
       k = (Nx -2)  + j*Nx;
       PP(k,k+1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k-2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = 3.0;

       val2[memplace] = 3.0;
       row2[memplace] = k;
       col2[memplace] = k-3;
       memplace++;
       
       k = (Nx - 3) + j*Nx;
       PP(k,k+1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k+1;
       memplace++;

       PP(k,k-1) = 1.0;

       val2[memplace] = 1.0;
       row2[memplace] = k;
       col2[memplace] = k-1;
       memplace++;

       PP(k,k+2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k+2;
       memplace++;

       PP(k,k-2) = 2.0;

       val2[memplace] = 2.0;
       row2[memplace] = k;
       col2[memplace] = k-2;
       memplace++;

       PP(k,k-3) = 3.0;

       val2[memplace] = 3.0;
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
          PP(k,k+1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k-1*Nx) = 4.0; //US-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

          PP(k,k+2*Nx) = 5.0; //DS-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;

          PP(k,k-2*Nx) = 5.0; //US-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k+3*Nx) = 6.0; //DS-C3

          val2[memplace] = 6.0;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;

          PP(k,k-3*Nx) = 6.0; //US-C3

          val2[memplace] = 6.0;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;

       }
   }         



    //  Y Boundary Values 
    for(int i = 0; i < Nx; i++){
       for(int j =  0; j < 3; j++){

          int k = i + j*Nx;
          PP(k,k+1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = 5.0; //DS-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;


          PP(k,k+3*Nx) = 6.0; //DS-C3

          val2[memplace] = 6.0;
          row2[memplace] = k;
          col2[memplace] = k+3*Nx;
          memplace++;
      }


      {
          int k = i + 1*Nx;
          PP(k,k-1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

      }  


      {
          int k = i + 2*Nx;
          PP(k,k-1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;

          PP(k,k-2*Nx) = 5.0; //DS-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;
      }  

      
   }         


    //  Y Boundary Values 
    for(int i = 0; i < Nx; i++){
       for(int j =  Ny-3 ; j < Ny; j++){

          int k = i + j*Nx;
          PP(k,k-1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k-1*Nx;
          memplace++;


          PP(k,k-2*Nx) = 5.0; //DS-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k-2*Nx;
          memplace++;

          PP(k,k-3*Nx) = 6.0; //DS-C3

          val2[memplace] = 6.0;
          row2[memplace] = k;
          col2[memplace] = k-3*Nx;
          memplace++;
       }


      {
          int k = i + (Ny-2)*Nx;
          PP(k,k+1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;
      }  


      {
          int k = i + (Ny-3)*Nx;
          PP(k,k+1*Nx) = 4.0; //DS-C1

          val2[memplace] = 4.0;
          row2[memplace] = k;
          col2[memplace] = k+1*Nx;
          memplace++;

          PP(k,k+2*Nx) = 5.0; //DS-C2

          val2[memplace] = 5.0;
          row2[memplace] = k;
          col2[memplace] = k+2*Nx;
          memplace++;
      }  
   }    


   for(int i = 0; i < elements; i++){
      TT(row2[i] , col2[i]) = val2[i];       
   }

    Coord_Mat_double SPARS(Nx,Ny,elements,val2,row2,col2);

    cout << A << H << endl << endl;
    cout << tol << "   "<< maxit << endl;
    cout << b << endl << endl;
    cout << x << endl << endl;
    cout << PP << endl << endl;
    cout << SPARS << endl << endl;
    cout << elements << memplace << endl;
    cout << TT << endl;

  return 0;
}
