/*
 * prueba.cpp
 *
 *  Created on: 20/03/2013
 *      Author: cipia
 */


   #include<iostream>
   #include<stdlib.h>
   #include<stdio.h>
   #include<math.h>

   #include "thomas_tem.hpp"
   #include<blitz/array.h>
   #include "NS_Array.hpp"
   #include "NS_Array.hpp"





//#include"CompactSchemes.hpp"
#include"Interpol.hpp"
#include"Fluid.hpp"
#include"Momentum.hpp"
#include"Derivatives.hpp"
#include"explicitDeriv.hpp"


   using namespace std;
   using namespace blitz;
   using namespace NS_Solver;

   template<typename Tprec>
   void imprimo(Array<Tprec,1>& Bobo)
   {
   	cout<< Bobo;
   }



int main(){





	  int nx = 100;
	  int ny = 6;

	//Array<double,3> A(nx,1,1),dek(nx,1,1);

	Array<double,1> B(2);

	Array<double,2> C(nx,ny-1);
	Array<double,2> D(nx-1,ny);
	Array<double,2> Dedo(nx-1,ny);
	Array<double,2> F(nx-1,ny-1);
	Array<double,2> V_X_Bound(2,nx), U_Y_Bound(2,ny);




	Array<double,1> deltas(2);

//    NS_Array<double,1>::Matrix x(nx), wsen(nx),wcos(nx),a(nx),b(nx),c(nx),r(nx);

    double Lx = 2.0*3.1415926;
    double deltaX, deltaY;

    deltas(0) = Lx / (double)(nx-1);
    deltas(1) = Lx / (double)(ny-1);


     for(int i = 0; i <= C.ubound(firstDim); i++){
    	 for(int j= 0; j <= C.ubound(secondDim); j++){
    		 C(i,j) = sin(deltas(0)*(double)i);
    	 }
    	 V_X_Bound(0,i) = sin(deltas(0)*(double)i);
    	 V_X_Bound(1,i) = sin(deltas(0)*(double)i);

     }

     for(int j= 0; j <= D.ubound(secondDim); j++){
         for(int i = 0; i <= D.ubound(firstDim); i++){
    		 D(i,j) = sin(deltas(1)*(double)j);
    		 Dedo(i,j) =2.0* sin(deltas(1)*(double)j) + cos(deltas(1)*(double)j);
    	 }
    	 U_Y_Bound(0,j) = sin(deltas(1)*(double)j);
    	 U_Y_Bound(1,j) = sin(deltas(1)*(double)j);

     }

     for(int j= 0; j <= F.ubound(secondDim); j++){
         for(int i = 0; i <= F.ubound(firstDim); i++){
    		 F(i,j) = sin(deltas(0)*(double)i);
    	 }
     }



     Interpolation<double,1> uCell(C,V_X_Bound);
     Interpolation<double,2> vCell(D,U_Y_Bound);
     Fluid<double> weight(C,C,C);
     //uCell.EW();
     Momentum<double,1> X(C,D,F,uCell,vCell,deltas);
     Momentum<double,2> Y(C,D,F,uCell,vCell,deltas);

     Compact<double> d(deltas);



        //uCell.C();


    // cout << C;

     //cout << D;
       //cout << uCell.P();
     //cout << uCell.IJ();
     //cout << vCell.IJ();
     /*
       cout << uCell.IJ();


      //cout << vCell.IJ();
      //cout << X.FluxF();
      //cout << X.FluxG();

      //cout << Y.FluxF();
      //cout << Y.FluxG();




       cout << X.DivFluxG();
      //cout << d.xS(C);

      //cout << X.DivFluxG();
     //cout << D;
     //cout << d.yM(vCell.Cell());

      //cout << X.DivFluxFv();
     //cout << Y.DivFluxGv();
     //cout << Y.DivPressure();
*/


    // cout << uCell.P();

     //uCell.P()*uCell.P();
     cout << X.DivFluxF();
     cout << X.DivPressure();

     cout << Y.DivFluxG();

     //cout << d.xS(C);



     return 0;

   }





