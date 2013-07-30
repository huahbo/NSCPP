/*
 * Derivatives.hpp
 *
 *  Created on: 22/04/2013
 *      Author: jrc-ubu
 */

#ifndef DERIVATIVES_HPP_
#define DERIVATIVES_HPP_

#include<iostream>
#include"NS_Array.hpp"
#include"thomas_tem.hpp"
using namespace NS_Solver;
using namespace std;

template<typename Tprec>
class Compact {
public:
	Compact(Array<Tprec,1> &Tdelta)
	: delta(Tdelta)
	{
	};

	inline void Resizes(int );
	inline void Resizes_2(int );
	inline void CenterVal();
	inline void MidVal();
	inline void CenteredInside(Array<Tprec,1> & ,Array<Tprec,1> &, Tprec);
	inline void MidPoints(Array<Tprec,1> & ,Array<Tprec,1> &, Tprec);
	inline void MidPoints_2(Array<Tprec,1> & ,Array<Tprec,1> &, Tprec);

	inline Array<Tprec,2> xS(const Array<Tprec,2> &); //Derivatives at the Integer points
	inline Array<Tprec,2> yS(const Array<Tprec,2> &);

	inline Array<Tprec,2> xM(const Array<Tprec,2> &); //Derivatives at the middle points
	inline Array<Tprec,2> yM(const Array<Tprec,2> &);


private:

	Tprec Alpha;
	Tprec alphaBound, alpha1, alphaN;
	Array<Tprec,1> &delta;
	Array<Tprec,1> mainBand, bandA,bandC;
	Tprec a, b, e1, e2, e3, e4, e5;
	Tprec alpha4, b4;

};


//Compute the Right side of Compact Schemes for Cell Centered
template<typename Tprec>
inline void Compact<Tprec>::CenteredInside(Array<Tprec,1> &A ,Array<Tprec,1> &B, Tprec dh){

	int size = B.size();
 	for (int i = 1; i < size - 1; i++) {
 		A(i) = (a * (B(i + 1) - B(i - 1)) + b * (B(i + 2) - B(i - 2))) / dh;
 	}

 	A(0) = (e1 * B(0) + e2*B(1) + e3*B(2) + e4* B(3)) / dh;
 	A(size - 1) = (-e4* B(size - 4) - e3*B(size - 3) - e2*B(size - 2) - e1* B(size - 1)) / dh;
 }



//Compute the Right side of Compact Schemes for Cell Centered
template<typename Tprec>
inline void Compact<Tprec>::MidPoints(Array<Tprec,1> &A ,Array<Tprec,1> &B, Tprec dh){

	int size = B.size() + 1;
    A.resize(size);

 	for (int i = 2; i <= size - 3; i++) {
 		A(i) = (a * (B(i) - B(i-1)) + b*(B(i+1) -B(i-2))) / dh;
 	}

 	A(1) = ((12.0/11.0)*(B(1) - B(0))) / dh;
 	A(size - 2) = ((12.0/11.0)*(B(size - 2) - B(size - 3))) / dh;


 	A(0) = (e1 *B(0) + e2*B(1) + e3*B(2) + e4* B(3) + e5*B(4)) /(dh);
 	A(size - 1) = (-e5*B(size -6)-(e4)* B(size- 5) - (e3)*B(size-4) - (e2)*B(size-3)  -(e1)* B(size -2)) / (dh);
 }


 template<typename Tprec>
 inline void Compact<Tprec>::Resizes(int N){
     mainBand.resize(N);
     bandA.resize(N);
     bandC.resize(N);


	 mainBand(Range::all()) = 1.0;
	 bandA(Range::all()) = Alpha;
	 bandC(Range::all()) = Alpha;
	 bandC(0) = alphaBound;
	 bandC(N-1) = 0.0;

	 bandA(0) = 0.0;
	 bandA(N-1) = alphaBound;
}

 template<typename Tprec>
  inline void Compact<Tprec>::Resizes_2(int N){
      mainBand.resize(N);
      bandA.resize(N);
      bandC.resize(N);


 	 mainBand(Range::all()) = 1.0;
 	 bandA(Range::all()) = Alpha;
 	 bandC(Range::all()) = Alpha;

 	 bandA(0) = 0.0;
 	 bandA(1) = alpha4;
 	 bandA(N-1) = alphaBound;
 	 bandA(N-2) = alpha4;

 	 bandC(0) = alphaBound;
 	 bandC(1) = alpha4;
  	 bandC(N-1) = 0.0;
  	 bandC(N-2) = alpha4;

  }



 template<typename Tprec>
 inline void Compact<Tprec>::CenterVal(){
		Alpha = 0.25;
	 	a = (3.0 / 2.0 ) / 2.0;
	 	b = (0.0) / 4.0;
	 	alphaBound = 3.0;
	 	e1 = - 17.0 / 6.0;
	 	e2 = 3.0 / 2.0;
	 	e3 = 3.0 / 2.0;
	 	e4 = -1.0 / 6.0;

 }

 template<typename Tprec>
 inline void Compact<Tprec>::MidVal(){

		Alpha = 9.0 / 62.0;
 	 	a = 63.0 / 62.0;
 	 	b = (17.0 / 186.0);
 	 	alphaBound = 331.0 / 15.0;
 	 	alpha4 = 1.0 / 22.0;
 	 	e1 = -8677.0 / 360.0;
 	 	e2 =  4531.0 / 180.0;
 	 	e3 = -11.0 / 10.0;
 	 	e4 = 1.0 / 36.0;
 	 	e5 = 1.0 / 360.0;
 }


 template<typename Tprec>
 Array<Tprec,2> Compact<Tprec>::xM(const Array<Tprec,2> &Field){

	 int Nx = Field.ubound(firstDim) + 1;
	 int Ny = Field.ubound(secondDim) + 1;
 	 Array<Tprec, 1> R(Nx), Rsolved(Nx), dev(Nx);
 	 Array<Tprec, 2> ResDev(Nx,Ny);
 	 CenterVal();
 	 Resizes(Nx);


 	 for(int j = 0; j < Ny; j++ ){

 			 R = Field(Range::all(), j);
 			 CenteredInside(Rsolved,R, delta(0));
 			 thomas(bandA,mainBand,bandC,Rsolved,dev);
 			 ResDev(Range::all(),j) = dev;
 	 }

 	 return ResDev;
 }


 template<typename Tprec>
 Array<Tprec,2> Compact<Tprec>::yM(const Array<Tprec,2> &Field){

	 int Nx = Field.ubound(firstDim) + 1;
	 int Ny = Field.ubound(secondDim) + 1;
 	 Array<Tprec, 1> R(Ny), Rsolved(Ny), dev(Ny);
 	 Array<Tprec, 2> ResDev(Nx,Ny);
 	 CenterVal();
 	 Resizes(Ny);


 	 for(int i = 0; i < Nx; i++ ){

 			 R = Field(i,Range::all());
 			 CenteredInside(Rsolved,R, delta(1));
 			 thomas(bandA,mainBand,bandC,Rsolved,dev);
 			 ResDev(i,Range::all()) = dev;
 	 }

 	 return ResDev;
 }

 template<typename Tprec>
 Array<Tprec,2> Compact<Tprec>::xS(const Array<Tprec,2> &Field){

	 int Nx = Field.ubound(firstDim) + 2;
	 int Ny = Field.ubound(secondDim) + 1;
 	 Array<Tprec, 1> R(Nx - 1), Rsolved(Nx), dev(Nx);
 	 Array<Tprec, 2> ResDev(Nx,Ny);
 	 MidVal();
 	 Resizes_2(Nx);


 	 for(int j = 0; j < Ny ; j++ ){

 			 R = Field(Range::all(), j);
 			 MidPoints(Rsolved,R, delta(0));
 			 thomas<Tprec>(bandA,mainBand,bandC,Rsolved,dev);
 			 ResDev(Range::all(),j) = dev(Range::all());
 	 }

 	 return ResDev;





 }

template<typename Tprec>
Array<Tprec,2> Compact<Tprec>::yS(const Array<Tprec,2> &Field){
	int Nx = Field.ubound(firstDim) + 1,
		Ny = Field.ubound(secondDim) + 2;

	Array<Tprec,2> ResDev(Nx,Ny);
	Array<Tprec,1> R(Ny - 1), Rsolved(Ny), dev(Ny);
	MidVal();
	Resizes_2(Ny);
	for(int i = 0; i < Nx;i++){
		R = Field(i,Range::all());
		MidPoints(Rsolved,R,delta(1));
		thomas(bandA,mainBand,bandC,Rsolved,dev);
		ResDev(i,Range::all()) = dev(Range::all());
	}



	return ResDev;
}


#endif /* DERIVATIVES_HPP_ */
