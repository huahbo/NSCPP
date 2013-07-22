/*
 * CompactSchemes.hpp
 *
 *  Created on: 22/03/2013
 *      Author: cipia
 */

//---------------------------------------------------------------------------//
// NS-Solver                                                                 //
// Numerical Solver to Navier-Stokes Equation for Incompressible Flow        //
// using Compact Schemes and the Poisson Equation to the Pressure Correction //
//---------------------------------------------------------------------------//
// Project built by Jorge Ramírez Cruz(JramirezCr@iingen.unam.mx             //
// Instituto de Ingeniería UNAM                                              //
//---------------------------------------------------------------------------//
// File Description:Class to derivate with Compact Schemes                   //                                                 //
//---------------------------------------------------------------------------//


#ifndef COMPACTSCHEMES_HPP_
#define COMPACTSCHEMES_HPP_
#include<iostream>
#include"NS_Array.hpp"
#include"thomas_tem.hpp"
using namespace NS_Solver;
using namespace std;

template<typename Tprec>
class Compact {
public:
	Compact(Array<Tprec, 3> &, Array<Tprec, 1> &,Tprec);

	inline void CenteredInside(Array<Tprec,1> & ,Array<Tprec,1> &, Tprec);  //Base
	inline void Resizes(int );                                              //Base
	const Array<Tprec,3> getDx(){return dX;}//Base 
	const Array<Tprec,3> getDy(){return dY;}//Base
	const Array<Tprec,3> getDz(){return dZ;}//Base
	inline void Centered();



protected:
	Array<Tprec, 3> Field, dX, dY, dZ;
	Array<Tprec, 1> deltas;
	Tprec Alpha;
	Tprec alphaBound;
	int Nx, Ny, Nz;
	Array<Tprec,1> mainBand, bandA,bandC;
	Tprec a, b, e1, e2, e3, e4;

};

//Constructor for 3d-Field
template<typename Tprec>
Compact<Tprec>::Compact(Array<Tprec, 3> &P, Array<Tprec, 1> &Ds,
			Tprec Aph){
		//Resize 3-dimensional Arrays
		Field.resize(P.shape());
		dX.resize(P.shape());
		dY.resize(P.shape());
		dZ.resize(P.shape());
		//Resize deltas' array
        deltas.resize(Ds.shape());
		//Array size
		Nx = Field.ubound(firstDim) + 1;
		Ny = Field.ubound(secondDim) + 1;
	    Nz = Field.ubound(thirdDim) + 1;

		Field = P;
		deltas = Ds;
		Alpha = Aph;

		Centered();
	}


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
 inline void Compact<Tprec>::Centered(){
	 	a = (3.0 / 2.0 ) / 2.0;
	 	b = (0.0) / 4.0;
	 	alphaBound = 3.0;
	 	e1 = - 16.0 / 6.0;
	 	e2 = 3.0 / 2.0;
	 	e3 = 3.0 / 2.0;
	 	e4 = -1.0 / 6.0;

 }


 /*------------------------------------------------------------------------------------------------
  * Sub Classes- Discretization(Staggered or Centered)
  * -----------------------------------------------------------------------------------------------
  * Sub Class: Staggered
  -----------------------------------------------------------------------------------------------*/


template<typename Tprec>
class Staggered:
	public Compact<Tprec>
	{
public:
	Staggered(Array<Tprec, 3> &P, Array<Tprec, 1> &Ds, Tprec Aph)
	: Compact<Tprec>(P,Ds,Aph)
	{
	}
	Staggered &Dx();
	Staggered &Dy();
	Staggered &Dz();
};


template<typename Tprec>
Staggered<Tprec> &Staggered<Tprec>::Dx(){
	 Array<Tprec,1> R(Compact<Tprec>::Nx), Rsolved(Compact<Tprec>::Nx),dev(Compact<Tprec>::Nx);

	 Compact<Tprec>::Resizes(Compact<Tprec>::Nx);
	 for(int j = 0; j < Compact<Tprec>::Ny; j++ ){
		 for(int k = 0; k < Compact<Tprec>::Nz; k++){
			 R = Compact<Tprec>::Field(Range::all(), j, k);
			 Compact<Tprec>::CenteredInside(Rsolved,R,Compact<Tprec>::deltas(0));
			 thomas(Compact<Tprec>::bandA,Compact<Tprec>::mainBand,Compact<Tprec>::bandC,Rsolved,dev);
			 Compact<Tprec>::dX(Range::all(),j,k) =  dev(Range::all());
		 }
	 }

	 return *this;
}

template<typename Tprec>
Staggered<Tprec> &Staggered<Tprec>::Dy(){
	 Array<Tprec,1> R(Compact<Tprec>::Ny), Rsolved(Compact<Tprec>::Ny),dev(Compact<Tprec>::Ny);

	 Compact<Tprec>::Resizes(Compact<Tprec>::Ny);

	 for(int i = 0; i < Compact<Tprec>::Nx; i++ ){
		 for(int k = 0; k < Compact<Tprec>::Nz; k++){
			 R = Compact<Tprec>::Field(i,Range::all(), k);
			 Compact<Tprec>::CenteredInside(Rsolved,R,Compact<Tprec>::deltas(0));
			 thomas(Compact<Tprec>::bandA,Compact<Tprec>::mainBand,Compact<Tprec>::bandC,Rsolved,dev);
			 Compact<Tprec>::dY(i,Range::all(),k) =  dev(Range::all());
		 }
	 }
	 return *this;
}

template<typename Tprec>
Staggered<Tprec> &Staggered<Tprec>::Dz(){
	 Array<Tprec,1> R(Compact<Tprec>::Nz), Rsolved(Compact<Tprec>::Nz),dev(Compact<Tprec>::Nz);

	 Compact<Tprec>::Resizes(Compact<Tprec>::Nz);
	 for(int i = 0; i < Compact<Tprec>::Nx; i++ ){
		 for(int j = 0; j < Compact<Tprec>::Ny; j++){
			 R = Compact<Tprec>::Field(i,j,Range::all());
			 Compact<Tprec>::CenteredInside(Rsolved,R,Compact<Tprec>::deltas(0));
			 thomas(Compact<Tprec>::bandA,Compact<Tprec>::mainBand,Compact<Tprec>::bandC,Rsolved,dev);
			 Compact<Tprec>::dZ(i,j,Range::all()) =  dev(Range::all());
		 }
	 }
	 return *this;
}

/*-------------------------------------------------------------------------------------------------
 * Sub Class: Centered
 ------------------------------------------------------------------------------------------------*/
template<typename Tprec>
class Centered: public Compact<Tprec>{

public:
	Centered(Array<Tprec,3> &P, Array<Tprec,1> &Ds, Tprec Aph, Tprec AphMean)
	: Compact<Tprec>(P,Ds,Aph)
	  {
		AlphaMean = AphMean;
		Center();
	  }

	inline void Mean();
	inline void Center();
	inline void ResizesMean(int);
	inline void CenteredMean(Array<Tprec,1> & ,Array<Tprec,1> &);
	Centered &Dx();

private:
	Tprec AlphaMean;
	Tprec alphaBoundMean;
	Array<Tprec,1> mainBandMean, bandAMean, bandCMean;
	Tprec aMean, bMean, e1Mean, e2Mean, e3Mean, e4Mean;

};

template<typename Tprec>
inline void Centered<Tprec>::Center(){

	 	Compact<Tprec>::a = 12.0 / 11.0;
	 	Compact<Tprec>::b = (0.0);
	 	Compact<Tprec>::alphaBound = 22.0;
	 	Compact<Tprec>::e1 = -577.0/24.0;
	 	Compact<Tprec>::e2 = 603.0/24.0;
	 	Compact<Tprec>::e3 = -27.0/24.0;
	 	Compact<Tprec>::e4 = 577.0/24.0;

}

template<typename Tprec>
inline void Centered<Tprec>::Mean(){

	 	aMean = 12.0 / 11.0;
	 	bMean = (0.0);
	 	alphaBoundMean = 22.0;
	 	e1Mean = 0.25;
	 	e2Mean = 3.0 / 2.0;
	 	e3Mean = 0.25;
	 	e4Mean = 0.0;

}

template<typename Tprec>
inline void Centered<Tprec>::ResizesMean(int N){
    mainBandMean.resize(N);
    bandAMean.resize(N);
    bandCMean.resize(N);


	 mainBandMean(Range::all()) = 1.0;
	 bandAMean(Range::all()) = AlphaMean;
	 bandCMean(Range::all()) = AlphaMean;
	 bandCMean(0) = alphaBoundMean;
	 bandCMean(N-1) = 0.0;
	 bandAMean(0) = 0.0;
	 bandAMean(N-1) = alphaBoundMean;
}

template<typename Tprec>
inline void Centered<Tprec>::CenteredMean(Array<Tprec,1> &A ,Array<Tprec,1> &B){

	int size = B.size() - 1;
	A.resize(size);
 	for (int i = 1; i < size - 1; i++) {
 		A(i) = aMean * (B(i + 1) + B(i - 1)) + bMean * (B(i + 2) + B(i - 2));
 	}

 	A(0) = e1Mean*B(0) + e2Mean*B(1) + e3Mean*B(2) + e4Mean* B(3);
 	A(size - 1) = e4Mean* B(size - 3) + e3Mean*B(size - 2) + e2Mean*B(size - 1) + e1Mean* B(size);
 }

template<typename Tprec>
Centered<Tprec> &Centered<Tprec>::Dx(){
	 Array<Tprec,1> R(Compact<Tprec>::Nx), Rsolved(Compact<Tprec>::Nx),dev(Compact<Tprec>::Nx);
	 Array<Tprec,1> Rmean(Compact<Tprec>::Nx);
	 ResizesMean(Compact<Tprec>::Nx);
	 Compact<Tprec>::Resizes(Compact<Tprec>::Nx);

	 for(int j = 0; j < Compact<Tprec>::Ny; j++ ){
		 for(int k = 0; k < Compact<Tprec>::Nz; k++){

			 Rmean = Compact<Tprec>::Field(Range::all(), j, k);
			 CenteredMean(R,Rmean);

			 thomas(bandAMean,mainBandMean,bandCMean,Rmean,R);

			 Compact<Tprec>::CenteredInside(Rsolved,Rmean,Compact<Tprec>::deltas(0));
			 thomas(Compact<Tprec>::bandA,Compact<Tprec>::mainBand,Compact<Tprec>::bandC,Rsolved,dev);
			 Compact<Tprec>::dX(Range::all(),j,k) =  dev(Range::all());
		 }
	 }

	 return *this;
}


#endif /* COMPACTSCHEMES_HPP_ */
