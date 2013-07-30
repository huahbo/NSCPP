/*
 * Interpol.hpp
 *
 *  Created on: 10/04/2013
 *      Author: jrc-ubu
 */

#ifndef INTERPOL_HPP_
#define INTERPOL_HPP_

#include<blitz/array.h>
#include"thomas_tem.hpp"

using namespace blitz;

template<int Zone>
struct typeDim
{
	enum{staggeredZone = Zone};
};

template<typename Tprec, int Zone>
class Interpolation{

private: //This Private is for the Partial Specialization

	//CELL CENTER
	inline void P_Ground(typeDim<1>){
		F_P.resize(F.ubound(firstDim),F.ubound(secondDim)+1);
		F_P = I(F);
	}

	inline void P_Ground(typeDim<2>){
		F_P.resize(F.ubound(firstDim)+1,F.ubound(secondDim));
		F_P = J(F);
	}


	//CORNERS
	inline void IJ_Ground(typeDim<1>){
		F_IJ.resize(F.ubound(firstDim)+1,F.ubound(secondDim)+2);
		F_IJ(Range::all(),Range(1,F.ubound(secondDim))) = J(F);
	    F_IJ(Range::all(),0) = UnkBound(0,Range::all());
		F_IJ(Range::all(), F_IJ.ubound(secondDim)) = UnkBound(1,Range::all());
	}

	inline void IJ_Ground(typeDim<2>){
		F_IJ.resize(F.ubound(firstDim)+2,F.ubound(secondDim)+1);
		F_IJ(Range(1,F.ubound(firstDim)),Range::all()) = I(F);
		//Cambiar lo siguente ppor las fronteras reales
		F_IJ(0,Range::all()) = UnkBound(0,Range::all());
		F_IJ(F_IJ.ubound(firstDim),Range::all()) = UnkBound(1,Range::all());
	}

/*
	inline void Init_Boundary(typeDim<1>,Array<Tprec,2> Unk_Var){
		UnkBound.resize(2,Unk_Var.ubound(secondDim) +1 );
		UnkBound = Unk_Var;
	}

	inline void Init_Boundary(typeDim<2>,Array<Tprec,2> Unk_Var){
			UnkBound.resize(2,Unk_Var.ubound(secondDim) +1 );
			UnkBound = Unk_Var;
	}

*/
public:
	Interpolation(Array<Tprec,2> &Atemp, Array<Tprec,2> Unk_Var)
	: F(Atemp), UnkBound(Unk_Var)   //Pass by reference the Field
	{
		setValues();
		//Init_Boundary(typeDim<Zone>(),Unk_Var);
	}

	inline void setValues();
	inline void Resizes(int N);
	inline void computeInt(Array<Tprec,1> & ,Array<Tprec,1> &);
	inline Array<Tprec,2> P();
	inline Array<Tprec,2> IJ();
	inline Array<Tprec,2> J(Array<Tprec,2>& );
	inline Array<Tprec,2> I(Array<Tprec,2>& );

	inline Array<Tprec,2> Cell();
	inline Array<Tprec,2> getP();
	inline Array<Tprec,2> getIJ();
//	inline void Init_Boundary(typeDim<Zone>(),Array<Tprec,2> );
	inline void setField(Array<Tprec,2> &A)
	{
        	F = A;
	}

	//Specialization Grounds
	inline void P_Ground(typeDim<Zone>());
	inline void IJ_Ground(typeDim<Zone>());
private:

	Array<Tprec,2> &F,F_P , F_IJ;
	Array<Tprec,2> &UnkBound;
	Tprec Alpha;
	Tprec alphaBound;
	Array<Tprec,1> mainBand, bandA, bandC;
	Tprec a, b, e1, e2, e3, e4;
};



//Member Functions Class


template<typename Tprec, int Zone>
inline void Interpolation<Tprec, Zone>::setValues(){
		Alpha = 1.0 / 6.0;
	 	a = 64.0 / 96.0;
	 	b = (0.0);
	 	alphaBound = 1.0;
	 	e1 = 0.25;
	 	e2 = 3.0 / 2.0;
	 	e3 = 0.25;
	 	e4 = 0.0;

}

template<typename Tprec, int Zone>
inline void Interpolation<Tprec,Zone>::Resizes(int N){
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

template<typename Tprec, int Zone>
inline void Interpolation<Tprec, Zone>::computeInt(Array<Tprec,1> &A ,Array<Tprec,1> &B){

	int size = B.size() - 1;
	A.resize(size);
 	for (int i = 1; i < size -1; i++) {
 		A(i) = a * (B(i + 1) + B(i )) + b * (B(i + 2) + B(i - 1));
 	}

 	A(0) = e1*B(0) + e2*B(1) + e3*B(2) + e4* B(3);
 	A(size - 1) = e4* B(size - 3) + e3*B(size - 2) + e2*B(size - 1) + e1* B(size);

 }




template<typename Tprec, int Zone>
inline Array<Tprec,2> Interpolation<Tprec,Zone>::J(Array<Tprec,2> &Base){
	Array<Tprec,2> Base_I(Base.ubound(firstDim)+1,Base.ubound(secondDim));
	Array<Tprec,1> A, A_I, B;

	A.resize(Base.ubound(secondDim));
	A_I.resize(Base.ubound(secondDim));
	B.resize(Base.ubound(secondDim)+1);

	Resizes(Base.ubound(secondDim));

	for(int i = 0; i <= Base.ubound(firstDim); i++){
		B(Range::all()) = Base(i,Range::all());
		computeInt(A,B);
		thomas(bandA,mainBand,bandC,A,A_I);
		Base_I(i,Range::all()) = A_I(Range::all());
	}

	return Base_I;
}




template<typename Tprec, int Zone>
inline Array<Tprec,2> Interpolation<Tprec,Zone>::I(Array<Tprec,2> &Base){
	Array<Tprec,2> Base_I(Base.ubound(firstDim),Base.ubound(secondDim)+1);
	Array<Tprec,1> A, A_I, B;

	A.resize(Base.ubound(firstDim));
	A_I.resize(Base.ubound(firstDim));
	B.resize(Base.ubound(firstDim)+1);

	Resizes(Base.ubound(firstDim));
	for(int j = 0; j <= Base.ubound(secondDim); j++){
		B(Range::all()) = Base(Range::all(),j);
		computeInt(A,B);
		thomas(bandA,mainBand,bandC,A,A_I);
		Base_I(Range::all(),j) = A_I(Range::all());
	}
	return Base_I;
}

template<typename Tprec, int Zone>
inline Array<Tprec,2>Interpolation<Tprec,Zone>::IJ(){
	IJ_Ground(typeDim<Zone>());
	return F_IJ;
}

template<typename Tprec, int Zone>
inline Array<Tprec,2>Interpolation<Tprec,Zone>::P(){
	P_Ground(typeDim<Zone>());
	return F_P;
}

template<typename Tprec, int Zone>
inline Array<Tprec,2>Interpolation<Tprec,Zone>::Cell()
{ return F;}

template<typename Tprec, int Zone>
inline Array<Tprec,2>Interpolation<Tprec,Zone>::getIJ(){
	return F_IJ;
}

template<typename Tprec, int Zone>
inline Array<Tprec,2>Interpolation<Tprec,Zone>::getP(){
	return F_P;
}



#endif /* INTERPOL_HPP_ */
