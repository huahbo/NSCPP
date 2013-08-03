#include<blitz/array.h>
#include"thomas_tem.hpp"

using namespace blitz;

	inline void setValues();
	inline void Resizes(int N);
	inline void computeInt(Array<Tprec,1> & ,Array<Tprec,1> &);

	inline Array<Tprec,2> J(Array<Tprec,2>& );
	inline Array<Tprec,2> I(Array<Tprec,2>& );


	Array<Tprec,2> &F,F_P , F_IJ;
	Array<Tprec,2> UnkBound;
	Tprec Alpha;
	Tprec alphaBound;
	Array<Tprec,1> mainBand, bandA, bandC;
	Tprec a, b, e1, e2, e3, e4;


inline MathInterpolation::setValues(){
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
