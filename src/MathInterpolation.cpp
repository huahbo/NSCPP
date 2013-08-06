#include"thomas_tem.hpp"
#include "MathInterpolation.h"


void MathInterpolation::setVal(){
		Alpha = 1.0 / 6.0;
	 	a = 64.0 / 96.0;
	 	b = (0.0);
	 	alphaBound = 1.0;
	 	e1 = 0.25;
	 	e2 = 3.0 / 2.0;
	 	e3 = 0.25;
	 	e4 = 0.0;

}

void MathInterpolation::Resizes(int N){
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


void MathInterpolation::computeInt(Array<real,1> &A ,Array<real,1> &B){

	int size = B.size() - 1;
	A.resize(size);
 	for (int i = 1; i < size -1; i++) {
 		A(i) = a * (B(i + 1) + B(i )) + b * (B(i + 2) + B(i - 1));
 	}

 	A(0) = e1*B(0) + e2*B(1) + e3*B(2) + e4* B(3);
 	A(size - 1) = e4* B(size - 3) + e3*B(size - 2) + e2*B(size - 1) + e1* B(size);

 }


Array<real,2> MathInterpolation::J(){
	Array<real,1> A, A_I, B;

	A.resize(Field.ubound(secondDim));
	A_I.resize(Field.ubound(secondDim));
	B.resize(Field.ubound(secondDim)+1);

	Resizes(Field.ubound(secondDim));

	for(int i = 0; i <= Field.ubound(firstDim); i++){
		B(Range::all()) = Field(i,Range::all());
		computeInt(A,B);
		thomas(bandA,mainBand,bandC,A,A_I);
		Field_J(i,Range::all()) = A_I(Range::all());
	}

	return Field_J;
}




Array<real,2> MathInterpolation::I(){
	Array<real,1> A, A_I, B;

	A.resize(Field.ubound(firstDim));
	A_I.resize(Field.ubound(firstDim));
	B.resize(Field.ubound(firstDim)+1);

	Resizes(Field.ubound(firstDim));
	for(int j = 0; j <= Field.ubound(secondDim); j++){
		B(Range::all()) = Field(Range::all(),j);
		computeInt(A,B);
		thomas(bandA,mainBand,bandC,A,A_I);
		Field_I(Range::all(),j) = A_I(Range::all());
	}
	return Field_I;
}
