#ifndef _MATH_INT_H
#define _MATH_INT_H

#include "iostream"
#include "presicion.h"
#include "blitz/array.h"
#include <iostream>

using namespace blitz;

//Array<real,2> AxisXtoCell(const Array<real,2>&, const Array<real,2>&);
//Array<real,2> AxisYtoCell(const Array<real,2>&, const Array<real,2>& );


class MathInterpolation
{
public:
    MathInterpolation(Array<real,2>& pf)
    : Field(pf)
    {   
        /**
        Dimension Field Array
        **/
        xDim = Field.ubound(firstDim) + 1; 
        yDim = Field.ubound(secondDim) + 1; 
        Field_J(Field.ubound(firstDim)+1,Field.ubound(secondDim));
        Field_I(Field.ubound(firstDim),Field.ubound(secondDim)+1);
        setVal();  

    }

    /*
    Sets the Schemes values (order of interpolation)
    */
    void setVal(); 
    inline void Resizes(int N);
    inline void computeInt(Array<real,1> & ,Array<real,1> &);
    inline Array<real,2> I();
    inline Array<real,2> J();

private:     

    Array<real,2> &Field;
    Array<real,2> Field_J,
                  Field_I;
    int xDim, yDim;
    real Alpha, alphaBound, a, b, e1, e2, e3, e4; 
    Array<real,1> mainBand, bandA, bandC;

};



#endif
