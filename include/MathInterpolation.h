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
        setValues();  

    }

    /*
    Sets the Schemes values (order of interpolation)
    */
    inline void setValues(); 


private:     

    Array<real,2>& Field;
    int xDim, yDim;
    real Alpha, alphaBound, a, b, e1, e2, e3, e4; 

};



#endif
