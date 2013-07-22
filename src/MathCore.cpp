
#include "../include/MathCore.h"

Array<real,2> AxisXtoCell(
                         Array<real,2> &Field
                         Array<real,2> &Boundary //Ghost Boundary
                         )
{
    /***
    Dimension Field (staggered)
    ***/
    int sNx = Field.ubound(firstDim) + 1;
    int sNy = Field.ubound(secondDim) + 1;


    /**
    Has one dimension more in the Axi-X
    **/
    Array<real,2> fixedField(sNx + 1, sNy);

    for(int i = 1; i < sNx ; i++){
        for(int j = 1; j < sNy ; j++){
 
            fixedField(i,j) = (- 1.0*Field(i+2,j) + 9.0*Field(i+1,j) 
                               + 9.0*Field(i,j) - Field(i-1,j)) / 16.0
        }                                           
    }

    /**
    Compute the Boundaries using Ghost Nodes
    **/
    for(int j = 1; j < sNy ; j++){
        fixedField(0,j) = -1.0*Field(2,j) + 9.0*Field(1,j) 
                          +9.0*Field(
    }

        
 
    return fixedField;

}
