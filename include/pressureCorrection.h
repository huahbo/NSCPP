#ifndef _PRESSURE_HPP
#define _PRESSURE_HPP

#include <blitz/array.h>

#include "../include/Mesh.h"
#include       "presicion.h"

using blitz::Array;


class PressureCorrect
{
public:
    PressureCorrect(
                    Array<real,2>& vel_u, 
                    Array<real,2>& vel_v,
                    Array<real,2>& press_p,
                    Array<real,2>& Bound_u,
                    Array<real,2>& Bound_v,
                    MeshTool::MeshBlock &Grid_set   
	                )
    : uS(vel_u), 
      vS(vel_v), 
      pressure(press_p),
      uBoundary(Bound_u),
      vBoundary(Bound_v),
      Grid(Grid_set)
    { 

    /**
     Gives the Pressure size to the Divergence
     all this equation is discretized in the Cell place.
    **/


        D.resize(pressure.shape()); 
    }


    void Div();

private:
    //Velocities located in the staggered place.
    const Array<real,2>& uS;
    const Array<real,2>& vS;

    const Array<real,2>& pressure; 

    //these arrays contain the Ghost Boundary nodes.   
    const Array<real,2>& uBoundary;    
    const Array<real,2>& vBoundary;    
    const MeshTool::MeshBlock & Grid; 

    //Divergence Array
    Array<real,2> D;
};

#endif
