#ifndef _PRESSURE_HPP
#define _PRESSURE_HPP

#include<iostream>
#include <blitz/array.h>

//Poisson Equation IML++
#include "gmres.h"
#include "coord_double.h"
#include "compcol_double.h"
#include "comprow_double.h"
#include "iotext_double.h"
#include "ilupre_double.h"
#include "icpre_double.h"
#include "diagpre_double.h"
#include "mgmres.hpp"
#include MATRIX_H

//NS++ 
#include "Mesh.h"
#include "presicion.h"

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

    std::cout << "Building Pressure Equation " << std::endl;
    /**
     Gives the Pressure size to the Divergence
     all this equation is discretized in the Cell place.
    **/

        D.resize(pressure.shape()); 
        elements = 13*Grid.cNx*Grid.cNy - 12*(Grid.cNx+Grid.cNy);

        val2 = new real[elements];
        RHS = new real[Grid.cNx*Grid.cNy];
        x_estimate = new real[Grid.cNx*Grid.cNy];
        row2 = new int[elements];
        col2 = new int[elements];
        ptr_col = new int[Grid.cNx*Grid.cNy + 1];
        id = new int[elements];

        std::cout << "elements: "<< elements << std::endl;

        maxit = 1;
        restart = 32;
        tol = 1.e-12;

    }

    void Div(real);
    void Make();
    void FillPoisson();
    void Solve(real);
    //void Update();

private:
    //Velocities located in the staggered place.
    const Array<real,2>& uS;
    const Array<real,2>& vS;

    Array<real,2>& pressure; 

    //these arrays contain the Ghost Boundary nodes.   
    const Array<real,2>& uBoundary;    
    const Array<real,2>& vBoundary;    
    const MeshTool::MeshBlock & Grid; 
    int elements; 
    
    //Poisson equation Matrix definition
    real *val2;
    int *row2;
    int *col2;
    int *ptr_col;
    int *id;
    real *RHS;
    real *x_estimate;

    int maxit, restart;
    real tol;


    //Divergence Array
    Array<real,2> D;
};

#endif
