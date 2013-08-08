#ifndef _CELLS_H_
#define _CELLS_H_

#include<blitz/array.h>

#include "presicion.h"
#include "MathInterpolation.h"
#include "Mesh.h"


using namespace blitz;

class StaggCellX
{
public:

    StaggCellX(
               Array<real,2> &F_tmp, 
               const Array<real,2>& Bound_tmp, 
               MeshTool::MeshBlock& Grid_tmp
              )
    : FCell(F_tmp), 
      Val(F_tmp),
      Grid(Grid_tmp), 
      IJ(F_IJ),
      Boundary(Bound_tmp)
    {
      F_IJ.resize(Grid.Nx, Grid.Ny);
    }


    void Solve();
    void Dirichlet();
    const Array<real,2>& IJ;
    Array<real,2> &Val; 
     
    

private:
    MathInterpolation FCell;
    MeshTool::MeshBlock& Grid;
    Array<real,2> F_IJ; 
    const Array<real,2>& Boundary; 
};


class StaggCellY
{
public:

    StaggCellY(
               Array<real,2> &F_tmp, 
               const Array<real,2>& Bound_tmp, 
               MeshTool::MeshBlock& Grid_tmp
              )
    : FCell(F_tmp), 
      Val(F_tmp),
      Grid(Grid_tmp), 
      IJ(F_IJ),
      Boundary(Bound_tmp)
    {
      F_IJ.resize(Grid.Nx, Grid.Ny);
    }


    void Solve();
    void Dirichlet();
    const Array<real,2>& IJ;
    Array<real,2> &Val; 

private:
    MathInterpolation FCell;
    MeshTool::MeshBlock& Grid;
    Array<real,2> F_IJ; 
    const Array<real,2>& Boundary; 
};

#endif
