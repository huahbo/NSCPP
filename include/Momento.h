#ifndef _MOMENTO_H_
#define _MOMENTO_H_

#include "blitz/array.h"
#include "presicion.h"
#include "Cells.h"
#include "Mesh.h"
#include<iostream>

using namespace blitz;

class Momentum
{
public:

    Momentum(
             const StaggCellX& U_tmp,
             const StaggCellY& V_tmp,
             const Array<real,2>& pressure_tmp,
             const MeshTool::MeshBlock& Grid_tmp  
            )
    : Ucell(U_tmp),
      Vcell(V_tmp),
      pressure(pressure_tmp),
      Grid(Grid_tmp)
    {
        std::cout << "Building Momentum Equation" << std::endl; 
    }

protected:
    const StaggCellX& Ucell;
    const StaggCellY& Vcell;
    const Array<real,2>& pressure;
    const MeshTool::MeshBlock& Grid;
    Array<real,2> Divergence,
                  FluxF,
                  FluxG,
                  FluxFv,
                  FluxGv,
                  DivP;
};



class MomentumX
: public Momentum 
{
public:

    MomentumX(
              const StaggCellX&,
              const StaggCellY&,
              const Array<real,2>&,
              const MeshTool::MeshBlock&  
             );

    Array<real,2> Solve();

};


class MomentumY
: public Momentum 
{
public:

    MomentumY(
              const StaggCellX& U_tmp,
              const StaggCellY& V_tmp,
              const Array<real,2>& pressure_tmp,
              const MeshTool::MeshBlock& Grid_tmp  
            );


    Array<real,2> Solve();

private:
    void SolveFluxF();
    void SolveFluxG();
    void SolveDivP();

};


#endif
