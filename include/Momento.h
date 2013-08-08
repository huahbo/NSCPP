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
: public 
{
public:

    MomentumX(
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
        Divergence.resize(Ucell.Val.shape());
        FluxF.resize(Ucell.Val.shape());
        FluxG.resize(Ucell.Val.shape());
        FluxFv.resize(Ucell.Val.shape());
        FluxGv.resize(Ucell.Val.shape());
        DivP.resize(Ucell.Val.shape());
        std::cout << "Building Momentum Equation Axis X" << std::endl; 
    }

    Array<real,2> Solve();

private:
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


#endif
