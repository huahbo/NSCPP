#include "Momento.h"

MomentumX::MomentumX(
                    const StaggCellX& U_tmp,
                    const StaggCellY& V_tmp,
                    const Array<real,2>& pressure_tmp,
                    const MeshTool::MeshBlock& Grid_tmp  
                    )   
    : Momentum::Momentum(U_tmp,V_tmp,pressure_tmp,Grid_tmp)
    {   
        Momentum::Divergence.resize(Ucell.Val.shape());
        Momentum::FluxF.resize(Ucell.Val.shape());
        Momentum::FluxG.resize(Ucell.IJ.shape());
        Momentum::FluxFv.resize(Ucell.Val.shape());
        Momentum::FluxGv.resize(Ucell.Val.shape());
        Momentum::DivP.resize(Ucell.Val.shape());
        std::cout << "Axis X" << std::endl; 
    } 


Array<real,2> MomentumX::Solve()
{
    return Divergence;
}

void MomentumX::SolveFluxF()
{
    FluxF = Momentum::Ucell.Val*Momentum::Ucell.Val;
    FluxFv =  
}




MomentumY::MomentumY(
                    const StaggCellX& U_tmp,
                    const StaggCellY& V_tmp,
                    const Array<real,2>& pressure_tmp,
                    const MeshTool::MeshBlock& Grid_tmp  
                    )   
: Momentum::Momentum(U_tmp,V_tmp,pressure_tmp,Grid_tmp)
{   
    Momentum::Divergence.resize(Vcell.Val.shape());
    Momentum::FluxF.resize(Vcell.Val.shape());
    Momentum::FluxG.resize(Vcell.Val.shape());
    Momentum::FluxFv.resize(Vcell.Val.shape());
    Momentum::FluxGv.resize(Vcell.Val.shape());
    Momentum::DivP.resize(Vcell.Val.shape());
    std::cout << "Axis Y" << std::endl; 
}


Array<real,2> MomentumY::Solve()
{
    return Divergence;
}
