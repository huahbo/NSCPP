#include "Cells.h"


void StaggCellX::Solve()
{
    F_IJ(Range::all(), Range(1,Grid.Ny-2)) = FCell.J();
    F_IJ(Range::all(),0) = Boundary(0,Range::all());
    F_IJ(Range::all(),Grid.Ny-1) = Boundary(1,Range::all());
}

void StaggCellX::Dirichlet()
{
    MathInterpolation F_Dirichlet(F_IJ);
    Val = F_Dirichlet.J();
}

void StaggCellY::Solve()
{
    F_IJ( Range(1,Grid.Nx-2),Range::all()) = FCell.I();
    F_IJ(0,Range::all()) = Boundary(0,Range::all());
    F_IJ(Grid.Nx-1,Range::all()) = Boundary(1,Range::all());
}

void StaggCellY::Dirichlet()
{
    MathInterpolation F_Dirichlet(F_IJ);
    Val = F_Dirichlet.I();
}

