#include "Mesh.h"

using namespace MeshTool;

/**
 Constructor 
**/
MeshBlock::MeshBlock(real lAxisX, real lAxisY, int nodeX, int nodeY)
: Nx(Nx_), 
  Ny(Ny_),
  cNx(cNx_), 
  cNy(cNy_), 
  Dx(Dx_), 
  Dy(Dy_),
  inverseDx(inverseDx_), 
  inverseDy(inverseDy_)
{

	lx = (lAxisX > 0.0) ? lAxisX : 1.0 ;
	ly = (lAxisY > 0.0) ? lAxisY : 1.0 ;
    
    Nx_ = (nodeX > 1 ) ? nodeX : 2;
    Ny_ = (nodeY > 1 ) ? nodeY : 2;

    cNx_ = Nx_ - 1;
    cNy_ = Ny_ - 1;

    startPosX = 0.0; //Default begining
    startPosY = 0.0;

}

void MeshBlock::cellDifference()
{
	Dx_ = lx  /  (float)(Nx_ - 1);
	Dy_ = ly /  (float)(Ny_ - 1);	
    inverseDx_ = 1.0 / Dx_;
    inverseDy_ = 1.0 / Dy_;
}



