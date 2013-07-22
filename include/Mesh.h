#ifndef _MESHTOOLS_H_
#define _MESHTOOLS_H_

#include "presicion.h"

namespace MeshTool{
class MeshBlock
{
public:

    MeshBlock(real lAxisX, real lAxisY, int nodeX, int nodeY);
    void cellDifference();

    //Reference only you can read same name plus _
    const int &Nx, &Ny; 
    const real &Dx, &Dy; 

    int &cNx; 
    int &cNy; 

    const real &inverseDx;
    const real &inverseDy;

private:
    real lx, ly; //longitud X-Axis Y-Axix
    real startPosX, startPosY; //Start position on the mesh


    int Nx_, Ny_; //nodes X and Y direction
    real Dx_, Dy_; //delta bettwen Cells X-Axis Y-Axis

    int cNx_; //cell nodes X direction
    int cNy_; //cell nodes Y direction

    real inverseDx_;//Gets the inverse of the delta X
    real inverseDy_;//Gets the inverse of the delta Y


};
}
#endif

