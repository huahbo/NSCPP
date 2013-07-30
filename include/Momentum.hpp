/*
 * Momentum.hpp
 *
 *  Created on: 12/04/2013
 *      Author: jrc-ubu
 */

#ifndef MOMENTUM_HPP_
#define MOMENTUM_HPP_

#include<blitz/array.h>

#include"GeneralEquation.hpp"
#include"Interpol.hpp"
#include"Derivatives.hpp"

using namespace blitz;

template<int Axis>
struct typeAxis{
	enum{StaggeredAxis = Axis};
};

template<typename Tprec, int Axis>
class Momentum
:public GeneralEquation<Tprec>
{

private:

	//Specialization Fluxes direccion Axis-X
	inline void FluxF_Ground(typeAxis<1>){
		F.resize(GeneralEquation<Tprec>::u.shape());
		F = u.Cell()*u.Cell();
	}

	inline void FluxG_Ground(typeAxis<1>){

		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1;
		int sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

		G.resize(sizeX,sizeY);
		G = u.IJ()*v.IJ();
	}

	inline void DivFluxF_Ground(typeAxis<1>){

		dF.resize(GeneralEquation<Tprec>::u.shape());
		dF = d.xM(FluxF());
	}

	inline void DivFluxG_Ground(typeAxis<1>){
		int Ny = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;
		int Nx = GeneralEquation<Tprec>::u.ubound(firstDim) + 1;
		Array<Tprec,2> A(Nx, Ny + 2);
		dG.resize(Nx,Ny);
		A = d.yS(FluxG());
		dG =  A(Range::all(),Range(1, Ny + 1));
	}

	inline void FluxFv_Ground(typeAxis<1>)
	{
		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1,
			sizeY = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;

		Fv.resize(sizeX + 1 ,sizeY);
		Fv = d.xS(u.Cell());
	}

	inline void FluxGv_Ground(typeAxis<1>)
	{
		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1,
			sizeY = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;

		Gv.resize(sizeX ,sizeY + 1);
		Gv = d.yS(u.Cell());
	}

    inline void DivFluxFv_Ground(typeAxis<1>)
    {
		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1,
			sizeY = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;

		Array<Tprec,2> A(sizeX + 2, sizeY);
		dFv.resize(sizeX,sizeY);

		A = d.xS(FluxFv());
		dFv = A(Range(1 , sizeX+1), Range::all());
    }

    inline void DivFluxGv_Ground(typeAxis<1>)
    {
		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1,
			sizeY = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;

		Array<Tprec,2> A(sizeX, sizeY + 2);
		dGv.resize(sizeX,sizeY);

		A = d.yS(FluxGv());
		dGv = A(Range::all(),Range(1 , sizeY + 1));
    }

    inline void DivPress_Ground(typeAxis<1>)
    {

		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1,
			sizeY = GeneralEquation<Tprec>::u.ubound(secondDim) + 1;

		dPress.resize(sizeX,sizeY);
		dPress = d.xS(GeneralEquation<Tprec>::p);

    }



	//Specializacion Fluxes direccion Axis-Y
	inline void FluxF_Ground(typeAxis<2>){

		int sizeX = GeneralEquation<Tprec>::u.ubound(firstDim) + 1;
		int sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

		F.resize(sizeX,sizeY);
		F = u.IJ()*v.IJ();
	}

		inline void FluxG_Ground(typeAxis<2>){

			G.resize(GeneralEquation<Tprec>::v.shape());
			G = v.Cell()*v.Cell();

		}

		inline void DivFluxF_Ground(typeAxis<2>){

			int Ny = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;
			int Nx = GeneralEquation<Tprec>::v.ubound(firstDim) + 1;
			Array<Tprec,2> A(Nx + 2, Ny);
			dF.resize(Nx,Ny);
			A = d.xS(FluxF());
			dF =  A(Range(1, Nx + 1),Range::all());

		}

		inline void DivFluxG_Ground(typeAxis<2>)
		{

			dG.resize(GeneralEquation<Tprec>::v.shape());
			dG = d.yM(FluxG());
		}

		inline void FluxFv_Ground(typeAxis<2>)
		{
				int sizeX = GeneralEquation<Tprec>::v.ubound(firstDim) + 1,
					sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

				Fv.resize(sizeX + 1 ,sizeY);
				Fv = d.xS(v.Cell());
		}

		inline void FluxGv_Ground(typeAxis<2>)
		{
			int sizeX = GeneralEquation<Tprec>::v.ubound(firstDim) + 1,
	        sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

			Gv.resize(sizeX ,sizeY + 1);
			Gv = d.yS(v.Cell());
		}

	    inline void DivFluxFv_Ground(typeAxis<2>)
	    {
			int sizeX = GeneralEquation<Tprec>::v.ubound(firstDim) + 1,
				sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

			Array<Tprec,2> A(sizeX + 2, sizeY);
			dFv.resize(sizeX,sizeY);

			A = d.xS(FluxFv());
			dFv = A(Range(1 , sizeX+1), Range::all());
	    }

	    inline void DivFluxGv_Ground(typeAxis<2>)
	    {
			int sizeX = GeneralEquation<Tprec>::v.ubound(firstDim) + 1,
				sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

			Array<Tprec,2> A(sizeX, sizeY + 2);
			dGv.resize(sizeX,sizeY);

			A = d.yS(FluxGv());
			dGv = A(Range::all(),Range(1 , sizeY + 1));
	    }

	    inline void DivPress_Ground(typeAxis<2>)
	     {

	 		int sizeX = GeneralEquation<Tprec>::v.ubound(firstDim) + 1,
	 			sizeY = GeneralEquation<Tprec>::v.ubound(secondDim) + 1;

	 		dPress.resize(sizeX,sizeY);
	 		dPress = d.yS(GeneralEquation<Tprec>::p);

	     }

       inline void DivInitial(typeAxis<1>)
       {
         Divergence_.resize(GeneralEquation<Tprec>::u.shape());
       }


       inline void DivInitial(typeAxis<2>)
       {
         Divergence_.resize(GeneralEquation<Tprec>::v.shape());
       }


public:

	Momentum(Array<Tprec, 2> &,
			 Array<Tprec, 2> &,
			 Array<Tprec, 2> &,
			 Interpolation<Tprec,1> &,
			 Interpolation<Tprec,2> &,
			 Array<Tprec,1> &);

	inline void FluxF_Ground(typeAxis<Axis>());
	inline void FluxG_Ground(typeAxis<Axis>());

	inline void FluxFv_Ground(typeAxis<Axis>());
	inline void FluxGv_Ground(typeAxis<Axis>());

	inline void DivFluxF_Ground(typeAxis<Axis>());
	inline void DivFluxG_Ground(typeAxis<Axis>());

	inline void DivFluxFv_Ground(typeAxis<Axis>());
	inline void DivFluxGv_Ground(typeAxis<Axis>());

	inline void DivPress_Ground(typeAxis<Axis>());

    inline void DivInitial(typeAxis<Axis>());

	inline Array<Tprec,2> FluxF();
	inline Array<Tprec,2> FluxG();
	inline Array<Tprec,2> DivFluxF();
	inline Array<Tprec,2> DivFluxG();
	inline Array<Tprec,2> FluxFv();
	inline Array<Tprec,2> FluxGv();
	inline Array<Tprec,2> DivFluxFv();
	inline Array<Tprec,2> DivFluxGv();
	inline Array<Tprec,2> DivPressure();
    
    inline void Solve();

    const Array<Tprec,2> & Divergence;

private:

	Interpolation<Tprec, 1> &u;
	Interpolation<Tprec, 2> &v;

	Compact<Tprec> d;

	Array<Tprec,2> F, G, Fv, Gv;
	Array<Tprec,2> dF, dG, dFv, dGv, dPress;
	Array<Tprec,2> Divergence_;
    

};

template<typename Tprec, int Axis>
Momentum<Tprec, Axis>::Momentum(Array<Tprec, 2> &A,
								Array<Tprec, 2> &B,
								Array<Tprec, 2> &C,
								Interpolation<Tprec,1> &UT,
								Interpolation<Tprec,2> &VT,
								Array<Tprec,1> &delta)

:GeneralEquation<Tprec>(A,B,C),u(UT),v(VT),d(delta), Divergence(Divergence_)
{
    DivInitial(typeAxis<Axis>());
}

//Compute the multiplication at the Cell Face in position 1
//we do not interpolate because the Velocity U is at the face.
template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::FluxF()
{
    FluxF_Ground(typeAxis<Axis>());
	return F;
}

//Compute the multiplication at the Cell Face in position 1
//we do not interpolate because the Velocity U is at the face.
template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::FluxG()
{
    FluxG_Ground(typeAxis<Axis>());
	return G;
}

template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::DivFluxF()
{
	DivFluxF_Ground(typeAxis<Axis>());
	return dF;
}


template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::DivFluxG()
{
	DivFluxG_Ground(typeAxis<Axis>());
	return dG;
}

template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::FluxFv(){
	FluxFv_Ground(typeAxis<Axis>());
	return Fv;
}

template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::FluxGv(){
	FluxGv_Ground(typeAxis<Axis>());
	return Gv;
}

template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::DivFluxFv()
{
	DivFluxFv_Ground(typeAxis<Axis>());
	return dFv;
}


template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec, Axis>::DivFluxGv()
{
	DivFluxGv_Ground(typeAxis<Axis>());
	return dGv;
}

template<typename Tprec, int Axis>
inline Array<Tprec,2> Momentum<Tprec,Axis>::DivPressure()
{
	DivPress_Ground(typeAxis<Axis>());
	return dPress;
}


template<typename Tprec, int Axis>
inline void Momentum<Tprec,Axis>::Solve()
{
	DivFluxF_Ground(typeAxis<Axis>());
	DivFluxG_Ground(typeAxis<Axis>());
	DivFluxFv_Ground(typeAxis<Axis>());
	DivFluxGv_Ground(typeAxis<Axis>());
	DivPress_Ground(typeAxis<Axis>());

    Divergence_ = dF + dG - dPress - dFv - dGv;

}


#endif /* MOMENTUM_HPP_ */
