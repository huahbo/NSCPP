/*
 * thomas_tem.hpp
 *
 *  Created on: 22/03/2013
 *      Author: cipia
 */

#ifndef THOMAS_TEM_HPP_
#define THOMAS_TEM_HPP_

//---------------------------------------------------------------------------//
// NS-Solver                                                                 //
// Numerical Solver to Navier-Stokes Equation for Incompressible Flow        //
// using Compact Schemes and the Poisson Equation to the Pressure Correction //
//---------------------------------------------------------------------------//
// Project built by Jorge Ramírez Cruz(JramirezCr@iingen.unam.mx             //
// Instituto de Ingeniería UNAM                                              //
//---------------------------------------------------------------------------//
// File Description:                                                         //
// thomas_tem.hpp                                                            //
// Tri-Band Matrix Solver                                                    //
//---------------------------------------------------------------------------//
//namespace NS_Solver{


#include<iostream>
#include<blitz/array.h>

using namespace std;
using namespace blitz;

//Solve Tri-Band Matrix with
template<typename Tprec>
void thomas( Array<Tprec,1> &a,  //1-Band
		     Array<Tprec,1> &b,  //2-Main Diagonal
		     Array<Tprec,1> &c,  //3-Band
		     Array<Tprec,1> &r,  //4-Know Matrix
		     Array<Tprec,1> &u ) //5-Result Matrix
   {


	   int size = u.size();
	   Array<Tprec,1> gam(size);
	   Tprec bet;

	   bet = b(0);
	   u(0) = r(0) / bet;

	   for(int i = 1; i < size; i++){
		   gam(i) = c(i-1) / bet;
		   bet = b(i) - a(i)*gam(i);
		   u(i) = (r(i) - a(i)*u(i-1)) / bet;
	   }


	   for(int i = size - 2; i >= 0; i-- ){
		   u(i) -= gam(i+1)*u(i+1);
	   }

   }

// Solve Tri-Band Matrix whit cyclic boundary conditions
template<typename Tprec>
void thomasCyclic(Array<Tprec,1> &a,  //1-Band
		          Array<Tprec,1> &b,  //2-Principal Band
		          Array<Tprec,1> &c,  //3-Band
		          Tprec alpha,        //4-Corner
		          Tprec beta ,        //5-Corner
		          Array<Tprec,1> &r,  //6-Know Matrix
		          Array<Tprec,1> &x ) //7-Results Matrix
{

	int size = r.size();
	Array<Tprec,1>bb(size),u(size),z(size);
	Tprec fact, gamma;


    gamma = -b(0);
    bb(0) = b(0) - gamma;
    bb(size-1) = b(size-1) - alpha*beta / gamma;

    for(int i = 1; i < size -1; i++) bb(i) = b(i);

    thomas(a, bb, c, r, x);

    u(0) = gamma;
    u(size-1) = alpha;


    for(int i=1; i < size -1; i++) u(i) = 0.0;

    thomas(a, bb, c, u, z);

    fact = (x(0) + beta*x(size-1) / gamma) / (1.0 + z(0) + beta*z(size-1) / gamma);

    for(int i=0; i < size; i++) x(i) -= fact*z(i);


}


#endif /* THOMAS_TEM_HPP_ */
