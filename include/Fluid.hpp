/*
 * Fluid.hpp
 *
 *  Created on: 12/04/2013
 *      Author: jrc-ubu
 */

#ifndef FLUID_HPP_
#define FLUID_HPP_

#include<iostream>
#include<blitz/array.h>


template<typename Tprec>
class Fluid{
public:
	Fluid(Array<Tprec,2> &, //1.-Velocity X-Axis
		  Array<Tprec,2> &, //2.-Velocity Y-Axis
		  Array<Tprec,2> &);//3.-Pressure
private:
	Array<Tprec,2> u, v, p;

};

template<typename Tprec>
	Fluid<Tprec>::Fluid(Array<Tprec,2> &A, Array<Tprec,2> &B, Array<Tprec,2> &C){
		u.resize(A.shape());
		v.resize(B.shape());
		p.resize(C.shape());
		u = A;
		v = B;
		p = C;
}

#endif /* FLUID_HPP_ */
