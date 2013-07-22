/*
 * GeneralEquation.cpp
 *
 *  Created on: 12/04/2013
 *      Author: jrc-ubu
 */

#ifndef GENERALEQUATION_CPP_
#define GENERALEQUATION_CPP_

#include<iostream>

#include<blitz/array.h>
using namespace blitz;

template<typename Tprec>
class GeneralEquation{
public:

	GeneralEquation(Array<Tprec, 2> &, Array<Tprec,2> &, Array<Tprec,2> &);

private:
	Array<Tprec,2> temporalTerm, convectiveTerm, diffusiveTerm, pressureTerm;

protected:
	Array<Tprec,2> &u,&v,&p;
};


template<typename Tprec>
GeneralEquation<Tprec>::GeneralEquation(Array<Tprec, 2> &A, Array<Tprec,2> &B, Array<Tprec,2> &C)
	: u(A), v(B), p(C)
{
}


#endif /* GENERALEQUATION_CPP_ */
