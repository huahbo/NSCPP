/*
 * NS_Array.hpp
 *
 *  Created on: 20/03/2013
 *      Author: cipia
 */


namespace NS_Solver
{

#ifndef NS_ARRAY_HPP_
#define NS_ARRAY_HPP_


#include<blitz/array.h>
using namespace blitz;

template< typename Tprec, int Dim>
struct NS_Array
{
	typedef Array<Tprec, Dim> Matrix;
};


#endif /* NS_ARRAY_HPP_ */
}
