#ifndef _EXPLICIT_DER_HPP
#define _EXPLICIT_DER_HPP

#include "../Interpol.hpp"
#include <blitz/array.h>
#include "presicion.h"
#include ".h"


class ExplicitDeriv
{
public:
	ExplicitDeriv(const Array<real,2>& vel_u, 
                      const Array<real,2>& vel_v,
                      const Array<real,2>& press_p)
    : u_Velocity(vel_u), v_Velocity(vel_v), Pressure(press_p)
    {}

    Array<real,2> XstaggToCenter(const Array<real,2>& );
    Array<real,2> YstaggToCenter(const Array<real,2>& );
	
private:
	const Array<real,2>& u_Velocity;
	const Array<real,2>& v_Velocity;
	const Array<real,2>& Pressure;    


};

#endif
