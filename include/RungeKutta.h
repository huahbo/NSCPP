#ifndef _RUNGE_KUTTA_H_
#define _RUNGE_KUTTA_H_

#include<blitz/array.h>

#include "presicion.h"

using namespace blitz;

class RungeKutta
{

public:

    RungeKutta()
    {
    }

    
    Array<real,2> Solve( const Array<real,2> ,const Array<real,2>, real );


};


#endif
