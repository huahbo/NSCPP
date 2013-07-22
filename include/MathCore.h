#ifndef _MATH_CORE_H
#define _MATH_CORE_H

#include "presicion.h"
#include "blitz/array.h"


Array<real,2> AxisXtoCell(const Array<real,2>&, const Array<real,2>&);
Array<real,2> AxisYtoCell(const Array<real,2>&, const Array<real,2>& );

#endif
