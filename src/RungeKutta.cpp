#include "RungeKutta.h"

Array<real,2> RungeKutta::Solve(
                                  const Array<real,2> A,
                                  const Array<real,2> B,
                                  real constant
                                  )
{
    Array<real,2> Result(A.shape());
    Result= A(Range::all(),Range::all())+ constant*B(Range::all(),Range::all());
    return Result;
}
