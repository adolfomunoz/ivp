#ifndef _IVP_BACKWARDEULER_H_
#define _IVP_BACKWARDEULER_H_

#include "method.h"
#include <cmath>

namespace IVP
{

class BackwardEuler : public Method<BackwardEuler>
{
	float tolerance;
public:
	BackwardEuler(float s, float tol) :   Method<BackwardEuler>(s),tolerance(tol)  { }
	BackwardEuler(unsigned int ns = 1, float tol = 1.e-3) : Method<BackwardEuler>(ns),tolerance(tol) { }

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, const real& ht) const
	{
		const unsigned int max_iterations = 10000;
		unsigned int i = 0;
		YType y_ti = y_t;
		YType y_ti1 = y_t+ht*f(t+ht,y_ti);
		while ((real( norm(y_ti1 - y_ti))>real(tolerance)) && (i<max_iterations))
		{
			y_ti = y_ti1; y_ti1 = y_t+ht*f(t+ht,y_ti); i++;
		}
		y_t = y_ti1;
		return t+ht;
	}
};

};

#endif
