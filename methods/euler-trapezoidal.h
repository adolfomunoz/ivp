#ifndef _IVP_EULERTRAPEZOIDAL_H_
#define _IVP_EULERTRAPEZOIDAL_H_

#include "method.h"

namespace IVP
{

class EulerTrapezoidal : public Method<EulerTrapezoidal>
{
	float tolerance;
public:
	EulerTrapezoidal(float s, float tol) :   Method<EulerTrapezoidal>(s),tolerance(tol)  { }
	EulerTrapezoidal(unsigned int ns = 1, float tol = 1.e-3) : Method<EulerTrapezoidal>(ns),tolerance(tol) { }

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, const real& ht) const
	{
		const unsigned int max_iterations = 10000;
		unsigned int i = 0;
		YType f_tyt = f(t,y_t);
		YType y_ti =  y_t + ht*f_tyt;
		real t_h1 = t+ht;
		YType y_ti1 = y_t + ht*real(0.5)*(f_tyt + f(t_h1,y_ti));

		while ( (real(norm(y_ti1 - y_ti))>real(tolerance))  && (i<max_iterations))
		{
			y_ti = y_ti1; 
			y_ti1 = y_t + ht*real(0.5)*(f_tyt + f(t_h1,y_ti));
			i++;
		}
		y_t = y_ti1;
		return t_h1;
	}
};

};

#endif
