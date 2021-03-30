#ifndef _IVP_ADAPTIVERUNGEKUTTA2_H_
#define _IVP_ADAPTIVERUNGEKUTTA2_H_

#include "method.h"
#include <cmath>

namespace IVP
{

class AdaptiveRungeKutta2 : public Method<AdaptiveRungeKutta2>
{
	float tolerance;
public:
	AdaptiveRungeKutta2(float s, float tol) :   Method<AdaptiveRungeKutta2>(s),tolerance(tol)  { }
	AdaptiveRungeKutta2(unsigned int ns = 1, float tol = 1.e-2) : Method<AdaptiveRungeKutta2>(ns),tolerance(tol) { }

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, real& ht) const
	{
  		YType k1=ht*f(t,y_t);
   		YType k2=ht*f(t+real(0.5)*ht,y_t+real(0.5)*k1);

		real d = real(norm(k2)+norm(k1));
		real error = real(norm(k2 - k1))/d;
		if (error > tolerance)
		{
			ht*=real(0.5);
			return next(f,t,y_t,ht);
		}
		else
		{
			real sol = t + ht;
			y_t=y_t+k2;
			if (d>0.00001) ht*=std::min(tolerance/error,real(2.0));
			else ht*=2.0;
			return sol;
		}
	}
};

};

#endif
