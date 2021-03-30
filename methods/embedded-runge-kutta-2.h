#ifndef _IVP_EMBEDDEDRUNGEKUTTA2_H_
#define _IVP_EMBEDDEDRUNGEKUTTA2_H_

#include "method.h"
#include <cmath>

namespace IVP
{

class EmbeddedRungeKutta2 : public MethodEmbedded<EmbeddedRungeKutta2>
{
public:
	EmbeddedRungeKutta2(float s, float tol) :   MethodEmbedded<EmbeddedRungeKutta2>(s) { }
	EmbeddedRungeKutta2(unsigned int ns = 1, float tol = 1.e-2) : MethodEmbedded<EmbeddedRungeKutta2>(ns) { }

	template<typename YType, typename Function, typename real>  
	real next_embedded(const Function& f, const real& t, YType& y_t, real& ht, YType& y_t_other) const
	{
  		YType k1 = ht*f(t,y_t);
   		y_t=y_t+ht*f(t+real(0.5)*ht,y_t+real(0.5)*k1);
		y_t_other = y_t+k1;
		return t+ht;
	}
};

};

#endif
