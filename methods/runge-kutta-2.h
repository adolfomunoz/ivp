#ifndef _IVP_RUNGEKUTTA2_H_
#define _IVP_RUNGEKUTTA2_H_

#include "method.h"


//This is midpoint's method, maybe we should change its name
namespace IVP
{

class RungeKutta2 : public Method<RungeKutta2>
{
public:
//	This works only on gcc >= 2.7
//		using Method<RungeKutta2<real>,real>::Method;
	RungeKutta2(float s) :   Method<RungeKutta2>(s)  { }
	RungeKutta2(unsigned int ns = 1) : Method<RungeKutta2>(ns) { }
	RungeKutta2(int ns) : Method<RungeKutta2>((unsigned int)ns) { }

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, const real& ht) const
	{
  		YType k1=ht*f(t,y_t);
   		YType k2=ht*f(t+real(0.5)*ht,y_t+real(0.5)*k1);
		y_t=y_t+k2;
		return t+ht;
	}
};

};

#endif
