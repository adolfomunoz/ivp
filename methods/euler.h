#ifndef _IVP_EULER_H_
#define _IVP_EULER_H_

#include "method.h"

namespace IVP
{

class Euler : public Method<Euler>
{
public:
//	This should work but doesn't. Tested on gcc 2.7
//	using Method<Euler>::Method;
	Euler(float s) :   Method<Euler>(s)  { }
	Euler(unsigned int ns) : Method<Euler>(ns) { }
	Euler(int ns = 1) : Method<Euler>((unsigned int)ns) { }


	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, const real& ht) const
	{
		y_t=y_t+ht*f(t,y_t);
		return t+ht;
	}
};

};

#endif
