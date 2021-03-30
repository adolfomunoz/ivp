#ifndef _IVP_BOGACKI_SHAMPINE_H_
#define _IVP_BOGACKI_SHAMPINE_H_

#include "method.h"
#include "problem.h"

namespace IVP
{

class BogackiShampine : public MethodEmbedded<BogackiShampine>
{
public:
	BogackiShampine(float initial_step) : MethodEmbedded<BogackiShampine>(initial_step) { }
	BogackiShampine(unsigned int ns = 1) : MethodEmbedded<BogackiShampine>(ns) { }
	BogackiShampine(int ns) : MethodEmbedded<BogackiShampine>((unsigned int)ns) { }

	IVP_FSAL;

	template<typename YType, typename Function, typename real>  
	real next_embedded(const Function& f, real t, YType& y_t, real& h, YType& f_t_yt, YType& other) const
	{
		YType k1 = h*f_t_yt;
		YType k2 = h*f(t+(1.0/2.0)*h, y_t + (1.0/2.0)*k1);
		YType k3 = h*f(t+(3.0/4.0)*h, y_t + (3.0/4.0)*k2);
		y_t = y_t + k1*(2.0/9.0) + k2*(1.0/3.0)   + k3*(4.0/9.0);
		f_t_yt = f(t+h         , y_t);
		YType k4 = h*f_t_yt;

		other = y_t + k1*(7.0/24.0)*k1 + (1.0/4.0)*k2 + (1.0/3.0)*k3 + (1.0/8.0)*k4;
		return t+h;
	}
};

};

#endif
