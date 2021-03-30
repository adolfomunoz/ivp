#ifndef _IVP_RUNGEKUTTA4_H_
#define _IVP_RUNGEKUTTA4_H_

#include "method.h"
#include "problem.h"
#include <tuple>

namespace IVP
{

/**
 * 
 *
 * Runge Kutta 4 is not FSAL, although we can precompute lots of stuff if we know that the function is linear.
 */
class RungeKutta4 : public Method<RungeKutta4>
{
public:
//	This works only on gcc >= 2.7
//		using Method<RungeKutta4<real>,real>::Method;
	RungeKutta4(float s) :   Method<RungeKutta4>(s)  { }
	RungeKutta4(unsigned int ns) : Method<RungeKutta4>(ns) { }
	RungeKutta4(int ns = 1) : Method<RungeKutta4>((unsigned int)ns) { }

	/* \brief Indicates which information is passed between steps (apart from the standard one).
         *
         * In this case, by default it is void, except when the function is linear, in which case we can take advantage
         * of data calculated at the previous step. 
         */
	template<typename YType, typename Function, typename real>
	void between_steps_first(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const
	{  }		

	template<typename YType, typename F0, typename F1, typename real>  
	auto between_steps_first(const LinearProblem<F0, F1>& f, const real& t_ini, const YType& y_ini, const real& t_end) const
		-> std::tuple<decltype(f.c0(t_ini)),decltype(f.c1(t_ini))>
	{
		return std::make_tuple(f.c0(t_ini),f.c1(t_ini));
	}

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, const real& ht) const
	{
   		YType k1=ht*f(t,y_t);
   		YType k2=ht*f(t+real(0.5)*ht,y_t+real(0.5)*k1);
   		YType k3=ht*f(t+real(0.5)*ht,y_t+real(0.5)*k2);
   		YType k4=ht*f(t+ht,y_t+k3);
		y_t= y_t + real(1.0/6.0)*k1 + real(1.0/3.0)*k2 + real(1.0/3.0)*k3+ real(1.0/6.0)*k4;
		return t+ht;
	}


	template<typename YType, typename F0, typename F1, typename real, typename TupleType>  
	real next(const LinearProblem<F0, F1>& f, const real& t, YType& y_t, const real& ht, TupleType& evals) const
	{
		YType k1 = ht*(std::get<1>(evals)*y_t + std::get<0>(evals));
		evals = std::make_tuple(f.c0(t+real(0.5)*ht),f.c1(t+real(0.5)*ht));
   		YType k2=ht*(std::get<1>(evals)*(y_t+real(0.5)*k1) + std::get<0>(evals));   
   		YType k3=ht*(std::get<1>(evals)*(y_t+real(0.5)*k2) + std::get<0>(evals));   
		evals = std::make_tuple(f.c0(t+ht),f.c1(t+ht));
  		YType k4=ht*(std::get<1>(evals)*(y_t+k3) + std::get<0>(evals));
		y_t = y_t + real(1.0/6.0)*k1 + real(1.0/3.0)*k2 + real(1.0/3.0)*k3+ real(1.0/6.0)*k4;
		return t+ht;
	}
};

};

#endif
