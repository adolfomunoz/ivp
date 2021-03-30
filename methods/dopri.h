#ifndef _IVP_DOPRI_H_
#define _IVP_DOPRI_H_

#include "method.h"
#include "problem.h"
#include <cmath>

namespace IVP
{

class DopriErrorStrategy
{
public:
	template<typename real>
	real new_step(const real& step, const real& error, const real& tolerance) const
	{
		return std::max(step/5.0, std::min(10.0*step, step*0.99*pow(tolerance/error,1.0/5.0) ) );
	}
};


class Dopri : public MethodEmbedded<Dopri>
{
public:
	Dopri(float initial_step) : MethodEmbedded<Dopri>(initial_step) { }
	Dopri(unsigned int ns = 1) : MethodEmbedded<Dopri>(ns) { }
	Dopri(int ns) : MethodEmbedded<Dopri>((unsigned int)ns) { }

	IVP_FSAL;

	template<typename YType, typename Function, typename real>  
	real next_embedded(const Function& f, real t, YType& y_t, real& h, YType& f_t_yt, YType& other) const
	{
		YType k0 = h*f_t_yt;
		YType k1 = h*f(t+h/5.0     , y_t + k0/5.0);
		YType k2 = h*f(t+h*3.0/10.0, y_t + k0*3.0/40.0       + k1*9.0/40.0);
		YType k3 = h*f(t+h*4.0/5.0 , y_t + k0*44.0/45.0      - k1*56.0/15.0      + k2*32.0/9.0);
		YType k4 = h*f(t+h*8.0/9.0 , y_t + k0*19372.0/6561.0 - k1*25360.0/2187.0 + k2*64448.0/6561.0 - k3*212.0/729.0);
		YType k5 = h*f(t+h         , y_t + k0*9017.0/3168.0  - k1*355.0/33.0     + k2*46732.0/5247.0 + k3*49.0/176.0 - k4*5103.0/18656.0);
		y_t = y_t + k0*35.0/384.0     + k2*500.0/1113.0   + k3*125.0/192.0   - k4*2187.0/6784.0 + k5*11.0/84.0;
		f_t_yt = f(t+h         , y_t);
		YType k6 = h*f_t_yt;

		other = y_t + (5179.0/57600.0)*k0 + (7571.0/16695.0)*k2 + (393.0/640.0)*k3 - (92097.0/339200.0)*k4
					+ (187.0/2100.0)*k5 + (1.0/40.0)*k6;
		return t+h;
	}

	template<typename YType, typename F0, typename F1, typename real>  
	real next_embedded(const LinearProblem<F0, F1>& f, real t, YType& y_t, real& h, YType& f_t_yt, YType& other) const
	{
		YType k0 = h*f_t_yt;
		YType k1 = h*f(t+h/5.0     , y_t + k0/5.0);
		YType k2 = h*f(t+h*3.0/10.0, y_t + k0*3.0/40.0       + k1*9.0/40.0);
		YType k3 = h*f(t+h*4.0/5.0 , y_t + k0*44.0/45.0      - k1*56.0/15.0      + k2*32.0/9.0);
		YType k4 = h*f(t+h*8.0/9.0 , y_t + k0*19372.0/6561.0 - k1*25360.0/2187.0 + k2*64448.0/6561.0 - k3*212.0/729.0);
		auto c05 = f.c0(t+h); auto c15 = f.c1(t+h);
		YType k5 = h*(c15
			*(y_t + k0*9017.0/3168.0  - k1*355.0/33.0     + k2*46732.0/5247.0 + k3*49.0/176.0 - k4*5103.0/18656.0) + c05);
		y_t = y_t + k0*35.0/384.0     + k2*500.0/1113.0   + k3*125.0/192.0   - k4*2187.0/6784.0 + k5*11.0/84.0;
		f_t_yt = c15*y_t + c05; 
		YType k6 = h*f_t_yt;

		other = y_t + (5179.0/57600.0)*k0 + (7571.0/16695.0)*k2 + (393.0/640.0)*k3 - (92097.0/339200.0)*k4
					+ (187.0/2100.0)*k5 + (1.0/40.0)*k6;
		return t+h;
	}
};


/*
class DopriOld : public Method<DopriOld>
{
	float atol, rtol;
public:
	DopriOld(float _atol, float _rtol, float initial_step):Method<DopriOld>(initial_step),atol(_atol),rtol(_rtol)
	{  }

	DopriOld(float tol = 1.e-6, unsigned int nsteps = 1):Method<DopriOld>(nsteps),atol(tol), rtol(tol)
	{  }

	DopriOld(float _atol, float _rtol, unsigned int nsteps = 1):Method<DopriOld>(nsteps),atol(_atol), rtol(_rtol)
	{  }

private:
	float relative_error_estimate(float err, float y, float next_y) const
	{
		float scale = atol + std::max(fabs(y),fabs(next_y))*rtol;
		return err / scale;	
	}

	double relative_error_estimate(double err, double y, double next_y) const
	{
		double scale = atol + std::max(fabs(y),fabs(next_y))*rtol;
		return err / scale;	
	}

	template<typename V>
	auto relative_error_estimate(const V& err, const V& y, const V& next_y) const -> typename std::remove_reference<decltype(*(y.begin()))>::type
	{
		unsigned int n = 0;
		typename std::remove_const<typename std::remove_reference<decltype(*(y.begin()))>::type>::type sum = 0.0;
		typename V::const_iterator ie, iy, iny;
		for (ie = err.begin(), iy = y.begin(), iny = next_y.begin() ; ie != err.end(); ie++, iy++, iny++)
		{
			auto v = relative_error_estimate(*ie, *iy, *iny);
			sum += v*v;
			n++;
		} 
		return sqrt(sum/(typename std::remove_reference<decltype(*(y.begin()))>::type)(n));
	}

public:
	template<typename YType, typename Function, typename real>  
	YType step_data(const Function& f, real t_ini, const YType& y_ini, real t_end) const
	{
		return f(t_ini, y_ini);
	}	

	template<typename YType, typename F0, typename F1, typename real>  
	real next(const LinearProblem<F0, F1>& f, real t, YType& y_t, real& h, YType& f_t_yt) const
	{
		real sol = t+h;

		//If it is linear, its coefficients are exactly the same for k5 and k6, so we can precalculate them
		YType k0 = h*f_t_yt;
		YType k1 = h*f(t+h/5.0     , y_t + k0/5.0);
		YType k2 = h*f(t+h*3.0/10.0, y_t + k0*3.0/40.0       + k1*9.0/40.0);
		YType k3 = h*f(t+h*4.0/5.0 , y_t + k0*44.0/45.0      - k1*56.0/15.0      + k2*32.0/9.0);
		YType k4 = h*f(t+h*8.0/9.0 , y_t + k0*19372.0/6561.0 - k1*25360.0/2187.0 + k2*64448.0/6561.0 - k3*212.0/729.0);
		auto c05 = f.c0(t+h); auto c15 = f.c1(t+h);
		YType k5 = h*(c15*(y_t + k0*9017.0/3168.0  - k1*355.0/33.0     + k2*46732.0/5247.0 + k3*49.0/176.0 - k4*5103.0/18656.0) + c05);
//		YType k5 = h*f(t+h         , y_t + k0*9017.0/3168.0  - k1*355.0/33.0     + k2*46732.0/5247.0 + k3*49.0/176.0 - k4*5103.0/18656.0);
		YType next_y = y_t + k0*35.0/384.0     + k2*500.0/1113.0   + k3*125.0/192.0   - k4*2187.0/6784.0 + k5*11.0/84.0;
		YType next_f_t_yt = c15*next_y + c05; 
//		YType next_f_t_yt = f(t+h         , next_y);
		YType k6 = h*next_f_t_yt;

		YType err = (35.0/384.0 - 5179.0/57600.0)*k0 
					+ (500.0/1113.0 - 7571.0/16695.0)*k2
					+ (125.0/192.0 - 393.0/640.0)*k3
					+ (-2187.0/6784.0 + 92097.0/339200.0)*k4
					+ (11.0/84.0 - 187.0/2100.0)*k5
					- (1.0/40.0)*k6;


		//Improved error estimate, from Numerical Recipes
		real error=relative_error_estimate(err, y_t, next_y);

		h = std::max(h/5.0, std::min(10.0*h, h*0.99*pow(1.0/error,1.0/5.0) ) );
		if (error >= 1.0) return next(f,t, y_t, h, f_t_yt);
		else
		{
			y_t = next_y;
			f_t_yt = next_f_t_yt; //First equals last, less calculations!
			return sol;
		}
	}

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, real t, YType& y_t, real& h, YType& f_t_yt) const
	{
 		real sol = t+h;

		YType k0 = h*f_t_yt;
		YType k1 = h*f(t+h/5.0     , y_t + k0/5.0);
		YType k2 = h*f(t+h*3.0/10.0, y_t + k0*3.0/40.0       + k1*9.0/40.0);
		YType k3 = h*f(t+h*4.0/5.0 , y_t + k0*44.0/45.0      - k1*56.0/15.0      + k2*32.0/9.0);
		YType k4 = h*f(t+h*8.0/9.0 , y_t + k0*19372.0/6561.0 - k1*25360.0/2187.0 + k2*64448.0/6561.0 - k3*212.0/729.0);
		YType k5 = h*f(t+h         , y_t + k0*9017.0/3168.0  - k1*355.0/33.0     + k2*46732.0/5247.0 + k3*49.0/176.0 - k4*5103.0/18656.0);
		YType next_y = y_t + k0*35.0/384.0     + k2*500.0/1113.0   + k3*125.0/192.0   - k4*2187.0/6784.0 + k5*11.0/84.0;
		YType next_f_t_yt = f(t+h         , next_y);
		YType k6 = h*next_f_t_yt;

		YType err = (35.0/384.0 - 5179.0/57600.0)*k0 
					+ (500.0/1113.0 - 7571.0/16695.0)*k2
					+ (125.0/192.0 - 393.0/640.0)*k3
					+ (-2187.0/6784.0 + 92097.0/339200.0)*k4
					+ (11.0/84.0 - 187.0/2100.0)*k5
					- (1.0/40.0)*k6;


		//Improved error estimate, from Numerical Recipes
		real error = relative_error_estimate(err, y_t, next_y);

		h = std::max(h/5.0, std::min(10.0*h, h*0.99*pow(1.0/error,1.0/5.0) ) );
		if (error >= 1.0) return next(f,t, y_t, h, f_t_yt);
		else
		{
			y_t = next_y;
			f_t_yt = next_f_t_yt; //First equals last, less calculations!
			return sol;
		}
	}
};
*/


};

#endif
