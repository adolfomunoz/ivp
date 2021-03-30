#ifndef _IVP_ADAPTIVE_H_
#define _IVP_ADAPTIVE_H_

#include "method.h"
//#include <iostream>

namespace IVP {

class ErrorEstimator
{
public:
	template<typename T>
	T estimate_error(const T& e1, const T& e2);

private:
	template<typename T>
	auto estimate_error_vector(const T& e1, const T& e2) const -> typename std::remove_reference<decltype(*(e1.begin()))>::type
	{
//		typename std::remove_const<typename std::remove_reference<decltype(*(e1.begin()))>::type>::type sol 
//										= decltype(*(e1.begin()))(0.0);
		auto sol = decltype(*(e1.begin()))(0.0);
		auto i1 = e1.begin();
		auto i2 = e2.begin();
		for(; (i1 != e1.end()) && (i2!=e2.end());i1++, i2++)
		{
			auto err = estimate_error(*i1, *i2);
			if (err>sol) sol = err;
		} 
		return sol;
	}

	template<typename T>
	T estimate_error_scalar(const T& e1, const T& e2) const 
	{
		if (std::max(std::abs(e1),std::abs(e2)) < 1.e-6) return std::abs(e1 - e2);
		else return std::abs(e1-e2)/std::max(std::abs(e1),std::abs(e2));
	}

public:
	template<typename T>
	auto estimate_error(const T& e1, const T& e2) const -> decltype(estimate_error_vector(e1,e2))
	{ return estimate_error_vector(e1,e2);	};

	int estimate_error(int e1, int e2) const { return estimate_error_scalar(e1,e2); }
	float estimate_error(float e1, float e2) const { return estimate_error_scalar(e1,e2); }
	double estimate_error(double e1, double e2) const { return estimate_error_scalar(e1,e2); }
};

class StandardStrategy
{
public:
	template<typename real>
	real new_step(const real& step, const real& error, const real& tolerance) const
	{
		return step*
			((error > tolerance)? real(0.5) : 
			((error>1.e-5)? std::min(tolerance/error,real(2.0)) : real(2.0)));
	}
};

template<typename BaseMethod, typename Estimator = ErrorEstimator, typename AdaptationStrategy = StandardStrategy,
	bool embedded = BaseMethod::is_embedded> 
class Adaptive : public Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy,embedded> >
{
	float _tolerance; float min_step;
	BaseMethod _base_method;
	Estimator estimator;
	AdaptationStrategy adaptation;
public:
	Adaptive(unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step) { }
	Adaptive(const BaseMethod& bm, unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm) { }
	Adaptive(const BaseMethod& bm, const Estimator& _estimator, unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm), estimator(_estimator) { }
	Adaptive(const BaseMethod& bm, const Estimator& _estimator, const AdaptationStrategy& _adaptation, 
		 unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm), estimator(_estimator), adaptation(_adaptation) { }

	const BaseMethod& base_method() const { return _base_method; }
	float tolerance() const { return _tolerance; }
	void set_tolerance(float t) { _tolerance = t; }

	/* \brief Indicates which information is passed between steps (apart from the standard one).
         *
         * We need this in order to get the first data created by the base method we are making adaptive
         */
	template<typename YType, typename Function, typename real>
	auto between_steps_first(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const 
		-> decltype(base_method().between_steps_first(f,t_ini,y_ini,t_end))
	{  return base_method().between_steps_first(f,t_ini,y_ini,t_end); }		

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, real& ht) const
	{
		if (ht<=min_step) return base_method().next(f,t,y_t,ht);
		else
		{
			YType s1 = y_t;
			YType s2 = y_t;
			real full_step = ht;
			real half_step = 0.5*ht;
			base_method().next(f,t,s1,full_step);
			real t1 = base_method().next(f,t,s2,half_step);
			real t2 = base_method().next(f,t1,s2,half_step);
			real error = estimator.estimate_error(s1,s2);
			ht = adaptation.new_step(ht,error,real(tolerance()));
			if (error>tolerance()) return next(f,t,y_t,ht);
			else {   y_t = s2; return t2; }
		}
	}

};

template<typename BaseMethod, typename Estimator, typename AdaptationStrategy> 
class Adaptive<BaseMethod, Estimator, AdaptationStrategy, true> :
	 public Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy,true> >
{
	float _tolerance; float min_step;
	BaseMethod _base_method;
	Estimator estimator;
	AdaptationStrategy adaptation;
public:
	Adaptive(unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step) { }
	Adaptive(const BaseMethod& bm, unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm) { }
	Adaptive(const BaseMethod& bm, const Estimator& _estimator, unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm), estimator(_estimator) { }
	Adaptive(const BaseMethod& bm, const Estimator& _estimator, const AdaptationStrategy& _adaptation, 
		 unsigned int ns = 1, float tol = 1.e-5, float _min_step = 0.0) :
		Method<Adaptive<BaseMethod,Estimator,AdaptationStrategy> >(ns), _tolerance(tol), min_step(_min_step), _base_method(bm), estimator(_estimator), adaptation(_adaptation) { }

	const BaseMethod& base_method() const { return _base_method; }
	float tolerance() const { return _tolerance; }
	void set_tolerance(float t) { _tolerance = t; }

	/* \brief Indicates which information is passed between steps (apart from the standard one).
         *
         * We need this in order to get the first data created by the base method we are making adaptive
         */
	template<typename YType, typename Function, typename real>
	typename Type<BaseMethod,YType,Function,real>::BetweenSteps
		 between_steps_first(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const 
	{  return base_method().between_steps_first(f,t_ini,y_ini,t_end); }		

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, real& ht) const
	{
		YType s1 = y_t;
		if (ht<=min_step) return base_method().next_embedded(f,t,y_t,ht,s1);
		else
		{
			YType s2 = y_t;
			real t2 = base_method().next_embedded(f,t,s2,ht,s1);
			real error = estimator.estimate_error(s1,s2);
			ht = adaptation.new_step(ht,error,real(tolerance()));
//			std::cerr<<t<<" - "<<y_t<<" -> "<<s1<<","<<s2<<" - Err = "<<error<<" | Step = "<<ht<<std::endl;
			if (error>tolerance()) return next(f,t,y_t,ht);
			else {   y_t = s2; return t2; }
		}
	}

	template<typename YType, typename Function, typename real, typename BetweenSteps>  
	real next(const Function& f, const real& t, YType& y_t, real& ht, BetweenSteps& bs) const
	{
		YType s1 = y_t;
		if (ht<=min_step) return base_method().next_embedded(f,t,y_t,ht,bs,s1);
		else
		{
			YType s2 = y_t;
			real t2 = base_method().next_embedded(f,t,s2,ht,bs,s1);
			real error = estimator.estimate_error(s1,s2);
			ht = adaptation.new_step(ht,error,real(tolerance()));
//			std::cerr<<t<<" - "<<y_t<<" -> "<<s1<<","<<s2<<" - Err = "<<error<<" | Step = "<<ht<<std::endl;
			if (error>tolerance()) return next(f,t,y_t,ht,bs);
			else {   y_t = s2; return t2; }
		}
	}

};

}; //namespace IVP

#endif
