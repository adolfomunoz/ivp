#ifndef _IVP_METHOD_H_
#define _IVP_METHOD_H_

#include <cmath>
#include <functional>

namespace IVP {

template<typename I>
class ConstIteratorFacade
{
public:
	I& operator++() { static_cast<I*>(this)->inc(); return static_cast<I&>(*this); }
	bool operator==(const I& that) const { return static_cast<const I*>(this)->equals(that);  }
	bool operator!=(const I& that) const { return !static_cast<const I*>(this)->equals(that); }
};

template<typename M>
class MethodEmbedded;

/* \brief Unusable struct for incomplete type stuff. Gathering datatypes in a metaprogramming style
 *         
 */
template<typename M, typename YType, typename Function, typename real>
struct Type
{
	const Function& f; const M& m;
	using BetweenSteps = decltype(m.between_steps_first(f,real(),YType(),real()));
	using Me = M;
};


template<typename M, typename YType, typename Function, typename real>
struct Type<MethodEmbedded<M>, YType, Function, real>
{
	const Function& f; const M& m;
	using BetweenSteps = decltype(m.between_steps_first(f,real(),YType(),real()));
	using Me = M;
};

template<typename M, typename YType, typename Function, typename real, typename BS = typename Type<M,YType,Function,real>::BetweenSteps >
struct BetweenStepsCaller
{
	using T = Type<M,YType,Function,real>;

	static BS between_steps_first(const M& m, const Function& f, const real& ini, const YType& y_ini, const real& end)
	{
		return static_cast<const typename T::Me&>(m).between_steps_first(f,ini,y_ini,end);
	}
};

template<typename M, typename YType, typename Function, typename real>
auto wrapped_between_steps_first(const M& m, const Function& f, const real& ini, const YType& y_ini, const real& end) ->
	decltype(BetweenStepsCaller<M,YType,Function,real>::between_steps_first(m,f,ini,y_ini,end))
{	return BetweenStepsCaller<M,YType,Function,real>::between_steps_first(m,f,ini,y_ini,end); }


template<typename M>
class Method {
	float h;
	unsigned int nsteps;
public:
	Method(float s) :   h(s),nsteps(0) { }
	Method(unsigned int ns) : h(-1.0f),nsteps(ns) { }

	static const bool is_embedded = false;
	static const bool data_between_steps = false;

	int expected_steps()                  const { return nsteps>0?nsteps:int(1.0f/h); }
	template<typename real>
	real step(const real& total) const { return nsteps>0?(total/real(nsteps)):h; }


protected:
	/* For error estimation in methods that work with tolerances */
	static int norm(int t) { return std::abs(t); }
	static float norm(float t) { return std::fabs(t); }
	static double norm(double t) { return std::fabs(t); }

	template<typename V>
	static typename V::value_type norm(const V& v) 
	{
		typename V::value_type sol(0.0);
		typename V::const_iterator i;
		for (i = v.begin(); i != v.end(); i++) 
			if (sol < real(*i)) sol=real(*i);
		return sol;
	}

public:
	/* \brief Indicates which information is passed between steps (apart from the standard one).
         *
         * By default this returns void. Default methods, therefore, just use the y and t of the step, and the step size itself. However,
         * if the method is FSAL or multi-step this should be redefined. 
         */
	template<typename YType, typename Function, typename real>
	void between_steps_first(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const
	{  }		

	template<typename YType, typename real>
	class StepData {

		template<typename YType2, typename Function, typename real2, typename BS>
		friend class Steps;

		YType  _y;
		real   _t;
		real   _step;
		real next_t()	    const { return t() + step(); }
	public:
		StepData(const YType& y = YType(), const real& t = real(), const real& step = real()) :
			_y(y), _t(t), _step(step) { }
		const YType& y()    const { return _y; }
		const real& t()     const { return _t; }
		const real& step()  const { return _step; }
		      YType& y()          { return _y; }
		      real& t()           { return _t; }
		      real& step()        { return _step; }
		
		bool is_pre_last(const real& t_end) const 
		{ return (step()>0)?((t()<t_end) && (next_t()>=t_end)):((t()>t_end) && (next_t()<=t_end)); }
		bool is_last(const real& t_end) const 
		{ return (step()>0)?(t()>=t_end):(t()<=t_end); }
	};

	template<typename YType, typename Function, typename real, 
		typename BetweenSteps = typename Type<M,YType,Function,real>::BetweenSteps>
	class Steps {
		friend class const_iterator;
		const M& m;
		Function f;
		YType y_ini; 
		real  t_ini, t_end;
		BetweenSteps bs_ini;
	public:
		Steps(const M& _m, const Function& _f, const real& _t_ini, const YType& _y_ini, const real& _t_end) :
			m(_m), f(_f), y_ini(_y_ini), t_ini(_t_ini), t_end(_t_end), 
//				bs_ini(m.between_steps_first(_f,_t_ini,_y_ini,_t_end)) { }
				bs_ini(wrapped_between_steps_first(m,_f,_t_ini,_y_ini,_t_end)) { }

		class const_iterator : public ConstIteratorFacade<const_iterator>
		{
			friend class Steps<YType,Function, real>;
			StepData<YType,real> step_data;
			BetweenSteps bs;
			const Steps& steps; bool done;
			const_iterator(const real& step, const Steps& _steps) : 
				step_data(_steps.y_ini,_steps.t_ini,step),bs(_steps.bs_ini),steps(_steps),done(false) { }
			const_iterator(const Steps& _steps) : steps(_steps),done(true) { }
		public:
			void inc() 
			{ 
			   if (step_data.is_last(steps.t_end)) done = true;
			   else {
			      if (step_data.is_pre_last(steps.t_end)) step_data.step() = steps.t_end - step_data.t();
                              step_data.t() = steps.m.next(steps.f, step_data.t(), step_data.y(), step_data.step(), bs); 
			   }
			}
			bool equals(const const_iterator& that) const 
			{       return (this->done == that.done); }
			const StepData<YType,real>& operator*() const { return step_data; } 		
		};

		const_iterator begin() const { return const_iterator(m.step(t_end - t_ini), *this);  }
		const_iterator end()   const { return const_iterator(*this); }
	};
public:
	template<typename YType, typename Function, typename real>
	Steps<YType, Function, real> steps(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const
	{
		return Steps<YType, Function, real>(static_cast<const M&>(*this),f,t_ini,y_ini,t_end);
	}

private:
	template<typename Function, typename VChange, typename DVChange>
	class VariableChangedFunction
	{
		Function f; VChange vc; DVChange dvc;
	public:
		VariableChangedFunction(const Function& _f, const VChange& _vc, const DVChange& _dvc) :
			f(_f), vc(_vc), dvc(_dvc) { }

		template<typename real>
		auto operator()(real t, real y) const -> decltype(f(vc(t),y)*dvc(t))
		{	return f(vc(t),y)*dvc(t); }
	};
public:

	template<typename YType, typename Function, typename real, typename VChange, typename InvVChange, typename DVChange>
	Steps<YType, VariableChangedFunction<Function, VChange, DVChange>, real> 
		steps_change_of_variable(const Function& f, real t_ini, const YType& y_ini, real t_end, 
			const VChange& vc, const InvVChange& vc_inv, const DVChange& dvc) const
	{
		return steps(VariableChangedFunction<Function,VChange,DVChange>(f,vc,dvc), vc_inv(t_ini),y_ini, vc_inv(t_end));
	}

	template<typename YType, typename Function, typename real>
	YType solve(const Function& f, real a, const YType& y_a, real b) const
	{
		YType y = y_a;
		for (auto s : steps(f,a,y_a,b)) { y = s.y(); }
		return y;
	}


	template<typename YType, typename Function, typename real, typename VChange, typename InvVChange, typename DInvVChange>
	YType solve_change_of_variable(const Function& f, real a, const YType& y_a, real b,
		const VChange& vc, const InvVChange& vc_inv, const DInvVChange& d_inv_vc) const
	{
		YType y = y_a;
		for (auto s : steps_change_of_variable(f,a,y_a,b,vc,vc_inv, d_inv_vc)) { y = s.y(); }
		return y;
	}

};

template<typename M>
template<typename YType, typename Function, typename real>
class Method<M>::Steps<YType,Function,real,void> {
	friend class const_iterator;
	const M& m;
	Function f;
	YType y_ini; 
	real  t_ini, t_end;
public:
	Steps(const M& _m, const Function& _f, const real& _t_ini, const YType& _y_ini, const real& _t_end) :
		m(_m), f(_f), y_ini(_y_ini), t_ini(_t_ini), t_end(_t_end) { }

	class const_iterator : public ConstIteratorFacade<const_iterator>
	{
		friend class Steps<YType,Function, real>;
		Method<M>::StepData<YType,real> step_data;
		const Steps& steps; bool done;
		const_iterator(const real& step, const Steps& _steps) : 
			step_data(_steps.y_ini,_steps.t_ini,step),steps(_steps),done(false) { }
		const_iterator(const Steps& _steps) : steps(_steps),done(true) { }
	public:
		void inc() 
		{ 
		   if (step_data.is_last(steps.t_end)) done = true;
		   else {
		      if (step_data.is_pre_last(steps.t_end)) step_data.step() = steps.t_end - step_data.t();
                      step_data.t() = steps.m.next(steps.f, step_data.t(), step_data.y(), step_data.step()); 
		   }
		}
		bool equals(const const_iterator& that) const 
		{       return (this->done == that.done); }
		const Method<M>::StepData<YType,real>& operator*() const { return step_data; } 		
	};


	const_iterator begin() const { return const_iterator(m.step(t_end - t_ini), *this);  }
	const_iterator end()   const { return const_iterator(*this); }


};

template<typename M>
class MethodEmbedded : public Method<MethodEmbedded<M>>
{
public:
	static const bool is_embedded = true;

//	The next line is for gcc 4.8
//      using Method<MethodEmbedded<M>>::Method<MethodEmbedded<M>>;
	MethodEmbedded(float s) :   Method<MethodEmbedded<M>>(s)  { }
	MethodEmbedded(unsigned int ns = 1) : Method<MethodEmbedded<M>>(ns) { }

	template<typename YType, typename Function, typename real>  
	real next(const Function& f, const real& t, YType& y_t, real& ht) const
	{
		YType y_t_other;
		return static_cast<const M&>(*this).next_embedded(f,t,y_t,ht,y_t_other);
	}


	template<typename YType, typename Function, typename real, typename BS>  
	real next(const Function& f, const real& t, YType& y_t, real& ht, BS& bs) const
	{
		YType y_t_other;
		return static_cast<const M&>(*this).next_embedded(f,t,y_t,ht,bs,y_t_other);
	}
};

#define IVP_FSAL template<typename YType, typename Function, typename real>\
	YType between_steps_first(const Function& f, const real& t_ini, const YType& y_ini, const real& t_end) const\
	{ return f(t_ini, y_ini); }

}; // namespace IVP

#endif
