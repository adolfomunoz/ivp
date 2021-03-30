#ifndef _IVP_TEST_TEST_H_
#define _IVP_TEST_TEST_H_

#include <iostream>
#include <ivp/ivp.h>
#include <timer/timer.h>
#include "math.h"

namespace IVP {

std::ostream& operator<<(std::ostream& os, const Euler& method)
{   os<<"Euler-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

std::ostream& operator<<(std::ostream& os, const RungeKutta2& method)
{   os<<"RK2-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

std::ostream& operator<<(std::ostream& os, const RungeKutta4& method)
{   os<<"RK4-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

std::ostream& operator<<(std::ostream& os, const EulerTrapezoidal& method)
{   os<<"EulerTrapezoidal-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

std::ostream& operator<<(std::ostream& os, const BackwardEuler& method)
{   os<<"Backwardeuler-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

std::ostream& operator<<(std::ostream& os, const EmbeddedRungeKutta2& method)
{   os<<"EmbeddedRK2-"<<std::setfill('0')<<std::setw(5)<<method.expected_steps(); return os; }

template<typename T>
std::ostream& operator<<(std::ostream& os, const Adaptive<T>& method)
{   os<<"Adapt-"<<method.base_method()<<"("<<std::setfill('0')<<std::setw(4)<<method.tolerance()<<")"; return os; }

template<typename M>
std::string traits_string(const M& method)
{
    std::stringstream ss;
    ss<<(M::is_embedded?"[E]":"   ");
    return ss.str();
}

template<typename T>
class TestResult
{
	float time;
	unsigned int nevals;
	T sol, approx;
public:
	   TestResult(float _time, unsigned int _nevals, const T& _sol, const T& _approx) :
		time(_time), nevals(_nevals), sol(_sol), approx(_approx) { }

	float absolute_error()     const { return fabs(sol - approx); }
	float relative_error()     const { return fabs(sol - approx)/std::max(fabs(sol),fabs(approx)); }
	float seconds()            const { return time;   } 
	unsigned int evaluations() const { return nevals; }
	T solution()               const { return sol;    }
        T approximation()          const { return approx; }
};

template<typename T>
std::ostream& operator<<(std::ostream& os, const TestResult<T>& tr)
{
   os<<std::fixed<<std::setw(9)<<std::setprecision(3)<<1.e6*tr.seconds()<<"us\t"<<std::setw(6)<<tr.evaluations()<<"ev\t";
   os<<std::fixed<<std::setw(9)<<std::setprecision(3)<<1.e9*tr.seconds()/float(tr.evaluations())<<"ns\t";
   os<<"Error = "<<std::scientific<<std::setprecision(2)<<tr.absolute_error();
   os<<" ("<<std::fixed<<std::setw(6)<<std::setprecision(2)<<100.0*tr.relative_error()<<"%)";
   return os;
}

template<typename F>
class EvaluationCounter
{
	F f;
	static long int n;
public:
	EvaluationCounter(const F& _f) : f(_f) { n=0; }
	float operator()(float t, float y) const { n++; return f(t, y); }
	static long int count() { return n; }
};

template<typename F>
class EvaluationCounterMini
{
	F f;
	static long int n;
public:
	EvaluationCounterMini(const F& _f) : f(_f) { n=0; }
	float operator()(float t) const { n++; return f(t); }
	static long int count() { return n; }
};

template<typename F>
EvaluationCounter<F> evaluation_counter(const F& f) { return EvaluationCounter<F>(f); }

template<typename C1, typename C0>
LinearProblem<C1,EvaluationCounterMini<C0>> evaluation_counter(const LinearProblem<C1,C0>& f) 
{ return LinearProblem<C1,EvaluationCounterMini<C0>>(f.c1,EvaluationCounterMini<C0>(f.c0)); }

template<typename F>
long int count(const EvaluationCounter<F>& f) { return f.count(); }

template<typename C1, typename C0>
long int count(const LinearProblem<C1,EvaluationCounterMini<C0>>& f) 
{ return f.c0.count(); }

template<typename F>
long int EvaluationCounter<F>::n = 0;

template<typename F>
long int EvaluationCounterMini<F>::n = 0;

template<typename Function, typename YType, typename Method, typename A, typename B>
TestResult<YType> test(const Function& f, const YType& sol, const Method& m, const A& a, const YType& y_a, const B& b)
{
	auto fc = evaluation_counter(f);
	unsigned int n;
	YType approx = m.solve(fc,a,y_a,b);
	Timer timer; timer.restart();
	m.solve(f,a,y_a,b); n=1;
	while (timer.runned().seconds()<=0.01)
	{   m.solve(f,a,y_a,b); n++; }
	timer.stop();
	return TestResult<YType>(timer.runned().seconds()/float(n),count(fc),sol,approx);
}

template<typename F>
class EvaluationWriter
{
	F f;
public:
	EvaluationWriter(const F& _f) : f(_f) {  }
	template<typename YType, typename real>
	YType operator()(real t, const YType& y) const 
	{ YType sol = f(t,y); std::cout<<"f("<<t<<","<<y<<")="<<sol<<" "; return sol; }
};


template<typename Function, typename YType, typename Method, typename A, typename B>
void print_steps(const Function& f, const Method& m, const A& a, const YType& y_a, const B& b,const char* prefix = "")
{
	for (auto s: m.steps(f,a,y_a,b) )
		std::cout<<prefix<<s.t()<<" - "<<std::scientific<<s.y()<<" ("<<s.step()<<")"<<std::endl;
}

template<typename Function, typename YType, typename Method, typename A, typename B>
void print_evals(const Function& f, const Method& m, const A& a, const YType& y_a, const B& b,const char* prefix = "")
{
	EvaluationWriter<Function> fw(f);
	auto s = m.steps(fw,a,y_a,b);
	std::cout<<prefix;
	auto i = s.begin(); 
	std::cout<<std::endl;
	while(i!= s.end()) 
	{	std::cout<<prefix<<(*i).t()<<" - "<<std::scientific<<(*i).y()<<" ("<<(*i).step()<<")"<<std::endl;
		std::cout<<prefix;
		++i;
		std::cout<<std::endl;
	}
//	std::cout<<prefix<<(*i).t()<<" - "<<std::scientific<<(*i).y()<<std::endl;
}


}

#endif

