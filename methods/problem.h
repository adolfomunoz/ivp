#ifndef _IVP_PROBLEM_H_
#define _IVP_PROBLEM_H_

namespace IVP {

	template<typename Function1, typename Function0>
	class LinearProblem
	{
	public:
		Function1 c1; Function0 c0;
		LinearProblem(const Function1& _c1, const Function0& _c0) : c1(_c1), c0(_c0) { }

		template<typename real, typename YType>
		YType operator()(const real& t, const YType& y) const { return c1(t)*y + c0(t); }
	};

	template<typename Function1, typename Function0>
	LinearProblem<Function1,Function0> linear_problem(const Function1& c1, const Function0& c0)
	{	return LinearProblem<Function1,Function0>(c1,c0); }

};

#endif
