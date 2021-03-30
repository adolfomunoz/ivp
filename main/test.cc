#include <iostream>
#include <iomanip>
#include <ivp.h>
#include <math.h>
#include <functional>
#include <list>
#include <chrono>

struct ErrorTime
{
public:
        float time; //us
        float error;
        ErrorTime(float _t, float _e) : time(_t), error(_e) { }
};

template<typename Function, typename FunctionSol, typename Method>
ErrorTime test(const Method& m, const Function& f, const FunctionSol& fsol, float a, float b)
{

        float fb = m.solve(f,a,fsol(a),b); unsigned int n=1;
        std::chrono::time_point<std::chrono::system_clock> start = std::chrono::system_clock::now();
        while (std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()<=0.01)
        {       fb = m.solve(f,a,fsol(a),b); n++; }
        return ErrorTime(std::chrono::duration<double>(std::chrono::system_clock::now() - start).count()*1.e6/float(n), fabs(fb - fsol(b))/std::max(fb, fsol(b))*100);
}


using namespace std::placeholders;

template<typename Method>
void test_method(const char* id, const Method& m)
{
        std::list<ErrorTime> tests; std::list<ErrorTime> linear_tests;

        auto null_function = [](float t) { return 0.0f; };

        auto f1 = [](float t) { return t; };
        auto f2 = [](float t) { return t*t; };
        auto f3 = [](float t) { return float(log(t)); };
        auto f4 = [](float t) { return float(exp(t)); };
        auto f5 = [](float t) { return float(exp(-2.0f*t)); };
        auto f6 = [](float t) { return float(fabs(t)); };

        auto df1 = [](float t, float y) { return 1.0f; };
        auto df2 = [](float t, float y) { return 2.0f*t; };
        auto df3 = [](float t, float y) { return 1.0f/t; };
        auto df4 = [](float t, float y) { return y; };
        auto df5 = [](float t, float y) { return -2.0f*y; };
        auto df6 = [](float t, float y) { return (t<0)?-1.0f:1.0f; };

        auto lf1 = IVP::linear_problem(null_function, [](float t) { return 1.0f; } );
        auto lf2 = IVP::linear_problem(null_function, [](float t) { return 2.0*t; } );
        auto lf3 = IVP::linear_problem(null_function, [](float t) { return 1.0f/t; } );
        auto lf4 = IVP::linear_problem([](float t) { return 1.0f; }, null_function);
        auto lf5 = IVP::linear_problem([](float t) { return -2.0f; }, null_function);
        auto lf6 = IVP::linear_problem(null_function, [](float t) { return (t<0)?-1.0f:1.0f; } );

        tests.push_back(test(m, df1, f1, 0.0, 1.0));
        tests.push_back(test(m, df2, f2, 0.0, 1.0));
        tests.push_back(test(m, df4, f4, 0.0, 10.0));
        tests.push_back(test(m, df5, f5, 0.0, 10.0));
        tests.push_back(test(m, df3, f3, 0.1, 10.0));
        tests.push_back(test(m, df6, f6, -1, 1));
        linear_tests.push_back(test(m, lf1, f1, 0.0, 1.0));
        linear_tests.push_back(test(m, lf2, f2, 0.0, 1.0));
        linear_tests.push_back(test(m, lf4, f4, 0.0, 10.0));
        linear_tests.push_back(test(m, lf5, f5, 0.0, 10.0));
        linear_tests.push_back(test(m, lf3, f3, 0.1, 10.0));
        linear_tests.push_back(test(m, lf6, f6, -1, 1));

        std::cout<<id<<std::endl;
        std::cout<<"\t y' = 1 \t y' = x \t y' = y \t y'=-2y \t y'=1/x \t y'=-1|1"<<std::endl;
        std::cout<<"\t  [0,1] \t  [0,1] \t [0,10] \t [0,10] \t [.1,10]\t [-1,1] "<<std::endl<<"\t";
        for(auto t: tests) std::cout<<std::fixed<<std::setw(7)<<std::setprecision(2)<<t.error<<"%\t";
        std::cout<<std::endl<<"\t";
        for(auto t: tests) std::cout<<std::fixed<<std::setw(6)<<std::setprecision(2)<<t.time<<"us\t";
        std::cout<<std::endl<<"\t";
        for(auto t: linear_tests) std::cout<<std::fixed<<std::setw(7)<<std::setprecision(2)<<t.error<<"%\t";
        std::cout<<std::endl<<"\t";
        for(auto t: linear_tests) std::cout<<std::fixed<<std::setw(6)<<std::setprecision(2)<<t.time<<"us\t";
        std::cout<<std::endl<<std::endl;
}

int main(int argc, char** argv)
{
        test_method("Euler 0.1 ",IVP::Euler(0.1f));
        test_method("RK2   0.1 ",IVP::RungeKutta2(0.1f));
        test_method("RK4   0.1 ",IVP::RungeKutta4(0.1f));
        test_method("Backward Euler 0.1 1.e-3 ",IVP::BackwardEuler(0.1f,1.e-3f));
        test_method("Trapezoidal 0.1 1.e-3 ",IVP::EulerTrapezoidal(0.1f,1.e-3f));
//        test_method("AdaptiveRK2 0.1 .e-1",IVP::AdaptiveRungeKutta2(0.1f,1.e-1f));
//        test_method("AdaptiveRK2 0.1 .e-2",IVP::AdaptiveRungeKutta2(0.1f,1.e-2f));
        test_method("AdaptiveRK2 1 .e-3",IVP::Adaptive<IVP::EmbeddedRungeKutta2>(1,1.e-3));
        test_method("AdaptiveRK2 0.1 .e-5",IVP::Adaptive<IVP::EmbeddedRungeKutta2>(10,1.e-5));
        test_method("BogackiShampine 1 .e-3",IVP::Adaptive<IVP::BogackiShampine>(1,1.e-3));
        test_method("BogackiShampine 0.1 .e-5",IVP::Adaptive<IVP::BogackiShampine>(10,1.e-5));
        test_method("Dopri 1 .e-3",IVP::Adaptive<IVP::Dopri>(1,1.e-3));
        test_method("Dopri 0.1 .e-5",IVP::Adaptive<IVP::Dopri>(10,1.e-5));
}
