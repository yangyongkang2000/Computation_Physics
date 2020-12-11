//
//  main.cpp
//  Diff_Integrate
//
//  Created by 杨永康 on 2020/11/12.
//
#include"include.hpp"
#include"AdaptiveIntegration.hpp"

int main(int argc, const char * argv[]) {
    using namespace std;
    
    /*printf("simpson int(tan(x),x,0,1.5):=%.15f\n",Diff_Integrate::simpson_integrate([](double x){return tan(x);}, {0,1.5},1000));
    printf("romberg int(tan(x),x,0,1.5):=%.15f\n",Diff_Integrate::romberg_integrate([](double x){return tan(x);}, {0,1.5},0.000000001));
    return 0;*/
    Diff_Integrate::main();
    /*Diff_Integrate::Diff<long double> _class;
    cout<<_class.richason_diff<10>([](long double x){return sinl(x);}, M_PI/2)<<endl;*/
    
}
