//
//  main.cpp
//  ODESolve
//
//  Created by 杨永康 on 2020/11/28.
//
#include"include.hpp"
#include "odesolve.h"

int main(int argc, const char * argv[]) {
    // insert code here...
    using namespace std;
    int i=1;
    vector<function<double(const double& ,const vector<double>&)>> f;
    f.push_back([](const double &x,const vector<double>& v){return v[0];});
    f.push_back([&i](const double &x,const vector<double>& v){return v[i];});
    auto g=function<double(const double& ,const vector<double>&)>([](const double &x,const vector<double> &v){return x-v[1]*x;});
    auto h=odesolve::desolve(g, pair<double,vector<double>>({0.0,{0.0,2.0}}), 1.0,0.01);
    for(auto &v:h)
    {
        cout<<v.first<<","<<v.second[0]<<endl;
    }
    return 0;
}
