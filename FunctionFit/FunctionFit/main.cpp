//
//  main.cpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/11/7.
//

#include"include.hpp"
#include"LU_Dolittle.hpp"
#include"lagrange_fit.hpp"
#include "tri.hpp"
#include"cubic_spline.hpp"
int main(int argc, const char * argv[]) {
    using namespace std;
    vector<vector<double>>p{{1,2,3},{2,3,4}};
    auto ip=FunctionFit::cubic_spline_fit(p);
    cout<<ip(2)<<endl;
    return 0;
}
