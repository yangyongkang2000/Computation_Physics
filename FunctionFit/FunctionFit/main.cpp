//
//  main.cpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/11/7.
//

#include"linear_fit.h"
int main(int argc, const char * argv[]) {
    using namespace std;
    vector<double> v=FunctionFit::linear_fit(vector<vector<double>>
    {{1,2,3},
    {3,4,7},
        {5,6,11},
        
    });
    for(auto &v1:v) cout<<v1<<endl;
    return 0;
}
