#include<valarray>
#include<algorithm>
#include<cmath>
#include<array>
#include<vector>
#include<iostream>
#include<utility>
#include<functional>
namespace NonLinearSolve
{
using namespace std;
template<typename T>
pair<double, double> bitsetion_solve(const T &equ,pair<double,double>&x,double eps=1e-5)
{
    if(equ(x.first)*equ(x.second)>0||x.second<=x.first)
        return pair<double,double>{};
    int n=1+round(log2(x.second-x.first)-log2(eps));
    double y;
    for(int i=0;i<n;i++)
    {
        y=(x.first+x.second)*0.5;
        if(equ(x.second)*equ(y)<0)
            x.first=y;
        else
            x.second=y;
    }
    return x;
}
}
int main(int argc, const char * argv[]) {
    using namespace std;
    pair<double,double> zeros{1,2};
    NonLinearSolve::bitsetion_solve([](double x){return exp(x)*log(x)-x*x;},zeros);
}
