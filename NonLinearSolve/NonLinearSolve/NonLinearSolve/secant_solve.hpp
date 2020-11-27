//
//  secant_solve.hpp
//  NonLinearSolve
//
//  Created by 杨永康 on 2020/11/8.
//

#ifndef secant_solve_h
#define secant_solve_h
#include"bitsection_solve.hpp"
namespace NonLinearSolve {
using namespace std;
template<typename _equation>
inline double secant_solve(const _equation &equ,pair<double,double> x,double eps=1e-5,int nmax=1000)
{
    if(equ(x.first)*equ(x.second)>0||x.second<=x.first)
        return NAN;
    double y=x.second,_y=x.first;
    for(int n=0;fabs(y-_y)>eps&&n<nmax;n++)
    {
        _y=y;
        y=x.second-(x.second-x.first)*equ(x.second)/(equ(x.second)-equ(x.first));
        if(equ(y)*equ(x.second)<0)
            x.first=y;
        else
            x.second=y;
    }
    return (y+_y)*0.5;
}
template<typename _equation>
inline double secant_solve(const _equation &equ,double start,double eps=1e-5,int nmax=1000,double delta=1e-1)
{
    pair<double,double> x{start,start+delta};
    double y;
    for(int n=0;fabs(x.second-x.first)>eps&&n<nmax;n++)
    {
        y=x.second-(x.second-x.first)*equ(x.second)/(equ(x.second)-equ(x.first));
        x.first=x.second;
        x.second=y;
    }
    return (x.first+x.second)*0.5;
}
}

#endif /* secant_solve_h */
