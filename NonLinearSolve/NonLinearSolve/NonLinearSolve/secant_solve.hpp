//
//  secant_solve.hpp
//  NonLinearSolve
//
//  Created by 杨永康 on 2020/11/8.
//

#ifndef secant_solve_h
#define secant_solve_h
namespace NonLinearSolve {
using namespace std;
template<typename _equation,typename T=double>
inline T secant_solve(const _equation &equ,pair<T,T> x,T eps=1e-5,int nmax=1000)
{
    if(equ(x.first)*equ(x.second)>0||x.second<=x.first)
        return NAN;
    T y=x.second,_y=x.first;
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
template<typename _equation,typename T=double>
inline T secant_solve(const _equation &equ,T start,T eps=1e-5,int nmax=1000,T delta=1e-1)
{
    pair<T,T> x{start,start+delta};
    T y;
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
