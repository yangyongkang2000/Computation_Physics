//
//  bitsection_solve.hpp
//  NonLinearSolve
//
//  Created by 杨永康 on 2020/11/6.
//

#ifndef bitsection_solve_h
#define bitsection_solve_h
namespace NonLinearSolve
{
using namespace std;
template<typename  _equation,typename T=double>
inline T bitsetion_solve(const _equation& equ,pair<T,T>x,T eps=1e-5)
{
    if(equ(x.first)*equ(x.second)>0||x.second<=x.first)
        return NAN;
    int n=1+round(log2(x.second-x.first)-log2(eps));
    T y;
    for(int i=0;i<n;i++)
    {
        y=(x.first+x.second)*0.5;
        if(equ(x.second)*equ(y)<0)
            x.first=y;
        else
            x.second=y;
    }
    return (x.first+x.second)*0.5;
}
}

#endif /* bitsection_solve_h */
