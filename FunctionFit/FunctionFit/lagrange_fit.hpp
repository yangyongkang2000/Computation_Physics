//
//  lagrange_fit.hpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/11/7.
//

#ifndef lagrange_fit_h
#define lagrange_fit_h

#include"LU_Dolittle.hpp"
namespace FunctionFit {
using namespace std;
template<const unsigned int N,typename point_type,typename T=double>
inline vector<T> lagrange_fit(const point_type & list)
{
    array<array<double,N>,N> matrix{};
    for(int i=N-1;i>=0;i--)
    {
        matrix[i][N-1]=1;
        for(int j=N-2;j>=0;j--)
        matrix[i][j]=list[0][i]*matrix[i][j+1];
    }
    return LinearSolve::LU_LinearSolve<N>(matrix, list);
}
template<const unsigned int N,typename _coeff_type,typename T=double>
inline T coeff_point(const _coeff_type &_coeff_list,const T &x)
{
    array<T,N> _x{1};
    for(int i=1;i<N;i++)
    _x[i]=_x[i-1]*x;
    return inner_product(_x.begin(), _x.end(), _coeff_list.begin(), 0);
    
}
template<const unsigned int N,typename _coeff_type,typename T=double>
inline vector<array<T,2>> coeff_list(const _coeff_type &_coeff_list, pair<double,double>_p,double step=0.01)
{
    vector<array<T,2>> v;
    for(;_p.first<=_p.second;_p.first+=step)
    v.push_back({_p.first,coeff_point<N>(_coeff_list,_p.first)});
    return v;
}
}

#endif /* lagrange_fit_h */
