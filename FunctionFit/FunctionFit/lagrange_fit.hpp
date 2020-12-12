//
//  lagrange_fit.hpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/11/7.
//
#ifndef lagrange_fit_h
#define lagrange_fit_h
namespace FunctionFit {
using namespace std;
template<typename point_type,typename T=double>
inline vector<T> lagrange_fit(point_type &list)
{
    int N=static_cast<int>(list[0].size());
    vector<vector<T>> matrix(N,vector<T>(N));
    for(int i=N-1;i>=0;i--)
    {
        matrix[i][N-1]=1;
        for(int j=N-2;j>=0;j--)
        matrix[i][j]=list[0][i]*matrix[i][j+1];
    }
    return LinearSolve::LU_LinearSolve<decltype(matrix),decltype(list[1]),T>(matrix, list[1]);
}
template<typename _coeff_type,typename T=double>
inline T coeff_point(const _coeff_type &_coeff_list,const T &x)
{
    int N=static_cast<int>(_coeff_list.size());
   vector<T> _x(N);
    _x[N-1]=1;
    for(int i=N-2;i>=0;i--)
    _x[i]=_x[i+1]*x;
    return inner_product(_x.begin(), _x.end(), _coeff_list.begin(), static_cast<T>(0));
}
template<typename  _coeff_type,typename T=double>
inline vector<array<T,2>> coeff_list(const _coeff_type &_coeff_list, pair<double,double>_p,double step=0.01)
{
    vector<array<T,2>> v;
    for(;_p.first<=_p.second;_p.first+=step)
    v.push_back({_p.first,coeff_point<decltype(_coeff_list),T>(_coeff_list,_p.first)});
    return v;
}
}

#endif /* lagrange_fit_h */
