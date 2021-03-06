//
//  cubic_spline.hpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/12/7.
//

#ifndef cubic_spline_h
#define cubic_spline_h
namespace FunctionFit {
using namespace std;
template<typename point_type,typename T=double>
class Interpolation{
public:
    Interpolation(){};
    Interpolation(const point_type &list,const vector<array<T,4>>&result):_list_(list),_result_(result){};
    Interpolation( point_type &&list, vector<array<T,4>>&&result):_list_(list),_result_(result){};
     T operator()(const T&x)
    {
        if(_list_.size()==0)
            return 0.0;
        if(x<=_list_[0])
            return coeff_point<decltype(_result_[0]),T>(_result_[0], x);
        if(x>=_list_[_list_.size()-1])
            return coeff_point<decltype(_result_[0]),T>(_result_[_result_.size()-1],x);
        return coeff_point<decltype(_result_[0]),T>(_result_[distance(_list_.begin(), lower_bound(_list_.begin(), _list_.end(), x))-1], x);
    }
    T derivative(T &x,bool _b)
    {
        if(_list_.size()==0)
            return 0.0;
        auto _func=_b?[](const array<T,4> &_coeff,const T &x){  return 6*_coeff[0]*x+2*_coeff[1];}:[](const array<T,4> &_coeff,const T &x){return _coeff[0]*pow(x,2.0)*3+2*_coeff[1]*x+_coeff[2];};
        if(x<=_list_[0])
            return _func(_result_[0],x);
        if(x>=_list_[_list_.size()-1])
            return _func(_result_[_result_.size()-1],x);
        return  _func(_result_[distance(_list_.begin(), lower_bound(_list_.begin(), _list_.end(), x))-1], x);
    }
    point_type _list_;
    vector<array<T,4>> _result_;
};
template<typename point_type,typename T=double>
inline auto cubic_spline_fit(const point_type &list) ->Interpolation<decltype(list[0]),T>
{
    if(list[0].size()<=2)
        return Interpolation<decltype(list[0]),T>({},{});
    int N=static_cast<int>(list[0].size());
    vector<array<T,4>> result;
    T alpha=list[0][1]-list[0][0],_gamma;
    vector<array<T,3>> matrix {{2*(list[0][2]-list[0][0]),_gamma=(list[0][2]-list[0][1]),0}};
    vector<T> _list_ {6*((list[1][2]-list[1][1])/_gamma-(list[1][1]-list[1][0])/alpha)};
    for(int i=2;i<N-2;i++)
    {
        matrix.push_back({alpha=(list[0][i]-list[0][i-1]),2*(list[0][i+1]-list[0][i-1]),_gamma=(list[0][i+1]-list[0][i])});
        _list_.push_back(6*((list[1][i+1]-list[1][i])/_gamma-(list[1][i]-list[1][i-1])/alpha));
    }
    vector<T> _result {0};
    for(auto &_:LinearSolve::tri_solve<decltype(matrix),decltype(_list_),T>(matrix, _list_))
        _result.push_back(_);
    _result.push_back(0);
    for(int i=0;i<N-1;i++)
    {
        result.push_back({
            (_result[i+1]-_result[i])/(list[0][i+1]-list[0][i])/6.,
            (list[0][i+1]*_result[i]-list[0][i]*_result[i+1])/(list[0][i+1]-list[0][i])*0.5,
           ((pow(list[0][i],2)-2*pow(list[0][i+1],2)-2*list[0][i+1]*list[0][i])*_result[i]+(2*pow(list[0][i],2)-pow(list[0][i+1],2)+2*list[0][i+1]*list[0][i])*_result[i+1]+6*(list[1][i+1]-list[1][i]))/6.0/(list[0][i+1]-list[0][i]),
            -(list[0][i+1]*list[0][i]*(list[0][i]-2*list[0][i+1])*_result[i]+list[0][i]*list[0][i+1]*(2*list[0][i]-list[0][i+1])*_result[i+1]+6*(list[0][i]*list[1][i+1]-list[0][i+1]*list[1][i]))/6.0/(list[0][i+1]-list[0][i])
        });
    }
    return Interpolation<decltype(list[0]),T>(list[0],result);
}
}

#endif /* cubic_spline_h */
