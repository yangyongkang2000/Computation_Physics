//
//  odesolve.h
//  ODESolve
//
//  Created by 杨永康 on 2020/11/28.
//

#ifndef odesolve_h
#define odesolve_h
namespace odesolve {
using namespace std;
template<typename vector_type,typename T=double>
inline vector_type operator+(const vector_type &vec,T _val){
    vector_type _result(vec.size());
    for(int i=0;i<vec.size();i++)
    _result[i]=vec[i]+_val;
    return _result;
}
template<typename function_type,typename T=double>
vector<pair<T,T>> runge_kutta(const function_type&func,const pair<T,T> &_cond,const T &_end, T h=1e-2)
{
    if(h<=0)
        return {};
    T k1,k2,k3,k4;
    vector<pair<T,T>> result {{_cond.first,_cond.second}};
    h=_cond.first>=_end?-h:h;
    auto _f=(_cond.first>=_end)?[](const T &x,const T&_end){return x>=_end;}:[](const T &x,const T&_end){return x<=_end;};
    for(T x=_cond.first,y=_cond.second;_f(x,_end);x+=h)
     {
         k1=func(x,y);
         k2=func(x+h*0.5,y+0.5*h*k1);
         k3=func(x+0.5*h,y+h*k2*0.5);
         k4=func(x+h,y+h*k3);
         y+=h*(k1+2*k2+2*k3+k4)/6.0;
        result.push_back({x+h,y});
     }
    return result;
}
template<typename function_vector_type,typename vector_type,typename T=double>
vector<pair<T,vector_type>> runge_kutta(const function_vector_type &func_vector, pair<T,vector_type>_cond,const T& _end,const T &h=1e-2)
{
   if(h<=0)
       return {};
    vector<pair<T,vector_type>> result {_cond};
    vector_type mid;
    T k1,k2,k3,k4;
    h=_cond.first>=_end?-h:h;
    auto _f=(_cond.first>=_end)?[](const T &x,const T&_end){return x>=_end;}:[](const T &x,const T&_end){return x<=_end;};
    for(T x=_cond.first;_f(x,_end);x+=h)
     {
         for(int i=0;i<func_vector.size();i++)
         {
             k1=func_vector[i](x,_cond.second);
             k2=func_vector[i](x+h*0.5,_cond.second+0.5*h*k1);
             k3=func_vector[i](x+0.5*h,_cond.second+h*k2*0.5);
             k4=func_vector[i](x+h,_cond.second+h*k3);
             mid.push_back(h*(k1+2*k2+2*k3+k4)/6.0);
         }
         for(int i=0;i<_cond.second.size();i++)
         _cond.second[i]+=mid[i];
         result.push_back({x+h,_cond.second});
         mid.clear();
     }
    return result;
}
template<typename function_type,typename vector_type,typename T=double>
vector<pair<T, vector_type>> desolve(const function_type& func,pair<T,vector_type> _cond,const T& _end,const T &h=1e-2)
{
    vector<function_type> func_vector;
    for(int i=1;i<_cond.second.size();i++)
        func_vector.push_back([=](const T&x,const vector_type &v){return v[i];});
    func_vector.push_back(func);
    return runge_kutta<decltype(func_vector),vector_type,T>(func_vector, _cond, _end,h);
}
}
#endif /* odesolve_h */
