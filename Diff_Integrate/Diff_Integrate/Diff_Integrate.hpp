//
//  Diff_Integrate.hpp
//  Diff_Integrate
//
//  Created by 杨永康 on 2020/11/12.
//

#ifndef Diff_Integrate_h
#define Diff_Integrate_h
namespace Diff_Integrate {
using namespace std;
template<typename T=double>
class  Diff {
public:
    template<unsigned int N,typename function_type>
    T function_diff(const function_type& func,const T &x)
    {
        T  _diff_num_=0,eps=pow(10,-15.0/N),_t=pow(eps, N);
        _init_p<N>();
        for(int i=0;i<=2*N;i++)
        _diff_num_+=func(x+i*eps-N*eps)*(_p.first)[i]/_t;
        return _diff_num_;
    }
    template<unsigned int N,typename function_type>
    T richason_diff(const function_type&func,const T &x)
    {
        T eps=pow(10,-15.0/N);
        _init_p<N>();
        if(N<=12)
        {
            return (calc___<N>(func, {x,4.0*eps})-20*calc___<N>(func,{x,2.0*eps})+64*calc___<N>(func,{x,eps}))/45.;
        }
        return (4*calc___<N>(func, {x,eps})-calc___<N>(func, {x,2.0*eps}))/3.0;
    }
    template<unsigned int N,typename function_type>
    vector<T> function_diff(const function_type& func,const  vector<T>&x_list)
    {
        T  _diff_num_=0,eps=pow(10,-15/N),_t=pow(eps, N);
        _init_p<N>();
        vector<T> list;
        for(auto &x:x_list)
        {
            _diff_num_=0;
            for(int i=0;i<=2*N;i++)
              _diff_num_+=func(x+i*eps-N*eps)*(_p.first)[i]/_t;
            list.push(_diff_num_);
        }
        return _diff_num_;
    }
    template<unsigned int N,typename function_type>
    vector<T> richason_diff(const function_type&func,const vector<T>&x_list)
    {
        vector<T> list;
        T eps=pow(10,-15.0/N);
        _init_p<N>();
        if(N<=12)
        {
            for(auto &x:x_list)
             list.push_back((calc___<N>(func, {x,4.0*eps})-20*calc___<N>(func,{x,2.0*eps})+64*calc___<N>(func,{x,eps}))/45.);
        }
       else
           for(auto &x:x_list)
               list.push_back((4*calc___<N>(func, {x,eps})-calc___<N>(func, {x,2.0*eps}))/3.0);
        return list;
    }
    void print(int n)
    {
        for(int i=0;i<=2*n;cout<<(_p.first)[i++]<<endl);
    }
private:
    template<unsigned int N,typename function_type>
   inline T calc___(const function_type &func,const pair<T,T> &__p)
        {
            T x=__p.first,eps=__p.second,_t=pow(eps,N),_diff_num_=0;
            for(int i=0;i<=2*N;i++)
              _diff_num_+=func(x+i*eps-N*eps)*(_p.first)[i]/_t;
            return _diff_num_;
        }
    template<unsigned int N>
    constexpr void _init_p()
    {
            _init_p<N-1>();
            for(int i=0;i<=(N-1)<<1;i++)
              (_p.second)[i]=(_p.first)[i];
            for(int i=0;i<=1;i++)
              (_p.first)[i]*=-0.5;
            (_p.first)[2*N-1]=(_p.second)[2*N-3]*0.5;
            (_p.first)[2*N]=(_p.second)[2*N-2]*0.5;
            for(int i=2;i<=(N-1)<<1;i++)
              (_p.first)[i]=((_p.second)[i-2]-(_p.second)[i])*0.5;
    }
    template<>
   constexpr void _init_p<0>()
    {
        (_p.first)[0]=1.;
        return ;
    }
    template<>
   constexpr void _init_p<1>()
    {
        (_p.first)[0]=-0.5;
        (_p.first)[1]=0;
        (_p.first)[2]=0.5;
        return ;
    }
    pair<array<T,10000>,array<T,10000>> _p{{},{}};
};
template<typename function_type,typename T=double>
inline T simpson_integrate(const function_type &func,const pair<T, T>& _p,int n=100)
{
    T h=(_p.second-_p.first)/n,_sum=0;
    for(int i=0;i<=n-1;i++)
       _sum+=4*func(_p.first+h*i+0.5*h);
    for(int i=1;i<=n-1;i++)
     _sum+=2*func(_p.first+h*i);
    return (_sum+func(_p.first)+func(_p.second))*h/6.0;
}
template<typename  function_type,typename T=double>
inline vector<T> simpson_integrate(const vector<function_type> &func_vector,const pair<T, T>& _p,int n=100)
{
    T h=(_p.second-_p.first)/n,_sum=0;
    vector<T> result;
    for(auto &func:func_vector)
    {
        _sum=0;
        for(int i=0;i<=n-1;i++)
           _sum+=4*func(_p.first+h*i+0.5*h);
        for(int i=1;i<=n-1;i++)
         _sum+=2*func(_p.first+h*i);
        result.push_back((_sum+func(_p.first)+func(_p.second))*h/6.0);
    }
    return result;
}
template<typename function_type,typename T=double>
inline T romberg_integrate(const function_type &func,const pair<T,T>&_p,T eps=0.00001,int _max=20)
{
    vector<vector<T>> table(_max,vector<T>(_max,0));
    T _h=(_p.second-_p.first);
    table[0][0]=_h*0.5*(func(_p.first)+func(_p.second));
    for(int i=1,n=1;i<_max;i++)
    {
        for(int k=0;k<=n-1;k++)
          table[i][0]+=func(_p.first+k*_h+0.5*_h)*0.5*_h;
        table[i][0]+=0.5*table[i-1][0];
        n*=2;_h/=2.0;
        for(int j=1;j<=i;j++)
            table[i][j]=(pow(4.0,j)*table[i][j-1]-table[i-1][j-1])/(pow(4.0,j)-1);
        if(fabs(table[i][i]-table[i-1][i-1])<eps)
            return table[i][i];
    }
    return table[_max-1][_max-1];
}
template<typename function_type,typename T=double>
vector<T> romberg_integrate(const vector<function_type> &func_vector,const pair<T,T>&_p,T eps=0.00001,int _max=20)
{
    vector<T> result;
    for(auto &func:func_vector)
        result.push_back(romberg_integrate<T>(func,_p,eps,max));
    return result;
}
template<typename function_type,typename T=double>
inline T calc_gauss(const function_type &func,const pair<T,T> _p,const array<array<double,2>,15>&gauss_kronrod)
{
    T _sum=0,a=(_p.second-_p.first)*0.5,b=(_p.second-_p.first)*0.5,c=(_p.second+_p.first)*0.5;
    for(auto &[A_i,x_i]:gauss_kronrod)
        _sum+=A_i*a*func(b*x_i+c);
    return _sum;
}
template<typename function_type,typename T=double>
inline T _calc_gauss(const function_type &func,const pair<T,T> _p,const array<array<double,2>,15>&gauss_kronrod)
{
    T _sum=0,a=(_p.second-_p.first)*0.5,b=(_p.second-_p.first)*0.5,c=(_p.second+_p.first)*0.5;
    for(auto &[A_i,x_i]:gauss_kronrod)
        _sum+=A_i*a*func(tan(b*x_i+c))/pow(cos(b*x_i+c),2);
    return _sum;
}
template<typename function_type,typename T=double>
T rec_calc(const function_type&func,const pair<T,T>&_p, const array<array<double,2>,15>&gauss_kronrod,const T ans,const T eps=1e-5)
{
    T mid=(_p.second+_p.first)*0.5,left_sum=calc_gauss(func, {_p.first,mid}, gauss_kronrod),right_sum=calc_gauss(func, {mid,_p.second}, gauss_kronrod),_sum=left_sum+right_sum;
    if(fabs(_sum-ans)<eps)
        return _sum;
    else
        return rec_calc(func, {_p.first,mid}, gauss_kronrod, left_sum,eps)+rec_calc(func, {mid,_p.second}, gauss_kronrod, right_sum,eps);
}
template<typename function_type,typename T=double>
T _rec_calc(const function_type&func,const pair<T,T>&_p, const array<array<double,2>,15>&gauss_kronrod,const T ans,const T eps=1e-5)
{
    T mid=(_p.second+_p.first)*0.5,left_sum=_calc_gauss(func, {_p.first,mid}, gauss_kronrod),right_sum=_calc_gauss(func, {mid,_p.second}, gauss_kronrod),_sum=left_sum+right_sum;
    if(fabs(_sum-ans)<eps)
        return _sum;
    else
        return _rec_calc(func, {_p.first,mid}, gauss_kronrod, left_sum,eps)+_rec_calc(func, {mid,_p.second}, gauss_kronrod, right_sum,eps);
}
template<typename function_type,typename T=double>
inline T gauss_integrate(const function_type &func,const pair<T,T>&_p,T eps=1e-5)
{
     const array<array<double,2>,15> gauss_kronrod {0.0307532, -0.987993, 0.070366, -0.937273, 0.107159, -0.848207, \
         0.139571, -0.724418, 0.166269, -0.570972, 0.186161, -0.394151, \
         0.198431, -0.201194, 0.202578, 0., 0.198431, 0.201194, 0.186161, \
         0.394151, 0.166269, 0.570972, 0.139571, 0.724418, 0.107159, 0.848207, \
         0.070366, 0.937273, 0.0307532, 0.987993};
    return rec_calc<decltype(func),T>(func, _p, gauss_kronrod, calc_gauss(func,_p, gauss_kronrod),eps);
}
template<typename function_type,typename T=double>
inline T _gauss_integrate(const function_type &func,const pair<T,T>&_p,T eps=1e-5)
{
     const array<array<double,2>,15> gauss_kronrod {0.0307532, -0.987993, 0.070366, -0.937273, 0.107159, -0.848207, \
         0.139571, -0.724418, 0.166269, -0.570972, 0.186161, -0.394151, \
         0.198431, -0.201194, 0.202578, 0., 0.198431, 0.201194, 0.186161, \
         0.394151, 0.166269, 0.570972, 0.139571, 0.724418, 0.107159, 0.848207, \
         0.070366, 0.937273, 0.0307532, 0.987993};
    return _rec_calc<decltype(func),T>(func, _p, gauss_kronrod, _calc_gauss(func,_p, gauss_kronrod),eps);
}
template<typename function_type,typename T=double>
inline T gauss_integrate(const function_type &func,const pair<T,T>&_p,const array<array<double,2>,15>&gauss_kronrod,T eps=1e-5)
{
    return _rec_calc<decltype(func),T>(func, _p, gauss_kronrod, _calc_gauss(func,_p, gauss_kronrod),eps);
}
template<typename function_type,typename T=double>
inline T infinity_integrate(const function_type &func,char _mode='2',T _p=0.0,T eps=1e-5)
{
    switch (_mode) {
        case'l':case '0':case'L':case 0:
            return _gauss_integrate<decltype(func),T>(func, {atan(_p),M_PI*0.5},eps);
        case '1':case'u':case'U':case 1:
            return _gauss_integrate<decltype(func),T>(func,{-M_PI*0.5,atan(_p)},eps);
        default:
            return _gauss_integrate<decltype(func),T>(func, {-M_PI*0.5,M_PI*0.5},eps);;
    }
}
template<typename function_type,typename T=double>
inline vector<T> gauss_integrate(const vector<function_type> &func_vector, const pair<T, T>_p,T _eps=1e-5)
    {
        const array<array<double,2>,15> gauss_kronrod {0.0307532, -0.987993, 0.070366, -0.937273, 0.107159, -0.848207, \
            0.139571, -0.724418, 0.166269, -0.570972, 0.186161, -0.394151, \
            0.198431, -0.201194, 0.202578, 0., 0.198431, 0.201194, 0.186161, \
            0.394151, 0.166269, 0.570972, 0.139571, 0.724418, 0.107159, 0.848207, \
            0.070366, 0.937273, 0.0307532, 0.987993};
        vector<T> result;
        for(auto &func:func_vector)
            result.push_back(gauss_integrate<decltype(func),T>(func, _p, gauss_kronrod,_eps));
        return result;
    }
}
#endif /* Diff_Integrate_h */
