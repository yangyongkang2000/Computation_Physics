//
//  tri.hpp
//  LinearSolve
//
//  Created by 杨永康 on 2020/11/7.
//

#ifndef tri_h
#define tri_h
#include"gauss_seidel.hpp"
namespace LinearSolve
{
   using namespace std;
   template<typename matrix_type,typename list_type,typename T=double>
   vector<T> tri_solve(const matrix_type &matrix,const list_type &list)
{
    int N=static_cast<int>(matrix.size());
    valarray<pair<T,T>> ed(N);
    if(matrix[0][0]==0)
        return vector<T>{};
    ed[0].first=matrix[0][1]/matrix[0][0];
   ed[0].second=list[0]/matrix[0][0];
   double t,precision=1e-322;
   for(int i=1;i<N;i++)
    if(fabs(t=matrix[i][1]-matrix[i][0]*ed[i-1].first)<precision)
       return vector<T>{};
     else {
         ed[i].first=matrix[i][2]/t;
         ed[i].second=(list[i]-matrix[i][0]*ed[i-1].second)/t;
     }
     if(fabs(t=matrix[N-1][1]-matrix[N-1][0]*ed[N-2].first)<precision)
        return vector<T>{};
     vector<T> result(N);
     result[N-1]=ed[N-1].second;
     for(int i=N-2;i>=0;i--)
         result[i]=ed[i].second-ed[i].first*result[i+1];
     return result;
}
   template<unsigned int N,typename T=double>
   vector<T> tri_solve(const array<array<T,3>,N> &matrix,const array<T,N>&list)
   {
       array<pair<T,T>,N> ed;
       if(matrix[0][0]==0)
           return vector<T>{};
       ed[0].first=matrix[0][1]/matrix[0][0];
      ed[0].second=list[0]/matrix[0][0];
      double t,precision=1e-322;
      for(int i=1;i<N;i++)
       if(fabs(t=matrix[i][1]-matrix[i][0]*ed[i-1].first)<precision)
          return vector<T>{};
        else {
            ed[i].first=matrix[i][2]/t;
            ed[i].second=(list[i]-matrix[i][0]*ed[i-1].second)/t;
        }
        if(fabs(t=matrix[N-1][1]-matrix[N-1][0]*ed[N-2].first)<precision)
           return vector<T>{};
        vector<T> result(N);
        result[N-1]=ed[N-1].second;
        for(int i=N-2;i>=0;i--)
            result[i]=ed[i].second-ed[i].first*result[i+1];
        return result;
   }
}
#endif /* tri_h */
