//
//  LU_Dolittle.h
//  LinearSolve
//
//  Created by 杨永康 on 2020/11/4.
//

#ifndef LU_Dolittle_h
#define LU_Dolittle_h
#include"include.hpp"
#include<omp.h>
namespace LinearSolve
{
using namespace std;
   template<unsigned int N,typename T=double>
   vector<T> LU_LinearSolve(array<array<T,N>,N>& matrix,array<T,N>& list)
{
    T sum,precision=1e-322;
    for(int i=1;i<N;i++)
      matrix[i][0]/=matrix[0][0];
    #pragma omp parallel num_threads(4)  private(sum)
    for(int r=1;r<N;r++)
    {
        #pragma omp for
        for(int i=r;i<N;i++)
        {
            sum=0;
            for(int k=0;k<r;k++)
              sum+=matrix[r][k]*matrix[k][i];
            matrix[r][i]-=sum;
        }
        #pragma omp for
        for(int i=r+1;i<N;i++)
        {
            sum=0;
            for(int k=0;k<r;k++)
             sum+=matrix[i][k]*matrix[k][r];
            matrix[i][r]=(matrix[i][r]-sum)/matrix[r][r];
        }
    }
    for(int i=0;i<N;i++)
    {
        sum=0;
        for(int k=0;k<i;k++)
         sum+=matrix[i][k]*list[k];
        list[i]-=sum;
    }
    for(int i=0;i<N;i++)
      if(fabs(matrix[i][i])<precision)
          return vector<T>{};
    vector<T> result(N);
    for(int i=N-1;i>=0;i--)
    {
        sum=0;
        for(int j=i+1;j<=N-1;j++)
           sum+=matrix[i][j]*result[j];
        result[i]=(list[i]-sum)/matrix[i][i];
    }
    return result;
}
}
#endif /* LU_Dolittle_h */
