#include"LU_Dolittle.hpp"
namespace LinearSolve
{
     template<unsigned int N,typename T=double>
      vector<T> gauss_seidel(array<array<T,N>,N>&matrix,array<T,N>&list,T delta=1e-5,int maxl=1000,double w=1.46)
{
    T precision=1e-322,sum,_z;
    for(int i=0;i<N;i++)
        if(fabs(matrix[i][i])<precision)
        {
            for(int j=0;j<N&&j!=i;j++)
              if(fabs(matrix[j][i])>precision&&fabs(matrix[i][j])>precision)
              {
                  swap(matrix[i],matrix[j]);
                  swap(list[i],list[j]);
                  break;
              }
            if(fabs(matrix[i][i])<precision)
                return vector<T>{};
        }
    vector<T> result(N,0);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
           matrix[i][j]/=i!=j?-matrix[i][i]:1;
        matrix[i][i]=list[i]/matrix[i][i];
    }
    for(int _count=0;_count<maxl;_count++)
    {
        _z=result[0];
        for(int i=0;i<N;i++)
        {
            sum=0;
            for(int j=0;j<N;j++)
              sum+=i!=j?matrix[i][j]*result[j]:matrix[i][i];
            result[i]=sum*w+(1-w)*result[i];
        }
        if(fabs(_z-result[0])<delta)
            break;
    }
    return result;
}
}
