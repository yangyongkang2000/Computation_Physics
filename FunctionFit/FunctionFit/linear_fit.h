//
//  linear_fit.h
//  FunctionFit
//
//  Created by 杨永康 on 2020/11/19.
//

#ifndef linear_fit_h
#define linear_fit_h
 namespace FunctionFit {
 using namespace std;
template<typename array_type,typename T=double>
vector<T> linear_fit(array_type _array)
{
    size_t N=_array.size(),M=_array[0].size();
    vector<vector<T>> matrix{N,vector<T>(N,0)};
    vector<T> list(N,0);
    for(int k=0;k<N-1;k++)
    {
        for(int i=0;i<M;i++)
          for(int j=0;j<=N-2;j++)
                matrix[k][j]+=_array[i][k]*_array[i][j];
        for(int i=0;i<M;i++)
        {
            matrix[k][N-1]+=_array[i][k];
            list[k]+=_array[i][k]*_array[i][N-1];
        }
    }
    for(int i=0;i<M;i++)
      for(int j=0;j<=N-2;j++)
            matrix[N-1][j]+=_array[i][j];
    for(int i=0;i<M;i++)
    {
        matrix[N-1][N-1]++;
        list[N-1]+=_array[i][N-1];
    }
    return LinearSolve::LU_LinearSolve(matrix, list);
}
}
#endif /* linear_fit_h */
