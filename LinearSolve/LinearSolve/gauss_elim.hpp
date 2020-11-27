#include"include.hpp"
namespace  LinearSolve 
{
     using namespace std;
     template<unsigned int N,typename T=double>
   inline  vector<T> gauss_elim(array<valarray<T>,N>&matrix,array<T,N>&list)
{
    
    T precision=1e-322,k,sum,max;
    unsigned int max_position;
    for(int i=0;i<N-1;i++)
    {
            max=fabs(matrix[i][i]);max_position=i;
            for(int j=i+1;j<N;j++)
                 if(fabs(matrix[j][i])>max)
                 {
                     max=fabs(matrix[j][i]);
                     max_position=j;
                 }
            if(fabs(matrix[max_position][i])<precision)
                return vector<T>{};
            if(max_position!=i)
            {
                swap(matrix[i],matrix[max_position]);
                swap(list[i],list[max_position]);
            }
        for(int j=i+1;j<N;j++)
        {
            matrix[j]-=(k=matrix[j][i]/matrix[i][i])*matrix[i];
            list[j]-=k*list[i];
        }
    }
    if(fabs(matrix[N-1][N-1])<precision)
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

