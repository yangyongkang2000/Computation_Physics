//
//  newtonpro.hpp
//  NonLinearSolve
//
//  Created by 杨永康 on 2020/11/8.
//

#ifndef newtonpro_h
#define newtonpro_h
namespace NonLinearSolve
{
using namespace std;
template<typename zeros_type,typename equation_type,typename T=double>
zeros_type newton_solve(const equation_type &equation,zeros_type &zeros,T w=1e-9,double eps=1e-7,int nmax=1000)
{
    int N=static_cast<int>(equation.size());
    vector<vector<T>> A_matrix(N,vector<T>(N));
    vector<T> delta,B_list(N);
   T sum,L;
    for(int n=0;n<nmax;n++)
    {
        for(int i=0;i<N;B_list[i++]=-equation[i](zeros));
        for(int i=0;i<N;i++)
            for(int j=0;j<N;j++)
            {
                zeros[j]+=w;
                L=equation[i](zeros);
                zeros[j]-=2*w;
                A_matrix[i][j]=(L-equation[i](zeros))/(2*w);
                zeros[j]+=w;
            }
        delta=LinearSolve::LU_LinearSolve<decltype(A_matrix),decltype(B_list),T>(A_matrix, B_list);
        if(delta.size()==0)
            return zeros_type{};
        for(int i=0;i<N;i++)
         zeros[i]+=delta[i];
        sum=0;
        for(auto &element:zeros)
            sum+=fabs(element);
        if(sum<eps)
            break;
    }
    return zeros;
}
}
#endif /* newtonpro_h */
