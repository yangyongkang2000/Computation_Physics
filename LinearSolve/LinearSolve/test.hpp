//
//  test.hpp
//  LinearSolve
//
//  Created by 杨永康 on 2020/11/5.
//

#ifndef test_h
#define test_h
namespace LinearSolve {
template<unsigned int n,typename T=double>
vector<T> gaussian_elimination(array<array<T,n>,n>&a,array<T,n>&b)
{
    vector<double>x(n);
    vector<double>mi_k(n);
    double sum;
    for (int k = 0; k < n - 1; k++)
    {
        //求出第i次初等行变换系数
        for (int j = k + 1; j < n; j++)
        {
            mi_k[j] = a[j][k] / a[k][k];
        }
        for (int i = k + 1; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                a[i][j] = a[i][j] - mi_k[i] * a[k][j];
            }
            b[i] = b[i] - mi_k[i] * b[k];
        }
    }    //回代过程
    x[n - 1] = b[n - 1] / a[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        sum = 0;
        for (int j = i + 1; j < n; j++)
        {
            sum = sum + a[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / a[i][i];
    }
    return x;
}
}
#endif /* test_h */
