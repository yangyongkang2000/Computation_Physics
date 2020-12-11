//
//  main.cpp
//  LinearSolve
//
//  Created by 杨永康 on 2020/11/2.
//

#include "include.hpp"
#include "tri.hpp"
#include<ctime>
#include<random>
const int n=100;
using T= double;
int main()
{
    using namespace std;
    array<array<T,3>,n> matrix;
    array<T, n> list;
    for(int i=1;i<=n;i++)
    list[i-1]=i;
    for(int i=1;i<n-1;i++)
    matrix[i]={-1,2,-1};
    matrix[0]={2,-1,0};
    matrix[n-1]={0,-1,2};
    vector<T> result=LinearSolve::tri_solve(matrix, list);
    for(auto &e:result)
      cout<<e<<endl;
}
