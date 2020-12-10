//
//  main.cpp
//  LinearSolve
//
//  Created by 杨永康 on 2020/11/2.
//

#include "LU_Dolittle.hpp"
#include<iostream>
#include<random>
const int n=3;
using T= double;
int main()
{
    using namespace std;
    array<array<T,n>,n> matrix{1,-2,3,1,-1,2,3,-3,4};
    array<T, n> list{10,-5,3};
    /*random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> elem(-1,1);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
           matrix[i][j]=elem(gen);
        list[i]=elem(gen);
    }*/
    vector< T> result=LinearSolve::LU_LinearSolve<n,T>(matrix,list);
    for(auto &e:result)
      cout<<e<<endl;
}
