//
//  main.cpp
//  NonLinearSolve
//
//  Created by 杨永康 on 2020/11/6.
//
#include"include.hpp"
#include"LU_Dolittle.hpp"
#include"newtonpro.hpp"
int main(int argc, const char * argv[]) {
    using namespace std;
    vector<double>zeros {0.1,0.1,0.1};
    const array<function<double(vector<double>&)>,3> equ {
        [](vector<double> &v)
        {
            double c_HSO4=v[0],c_SO4=v[1];
            return c_HSO4+c_SO4-0.1;
        },
        [](vector<double> &v)
        {
            double c_HSO4=v[0],c_SO4=v[1],c_H=v[2];
            return  c_HSO4+2*c_SO4-c_H-1e-14/c_H;
        }, [](vector<double> &v)
        {
            double c_HSO4=v[0],c_SO4=v[1],c_H=v[2];
            return  c_SO4*c_H/c_HSO4-0.012;
        }
        
    };
    NonLinearSolve::newton_solve(equ, zeros);
    cout<<zeros[0]<<","<<zeros[1]<<","<<zeros[2]<<endl;
    cout<<"ph="<<-log10(zeros[2])<<endl;
    
}
