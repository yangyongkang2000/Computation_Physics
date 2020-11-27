#include"Diff_Integrate.hpp"
namespace Diff_Integrate{
class Integrate {
    std::function<double(double)> f;
    double eps;
    // 自适应积分的关键递归部分，基本原理是定积分的区间可加性
    // 首先假设我们有一个可以估算某区间内积分的函数inte和一个初始估计值ans
    // 则将所求区间等分，先用inte分别求出左右区间的估值l_int与r_int，
    // 相加后与ans比较，绝对值小于eps则返回ans
    // 否则，对左右区间递归地应用该函数，采用l_int与r_int作为新的估计值，得到结果后相加
    // 因为划分了两个区间，故递归时要将eps也除以2
    double rec_calc(const double l, const double r, const double ans, const double eps,
                    const std::function<double(double, double)>& inte) {
        double mid = (l + r) / 2;
        double l_int = inte(l, mid),
               r_int = inte(mid, r);
        if (fabs(l_int + r_int - ans) < eps) return ans;
        else return rec_calc(l, mid, l_int, eps/2, inte) 
                  + rec_calc(mid, r, r_int, eps/2, inte);
    }
    // 两种数值积分的方法：辛普森法和高斯法
    std::function<double(double, double)> _simpson;
    static const int n = 7;
    static const double kronrod[n + 1][2];
    std::function<double(double, double)> _gauss;
public:
    explicit Integrate(const std::function<double(double)>& func, const double e = 1e-5)
        :f(func), eps(e) {
        _simpson = [&](const double l, const double r) {
            return (r-l)/6 * ( f(l) + f(r) + f((l+r)/2) * 4 );
        };
        _gauss = [&](const double l, const double r) {
            // 构造g使得int(f, l, r) = int(g, -1, 1)
            auto g = [&](double t) {return (r-l)/2 * f( ( t*(r-l)+l+r )/2 );}; // 线性变换
            // auto g = [&](double t) {
            //     double _l = atan(l), _r = atan(r);
            //     double _t = (t * (_r - _l) + _l + _r) / 2;
            //     return f(tan(_t)) / pow(cos(_t), 2) * (_r - _l) / 2;
            // };                                                                // 利用正切的变换
            double sum = 0.0;
            for (int i = 0; i < n; ++i)
                sum += (g(kronrod[i][0]) + g(-kronrod[i][0])) * kronrod[i][1];
            return sum + g(kronrod[n][0]) * kronrod[n][1];
        };
    }
    void set_eps(const double e) {eps = e;}
    double simpson(const double l, const double r) {
        return rec_calc(l, r, _simpson(l, r), eps, _simpson);
    }
    double gauss(const double l, const double r) {
        return rec_calc(l, r, _gauss(l, r), eps, _gauss);
    }
};
// 15点高斯-克朗罗德法所用到的点及系数
const double Integrate::kronrod[Integrate::n + 1][2] = {
    {0.9914553711, 0.02293532201},
    {0.9491079123, 0.06309209263},
    {0.8648644234, 0.1047900103 },
    {0.7415311856, 0.1406532597 },
    {0.5860872355, 0.1690047266 },
    {0.4058451514, 0.1903505781 },
    {0.2077849550, 0.2044329401 },
    {0.0         , 0.2094821411 }
};
// const double Integrate::kronrod[1 + Integrate::n][2] = {
//     {0.9061798459, 0.2369268851},
//     {0.5384693101, 0.4786286705},
//     {0.0,          0.5688888889}
// };
int main();
}
int Diff_Integrate::main() {
    auto f = [](double x){return exp(-x * x / 2) / sqrt(2 * M_PI);}; 
    // 均值0标准差1的正态分布概率密度函数
    auto g=[](double x){return 1/(1+x*x);};
    Integrate test(f);
    /*const double A = 1;*/
    /*printf("%.15lf\n%.15lf\n%.15lf\n", test.simpson(-A, A), test.gauss(-A, A),romberg_integrate(f,{-A,A},0.000000001));*/
    cout<<Diff_Integrate::infinity_integrate(g)<<endl;
    // 把A改成2500，可以发现高斯法（fx-991CN X等所用积分法）出现了错误
    // 错误的主要来源原因是高斯法需将函数变换至[-1,1]区间计算
    // 而将[-2500,2500]上的正态分布概率密度函数线性变换到[-1,1]上时变得类似单位冲激函数，
    // 导致每次15点取样均取到接近0的值（一般来说，这在自适应算法中不会发生）
    // 导致最后算出接近0的结果
    // 将线性变换改为利用正切进行的变换，可以有效防止此类错误的发生（但是会增大计算量，增大计算误差）
    // 利用正切的变换已写在注释中，取消注释并注释线性变换代码即可使用
    return 0;
}

