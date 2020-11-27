# Computation_Physics
Numeric Function Programming With C++ 
## 前言
> 虽然名字是计算物理，但主要写的算法是数值计算或者计算方法里面的，当然会增加一些分子动力学或者数学物理方程的程序，但目前我完成了线性方程组，非线性方程组，函数拟合，微分与积分方面的工作，后续会持续更新。

## 简单介绍关于数值计算方面的程序（目前完成的工作）  



### 线性代数
+ Gauss消元解决线性方程组  
```C++
 template<unsigned int N,typename T=double>
   inline  vector<T> gauss_elim(array<valarray<T>,N>&matrix,array<T,N>&list)
```
+ LU分解解决线性方程方程组
```C++
inline vector<T> LU_LinearSolve(matrix_type &matrix,list_type &list)

template<typename matrix_type,typename list_type,typename T=double>
inline vector<T> LU_LinearSolve(matrix_type &matrix,list_type &list)
```
+ 高斯塞德尔迭代法解线性方程组
```C++
template<unsigned int N,typename T=double>
      vector<T> gauss_seidel(array<array<T,N>,N>&matrix,array<T,N>&list,T delta=1e-5,int maxl=1000,double w=1.46)
```
+ 三对角线性方程的解法
```C++
template<unsigned int N,typename T=double>
   vector<T> tri_solve(array<array<T,3>,N> &matrix,array<T,N>&list)
```

### 非线性方程
+ 二分法
```C++
template<typename  _equation>
inline double bitsetion_solve(const _equation& equ,pair<double,double>x,double eps=1e-5)
template<typename _equation>
inline double secant_solve(const _equation &equ,pair<double,double> x,double eps=1e-5,int nmax=1000)
```
+ 弦截法 
``` C++
template<typename _equation>
inline double secant_solve(const _equation &equ,double start,double eps=1e-5,int nmax=1000,double delta=1e-1)
```
+ 牛顿法
```C++
template<unsigned int N,typename zeros_type,typename equation_type>
zeros_type newton_solve(const equation_type &equation,zeros_type &zeros,double w=1e-9,double eps=1e-7,int nmax=1000)
```
### 函数拟合
+ 拉格朗日拟合 
```C++
template<const unsigned int N,typename point_type,typename T=double>
inline vector<T> lagrange_fit(const point_type & list)
```
+ 线性拟合
```C++
template<typename array_type,typename T=double>
vector<T> linear_fit(array_type _array)
```

### 微分与积分
+ 中间差商法求微分
+ 理查森外推法求微分
+ 辛普森法求积分
+ 罗贝格法求积分
+ 高斯法求积分