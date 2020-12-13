//
//  FunctionFit.cpp
//  FunctionFit
//
//  Created by 杨永康 on 2020/12/11.
//

#include "include.hpp"
#include "LU_Dolittle.hpp"
#include"linear_fit.h"
#include"tri.hpp"
#include "lagrange_fit.hpp"
#include "cubic_spline.hpp"
#include "WolframLibrary.h"
EXTERN_C DLLEXPORT int CubicSpline(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
{
    using namespace std;
    MTensor tensor=MArgument_getMTensor(Args[0]);
    double *v=libData->MTensor_getRealData(tensor);
    const mint *dims=libData->MTensor_getDimensions(tensor);
    vector<pair<double,double>> mid;
    for(int i=0;i<*dims;i++)
        mid.push_back({v[2*i],v[2*i+1]});
    sort(mid.begin(), mid.end(), [](auto &_,auto &__){return _.first<__.first;});
    vector<vector<double>> list(2);
    for(auto& p:mid)
    {
        list[0].push_back(p.first);
        list[1].push_back(p.second);
    }
    auto ipf=FunctionFit::cubic_spline_fit(list);
    auto& result=ipf._result_;
    mint len=static_cast<mint>(result.size());
    MTensor _tensor;
    libData->MTensor_new(MType_Real,2,(array<mint,2>{len,6}).data(),&_tensor);
    double *_v=libData->MTensor_getRealData(_tensor);
    for(int i=0;i<len;i++)
    {
        for(int j=0;j<4;j++)
          _v[6*i+j]=result[i][j];
        _v[6*i+4]=ipf._list_[i];
        _v[6*i+5]=ipf._list_[i+1];
    }
    MArgument_setMTensor(Res, _tensor);
    return LIBRARY_NO_ERROR;
}
EXTERN_C DLLEXPORT int CubicSplineCalc(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
{
    using namespace std;
    MTensor tensor1=MArgument_getMTensor(Args[0]),tensor2=MArgument_getMTensor(Args[1]),new_tensor;
    const mint *dims1=libData->MTensor_getDimensions(tensor1),*dims2=libData->MTensor_getDimensions(tensor2);
    auto ipf=FunctionFit::Interpolation<vector<double>>();
    double *v1=libData->MTensor_getRealData(tensor1),*v2=libData->MTensor_getRealData(tensor2);
    for(int i=0;i<*dims1;i++)
    {
        ipf._list_.push_back(v1[6*i+4]);
        ipf._result_.push_back({v1[6*i],v1[6*i+1],v1[6*i+2],v1[6*i+3]});
    }
    ipf._list_.push_back(v1[6*(*dims1-1)+5]);
    libData->MTensor_new(MType_Real,1,dims2,&new_tensor);
    double *new_v=libData->MTensor_getRealData(new_tensor);
    for(int i=0;i<*dims2;i++)
    new_v[i]=ipf(v2[i]);
    MArgument_setMTensor(Res, new_tensor);
    return LIBRARY_NO_ERROR;
}
EXTERN_C DLLEXPORT int LagrangeFit(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
{
    using namespace std;
    MTensor tensor=MArgument_getMTensor(Args[0]);
    double *v=libData->MTensor_getRealData(tensor);
    const mint *dims=libData->MTensor_getDimensions(tensor);
    vector<vector<double>> list(2);
    for(int i=0;i<*dims;i++)
    {
        list[0].push_back(v[2*i]);
        list[1].push_back(v[2*i+1]);
    }
   auto result= FunctionFit::lagrange_fit(list);
    MTensor _tensor;
    libData->MTensor_new(MType_Real,1,dims,&_tensor);
    double *_v=libData->MTensor_getRealData(_tensor);
    for(int i=0;i<*dims;i++)
      _v[i]=result[i];
    MArgument_setMTensor(Res, _tensor);
    return LIBRARY_NO_ERROR;
}
EXTERN_C DLLEXPORT int LinearFit(WolframLibraryData libData, mint Argc, MArgument *Args, MArgument Res)
{
    using namespace std;
    MTensor tensor=MArgument_getMTensor(Args[0]),new_tensor;
    double *v=libData->MTensor_getRealData(tensor);
    const mint *dims=libData->MTensor_getDimensions(tensor);
    vector<vector<double>> list(dims[0]);
    for(int i=0;i<*dims;i++)
    for(int j=0;j<dims[1];j++)
    {
        list[i].push_back(v[dims[0]*i+j]);
    }
    auto result=FunctionFit::linear_fit(list);
    libData->MTensor_new(MType_Real,1,dims,&new_tensor);
    double *new_v=libData->MTensor_getRealData(new_tensor);
    for(int i=0;i<*dims;i++)
    new_v[i]=result[i];
    MArgument_setMTensor(Res, new_tensor);
    return LIBRARY_NO_ERROR;
}
