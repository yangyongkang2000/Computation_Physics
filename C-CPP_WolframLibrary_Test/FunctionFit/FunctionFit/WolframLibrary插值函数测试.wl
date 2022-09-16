(* ::Package:: *)

CubicSpline=LibraryFunctionLoad["libFunctionFit.dylib","CubicSpline",{{Real,2}},{Real,2}];
CubicSplineFunciton[list_List]:=Function[x,Evaluate@Piecewise[{#1*x^3+#2*x^2+#3*x+#4,#5<=x<#6}&@@@CubicSpline[list]]];
CubicSplineCalc=LibraryFunctionLoad["libFunctionFit.dylib","CubicSplineCalc",{{Real,2},{Real,1}},{Real,1}];
DerivativeCubicSplineCalc=LibraryFunctionLoad["libFunctionFit.dylib","DerivativeCubicSplineCalc",{{Real,2},{Real,1},"Boolean"},{Real,1}];
LagrangeFit=LibraryFunctionLoad["libFunctionFit.dylib","LagrangeFit",{{Real,2}},{Real,1}];
LinearFit=LibraryFunctionLoad["libFunctionFit.dylib","LinearFit",{{Real,2}},{Real,1}];
