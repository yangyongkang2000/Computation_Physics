Dir=../../..
FunctionFit: FunctionFit.cpp
	clang++ -Ofast -I $(Dir)/LinearSolve/LinearSolve -I $(Dir)/FunctionFit/FunctionFit  FunctionFit.cpp -march=native -o libFunctionFit.dylib -std=c++20 -dynamiclib
clean:
	 rm *.dylib
