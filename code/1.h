#pragma once
#include"2.h"

class L_direct_solver{
public:
	int K;
	double (*p)(double);
	double (*q)(double);
	double (*f)(double);
	double zl0;
	double zl1;
	double zr0;
	double zr1;
	double Gl;
	double Gr;
	L_interval_binary_tree mesh;
	mat P;
	vec tau;
	vec sigma;
	vec u;

	L_direct_solver(L_interval_binary_tree _mesh,int _K,
					double (*_p)(double),double (*_q)(double),double (*_f)(double),
					double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr);
};