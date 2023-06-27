#pragma once
#include"2.h"

class L_function_operator_sum:public L_function_operator{
public:
	vector<L_function> functions;
	vector<double> coefficients;
	virtual double operator()(double x);
	virtual L_function_operator *copy();
	void add(L_function f,double a=1);
};
class L_function_operator_f:public L_function_operator{
public:
	double (*f)(double x,double u,double u_);
	L_function u;
	L_function u_;
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function_operator_f_u:public L_function_operator{
public:
	double (*f)(double x,double u,double u_);
	L_function u;
	L_function u_;
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function_operator_f_u_:public L_function_operator{
public:
	double (*f)(double x,double u,double u_);
	L_function u;
	L_function u_;
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function_operator_polynomial:public L_function_operator{
public:
	vector<double> coefficients;
	L_function_operator_polynomial(vector<double> c);
	virtual double operator()(double x);
	virtual L_function_operator *copy();
	L_function_operator_polynomial derivative();
};
class L_function_operator_mesh_u:public L_function_operator{
public:
	L_interval_binary_tree mesh;
	L_function_operator_mesh_u(L_interval_binary_tree m);
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function_operator_mesh_u_:public L_function_operator{
public:
	L_interval_binary_tree mesh;
	L_function_operator_mesh_u_(L_interval_binary_tree m);
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function_operator_mesh_sigma:public L_function_operator{
public:
	L_interval_binary_tree mesh;
	L_function_operator_mesh_sigma(L_interval_binary_tree m);
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};

class L_nonlinear_solver{
	int K;
	double (*f)(double x,double u,double u_);
	double zl0;
	double zl1;
	double zr0;
	double zr1;
	double Gl;
	double Gr;
	double C;
	double TOL;
public:
	L_function_operator_sum u;
	L_function_operator_sum u_;
	L_function_operator_sum u__;
	//initial 初值函数为满足边界条件的多项式, 向量中的第 i 项为 x^i 项系数
	//方程性质不好收敛不了的话可能会爆内存
	L_nonlinear_solver(L_interval_binary_tree _mesh,int _K,double _C,double _TOL,
					   double (*_f)(double x,double u,double u_),
					   double _zl0,double _zl1,double _zr0,double _zr1,vector<double> initial);
};





