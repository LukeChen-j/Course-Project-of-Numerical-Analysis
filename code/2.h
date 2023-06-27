#pragma once
#include<iostream>
#include<vector>
#include<cmath>
#include <armadillo>
using namespace std;
using namespace arma;

const double pi=3.14159265358979323;

class L_function_operator{
public:
	double (*func)(double);
	L_function_operator();
	L_function_operator(double (*f)(double));
	virtual double operator()(double x);
	virtual L_function_operator *copy();
};
class L_function{
public:
	L_function_operator *func;
	L_function();
	L_function(double (*f)(double));
	L_function(const L_function &lf);
	L_function(L_function_operator *lco);
	~L_function();
	L_function &operator=(const L_function &lf);
	double operator()(double x) const;
	L_function_operator *operator->();
};

void outputChebyshev(int k);
double Chebyshev(int k,double x);
class L_Chebyshev_polynomial{
public:
	vector<double> f;
	void set_coefficient(vector<double> c);
	void interpolation(vector<double> func);
	double operator()(double x);
};
class L_Chebyshev_integral{
	vector<double> f;
	vector<double> a;
	vector<double> b;
public:
	L_Chebyshev_integral(vector<double> func);
	double operator()(double x);
	double integral();
	double Fl(double x);
	double Fr(double x);
};

class L_interval_data{
public:
	double al;
	double ar;
	double bl;
	double br;
	double dl;
	double dr;
	double m;
	double ml;
	double mr;
	bool abdFlag;
	vec tau;
	vec P_f;
	vec P_psil;
	vec P_psir;
	vec gl;
	vec gr;
	vec dgl;
	vec dgr;
	double s;
	vec sigma;
	double Jl;
	double Jr;
	vec u;
	vec du;
	double S;

	L_interval_data();
};
class L_interval_binary_tree{
	vector<bool> pos;
public:
	L_interval_binary_tree *parent;
	L_interval_binary_tree *child;
	double left;
	double right;
	L_interval_data data;
	L_Chebyshev_polynomial value;
	L_Chebyshev_polynomial dvalue;
	L_Chebyshev_polynomial svalue;

	L_interval_binary_tree();
	L_interval_binary_tree(double l,double r);
	L_interval_binary_tree(const L_interval_binary_tree &l);
	L_interval_binary_tree(const vector<double> &v);
	~L_interval_binary_tree();
	L_interval_binary_tree &operator=(L_interval_binary_tree &l);
	double operator()(double x);
	double du(double x);
	double sigma(double x);

	L_interval_binary_tree *begin();
	L_interval_binary_tree *next();
	L_interval_binary_tree *previous();
	L_interval_binary_tree *end();
	bool isLeft();
	bool isRight();
	bool isRoot();
	bool isLeaf();
	bool add_node();
	bool add_node(double b);
	bool delete_node();

	void cleanData();
	void needUpdate();
	void evaluateABD(int K,L_function p,L_function q,L_function f,double a,double c,
					 double zl0,double zl1,double zr0,double zr1,double uia,double uik,double uib);
	bool evaluateABD();
	bool evaluateMu();
	void evaluateU(int K,double uia,double uik,double uib);
};
class L_solve_on_mesh{
	L_interval_binary_tree *mesh;
	int K;
	L_function p;
	L_function q;
	L_function f;
	double a;
	double c;
	double zl0;
	double zl1;
	double zr0;
	double zr1;
	double Gl;
	double Gr;
	double uia;
	double uik;
	double uib;
public:
	L_solve_on_mesh(L_interval_binary_tree *_mesh,int _K,L_function _p,L_function _q,L_function _f,
					double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr);
	vector<vector<double>> tau();
	vector<vector<double>> u();
	vector<vector<double>> sigma();
};
class L_adaptive_solver{
	int K;
	L_function p;
	L_function q;
	L_function f;
	double zl0;
	double zl1;
	double zr0;
	double zr1;
	double Gl;
	double Gr;
	double C;
	double TOL;
public:
	L_interval_binary_tree mesh;
	L_adaptive_solver(L_interval_binary_tree _mesh,int _K,double _C,double _TOL,
					  L_function _p,L_function _q,L_function _f,
					  double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr);
	vector<vector<double>> tau();
	vector<vector<double>> u();
	vector<vector<double>> sigma();
};




