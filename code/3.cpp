#include"3.h"

double L_function_operator_sum::operator()(double x)
{
	double sum=0;
	for(int i=0;i<functions.size();i++)
		sum+=coefficients[i]*functions[i](x);
	return sum;
}
L_function_operator *L_function_operator_sum::copy()
{
	return new L_function_operator_sum(*this);
}
void L_function_operator_sum::add(L_function f,double a)
{
	functions.push_back(f);
	coefficients.push_back(a);
}

double L_function_operator_f::operator()(double x)
{
	return f(x,u(x),u_(x));
}
L_function_operator *L_function_operator_f::copy()
{
	return new L_function_operator_f(*this);
}
double L_function_operator_f_u::operator()(double x)
{
	double U=u(x);
	double U_=u_(x);
	double h=(abs(x)+abs(f(x,U,U_)))*1e-8;
	auto F=[this,x,U,U_,h](double d)
	{
		return f(x,U+d*h,U_);
	};
	return (F(-2)-8*F(-1)+8*F(1)-F(2))/(12*h);
}
L_function_operator *L_function_operator_f_u::copy()
{
	return new L_function_operator_f_u(*this);
}
double L_function_operator_f_u_::operator()(double x)
{
	double U=u(x);
	double U_=u_(x);
	double h=(abs(x)+abs(f(x,U,U_)))*1e-8;
	auto F=[this,x,U,U_,h](double d)
	{
		return f(x,U,U_+d*h);
	};
	return (F(-2)-8*F(-1)+8*F(1)-F(2))/(12*h);
}
L_function_operator *L_function_operator_f_u_::copy()
{
	return new L_function_operator_f_u_(*this);
}

L_function_operator_polynomial::L_function_operator_polynomial(vector<double> c)
{
	coefficients=c;
}
double L_function_operator_polynomial::operator()(double x)
{
	if(coefficients.size()==0)
		return 0;
	double result=0;
	for(int i=coefficients.size()-1;i>0;i--)
	{
		result+=coefficients[i];
		result*=x;
	}
	result+=coefficients[0];
	return result;
}
L_function_operator *L_function_operator_polynomial::copy()
{
	return new L_function_operator_polynomial(*this);
}
L_function_operator_polynomial L_function_operator_polynomial::derivative()
{
	if(coefficients.size()<2)
		return L_function_operator_polynomial(vector<double>());
	vector<double> c(coefficients.size()-1);
	for(int i=1;i<coefficients.size();i++)
		c[i-1]=i*coefficients[i];
	return L_function_operator_polynomial(c);
}

L_function_operator_mesh_u::L_function_operator_mesh_u(L_interval_binary_tree m)
{
	mesh=m;
}
double L_function_operator_mesh_u::operator()(double x)
{
	return mesh(x);
}
L_function_operator *L_function_operator_mesh_u::copy()
{
	return new L_function_operator_mesh_u(*this);
}
L_function_operator_mesh_u_::L_function_operator_mesh_u_(L_interval_binary_tree m)
{
	mesh=m;
}
double L_function_operator_mesh_u_::operator()(double x)
{
	return mesh.du(x);
}
L_function_operator *L_function_operator_mesh_u_::copy()
{
	return new L_function_operator_mesh_u_(*this);
}
L_function_operator_mesh_sigma::L_function_operator_mesh_sigma(L_interval_binary_tree m)
{
	mesh=m;
}
double L_function_operator_mesh_sigma::operator()(double x)
{
	return mesh.sigma(x);
}
L_function_operator *L_function_operator_mesh_sigma::copy()
{
	return new L_function_operator_mesh_sigma(*this);
}

L_nonlinear_solver::L_nonlinear_solver(L_interval_binary_tree _mesh,int _K,double _C,double _TOL,
									   double(*_f)(double x,double u,double u_),
									   double _zl0,double _zl1,double _zr0,double _zr1,vector<double> initial)
{
	K=_K;
	f=_f;
	zl0=_zl0;
	zl1=_zl1;
	zr0=_zr0;
	zr1=_zr1;
	C=_C;
	TOL=_TOL;
	u=L_function_operator_sum();
	u_=L_function_operator_sum();
	u__=L_function_operator_sum();
	L_interval_binary_tree mesh=_mesh;
	double a=mesh.begin()->left;
	double c=mesh.end()->right;
	bool gtype=abs(zl0)<abs(zl1)&&abs(zr0)<abs(zr1);
	L_function_operator_polynomial ui(initial);
	u.add(ui.copy());
	u_.add(ui.derivative().copy());
	u__.add(ui.derivative().derivative().copy());
	for(int III=0;III<100;III++)
	{
		L_function_operator_f FFF;
		L_function_operator_f_u_ PPP;
		L_function_operator_f_u QQQ;
		QQQ.f=PPP.f=FFF.f=f;
		QQQ.u=PPP.u=FFF.u=L_function(u.copy());
		QQQ.u_=PPP.u_=FFF.u_=L_function(u_.copy());
		L_function_operator_sum QQ,PP,FF;
		QQ.add(QQQ.copy(),-1);
		PP.add(PPP.copy(),-1);
		FF.add(FFF.copy(),1);
		FF.add(u__.copy(),-1);
		L_adaptive_solver S(mesh,K,C,TOL*1e-2,PP.copy(),QQ.copy(),FF.copy(),zl0,zl1,zr0,zr1,0,0);
		u.add(new L_function_operator_mesh_u(S.mesh),1);
		u_.add(new L_function_operator_mesh_u_(S.mesh),1);
		u__.add(new L_function_operator_mesh_sigma(S.mesh),1);
		if(gtype)
			u__.add(new L_function_operator_mesh_u(S.mesh),1);
		double normu=0,normdeltau=0;
		vector<vector<double>> X=S.tau(),Y=S.u();
		for(int i=0;i<X.size();i++)
			for(int j=0;j<X[i].size();j++)
			{
				double uuu=abs(u(X[i][j])),ddduuu=abs(Y[i][j]);
				if(uuu>normu)
					normu=uuu;
				if(ddduuu>normdeltau)
					normdeltau=ddduuu;
			}
		cout<<"--------------\n   d | "<<normdeltau<<"\n   u | "<<normu<<"\n d/u | "<<normdeltau/normu<<"\n--------------\n";
		if(normdeltau/normu<TOL)
			break;
	}
}

