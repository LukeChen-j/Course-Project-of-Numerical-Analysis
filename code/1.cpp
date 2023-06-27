#include"1.h"

L_direct_solver::L_direct_solver(L_interval_binary_tree _mesh,int _K,double(*_p)(double),double(*_q)(double),double(*_f)(double),double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr)
{
	K=_K;
	p=_p;
	q=_q;
	f=_f;
	zl0=_zl0;
	zl1=_zl1;
	zr0=_zr0;
	zr1=_zr1;
	Gl=_Gl;
	Gr=_Gr;
	mesh=_mesh;
	double a=mesh.begin()->left;
	double c=mesh.end()->right;
	vec::fixed<2> _ui=solve(mat::fixed<2,2>({{zl0*a+zl1,zl0},{zr0*c+zr1,zr0}}),vec::fixed<2>({_Gl,_Gr}));
	double uik=_ui(0);
	double uib=_ui(1);
	bool gtype=abs(zl0)<abs(zl1)&&abs(zr0)<abs(zr1);
	auto gl=[gtype,this,a](double x)
	{
		if(gtype)
			return zl1*cosh(x-a)-zl0*sinh(x-a);
		return zl0*(x-a)-zl1;
	};
	auto gr=[gtype,this,c](double x)
	{
		if(gtype)
			return zr1*cosh(x-c)-zr0*sinh(x-c);
		return zr0*(x-c)-zr1;
	};
	auto dgl=[gtype,this,a](double x)
	{
		if(gtype)
			return zl1*sinh(x-a)-zl0*cosh(x-a);
		return zl0;
	};
	auto dgr=[gtype,this,c](double x)
	{
		if(gtype)
			return zr1*sinh(x-c)-zr0*cosh(x-c);
		return zr0;
	};
	auto p_=[gtype,this](double x)
	{
		return p(x);
	};
	auto q_=[gtype,this](double x)
	{
		if(gtype)
			return q(x)+1;
		return q(x);
	};
	auto f_=[this,uik,uib](double x)
	{
		return f(x)-p(x)*uik-q(x)*(uik*x+uib);
	};
	double s=gl(0)*dgr(0)-dgl(0)*gr(0);
	auto psil=[p_,q_,gl,gr,dgl,dgr,s](double x)
	{
		return (p_(x)*dgr(x)+q_(x)*gr(x))/s;
	};
	auto psir=[p_,q_,gl,gr,dgl,dgr,s](double x)
	{
		return (p_(x)*dgl(x)+q_(x)*gl(x))/s;
	};
	auto G0=[gl,gr,s](double x,double t)
	{
		if(x<t)
			return gl(x)*gr(t)/s;
		else
			return gl(t)*gr(x)/s;
	};
	int M=0;
	for(auto i=mesh.begin();i!=NULL;i=i->next())
		M++;
	P=mat(M*K,M*K,fill::eye);
	tau=vec(M*K);
	vec vf(M*K);
	auto ii=mesh.begin();
	for(int i=0;i<M;i++)
	{
		double center=(ii->right+ii->left)/2,radius=(ii->right-ii->left)/2;
		for(int j=0;j<K;j++)
		{
			double x=center+radius*cos(pi*(K-j-0.5)/K);
			tau(i*K+j)=x;
			vf(i*K+j)=f_(x);
		}
		ii=ii->next();
	}
	ii=mesh.begin();
	for(int i=0;i<M;i++)
	{
		double center=(ii->right+ii->left)/2,radius=(ii->right-ii->left)/2;
		for(int j=0;j<K;j++)
		{
			double x=tau(i*K+j);
			vec el(K,fill::zeros),er(K,fill::zeros);
			el(j)=gl(x)*radius;
			er(j)=gr(x)*radius;
			L_Chebyshev_integral Il(conv_to<vector<double>>::from(el));
			L_Chebyshev_integral Ir(conv_to<vector<double>>::from(er));
			double il=Il.integral();
			double ir=Ir.integral();
			for(int k=0;k<i*K;k++)
				P(k,i*K+j)+=psir(tau(k))*ir;
			for(int k=i*K;k<(i+1)*K;k++)
				P(k,i*K+j)+=psil(tau(k))*Il.Fl((tau(k)-center)/radius)+psir(tau(k))*Ir.Fr((tau(k)-center)/radius);
			for(int k=(i+1)*K;k<M*K;k++)
				P(k,i*K+j)+=psil(tau(k))*il;
		}
		ii=ii->next();
	}
	sigma=solve(P,vf);
//	cout<<"rank: "<<arma::rank(P)<<endl;
	vec vint(K);
	for(int i=0;i<K;i++)
	{
		vec v(K,fill::zeros);
		v(i)=1;
		vint(i)=L_Chebyshev_integral(conv_to<vector<double>>::from(v)).integral();
	}
	vec INT(M*K);
	ii=mesh.begin();
	for(int i=0;i<M;i++)
	{
		INT.subvec(i*K,size(vint))=(ii->right-ii->left)/2*vint;
		ii=ii->next();
	}
	u=vec(M*K);
	for(int i=0;i<M*K;i++)
	{
		vec v(M*K);
		for(int j=0;j<M*K;j++)
			v(j)=G0(tau(i),tau(j))*sigma(j);
		u(i)=dot(v,INT)+uik*(tau(i))+uib;
	}
}
