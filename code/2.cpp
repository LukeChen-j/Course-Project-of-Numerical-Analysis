#include "2.h"

L_function_operator::L_function_operator()
{
	func=[](double x)
	{
		return 0.0;
	};
}
L_function_operator::L_function_operator(double(*f)(double))
{
	func=f;
}
double L_function_operator::operator()(double x)
{
	return func(x);
}
L_function_operator *L_function_operator::copy()
{
	return new L_function_operator(*this);
}
L_function::L_function()
{
	func=NULL;
}
L_function::L_function(double (*f)(double))
{
	func=new L_function_operator(f);
}
L_function::L_function(const L_function &lf)
{
	func=lf.func->copy();
}
L_function::L_function(L_function_operator *lco)
{
	func=lco;
}
L_function::~L_function()
{
	if(func!=NULL)
		delete func;
}
L_function &L_function::operator=(const L_function &lf)
{
	if(func!=NULL)
		delete func;
	func=lf.func->copy();
	return *this;
}
double L_function::operator()(double x) const
{
	if(func==NULL)
		return NAN;
	return (*func)(x);
}
L_function_operator *L_function::operator->()
{
	return func;
}

void outputChebyshev(int k)
{
	vector<vector<int>> c;
	c.push_back(vector<int>{1});
	c.push_back(vector<int>{0,1});
	for(int i=1;i<k;i++)
	{
		vector<int> d;
		d.push_back(-c[i-1][0]);
		for(int j=1;j<i;j++)
			d.push_back(2*c[i][j-1]-c[i-1][j]);
		d.push_back(0);
		d.push_back(2*c[i][i]);
		c.push_back(d);
	}
	for(int i=0;i<k+1;i++)
	{
		cout<<"case "<<i<<":return ";
		if(i%2)
		cout<<";\n";
	}
}
vector<vector<int>> ChebyshevCoefficients{{1},{0,1}};
double Chebyshev(int k,double x)
{
	if(k<0)
		return Chebyshev(-k,x);
	if(k>=ChebyshevCoefficients.size())
		for(int i=ChebyshevCoefficients.size()-1;i<k;i++)
		{
			vector<int> d;
			d.push_back(-ChebyshevCoefficients[i-1][0]);
			for(int j=1;j<i;j++)
				d.push_back(2*ChebyshevCoefficients[i][j-1]-ChebyshevCoefficients[i-1][j]);
			d.push_back(0);
			d.push_back(2*ChebyshevCoefficients[i][i]);
			ChebyshevCoefficients.push_back(d);
		}
	double result=0;
	for(int i=k;i>0;i--)
	{
		result+=ChebyshevCoefficients[k][i];
		result*=x;
	}
	return result+ChebyshevCoefficients[k][0];
}
void L_Chebyshev_polynomial::set_coefficient(vector<double> c)
{
	f=c;
}
void L_Chebyshev_polynomial::interpolation(vector<double> func)
{
	int K=func.size();
	f=vector<double>(K);
	if(K)
	{
		double sum=0;
		for(int j=0;j<K;j++)
			sum+=func[j];
		f[0]=sum/K;
	}
	for(int i=1;i<K;i++)
	{
		double sum=0;
		for(int j=0;j<K;j++)
			sum+=func[j]*cos(i*pi*(K-j-0.5)/K);
		f[i]=sum/K*2;
	}
}
double L_Chebyshev_polynomial::operator()(double x)
{
	double sum=0;
	for(int i=0;i<f.size();i++)
		sum+=f[i]*Chebyshev(i,x);
	return sum;
}
L_Chebyshev_integral::L_Chebyshev_integral(vector<double> func)
{
	int K=func.size();
	f=vector<double>(K);
	if(K)
	{
		double sum=0;
		for(int j=0;j<K;j++)
			sum+=func[j];
		f[0]=sum/K;
	}
	for(int i=1;i<K;i++)
	{
		double sum=0;
		for(int j=0;j<K;j++)
			sum+=func[j]*cos(i*pi*(K-j-0.5)/K);
		f[i]=sum/K*2;
	}
	auto ff=[this,K](int i)
	{
		if(i<K)
			return this->f[i];
		return 0.0;
	};
	a=vector<double>(K+1);
	b=vector<double>(K+1);
	for(int i=2;i<K+1;i++)
	{
		a[i]=(ff(i-1)-ff(i+1))/2/i;
		b[i]=(ff(i+1)-ff(i-1))/2/i;
	}
	if(K)
	{
		a[1]=(2*ff(0)-ff(2))/2;
		b[1]=(ff(2)-2*ff(0))/2;
	}
	{
		double sum=0;
		for(int i=1;i<K+1;i+=2)
			sum+=a[i];
		for(int i=2;i<K+1;i+=2)
			sum-=a[i];
		a[0]=sum;
	}
	{
		double sum=0;
		for(int i=1;i<K+1;i++)
			sum+=b[i];
		b[0]=-sum;
	}
}
double L_Chebyshev_integral::operator()(double x)
{
	double sum=0;
	for(int i=0;i<f.size();i++)
		sum+=f[i]*Chebyshev(i,x);
	return sum;
}
double L_Chebyshev_integral::integral()
{
	double sum=0;
	for(int i=0;i<a.size();i++)
		sum+=a[i];
	return sum;
}
double L_Chebyshev_integral::Fl(double x)
{
	double sum=0;
	for(int i=0;i<a.size();i++)
		sum+=a[i]*Chebyshev(i,x);
	return sum;
}
double L_Chebyshev_integral::Fr(double x)
{
	double sum=0;
	for(int i=0;i<b.size();i++)
		sum+=b[i]*Chebyshev(i,x);
	return sum;
}

L_interval_data::L_interval_data()
{
	al=ar=bl=br=dl=dr=m=ml=mr=Jl=Jr=s=S=0;
	abdFlag=false;
}

L_interval_binary_tree::L_interval_binary_tree()
{
	left=0;
	right=0;
	parent=NULL;
	child=NULL;
}
L_interval_binary_tree::L_interval_binary_tree(double l,double r)
{
	left=l;
	right=r;
	parent=NULL;
	child=NULL;
}
L_interval_binary_tree::L_interval_binary_tree(const L_interval_binary_tree&l)
{
	left=l.left;
	right=l.right;
	pos=l.pos;
	parent=l.parent;
	data=l.data;
	value=l.value;
	dvalue=l.dvalue;
	svalue=l.svalue;
	if(l.child==NULL)
		child=NULL;
	else
	{
		child=new L_interval_binary_tree[2];
		child[0]=l.child[0];
		child[1]=l.child[1];
		child[0].parent=this;
		child[1].parent=this;
	}
}
L_interval_binary_tree::L_interval_binary_tree(const vector<double> &v)
{
	left=0;
	right=0;
	parent=NULL;
	child=NULL;
	int N=v.size();
	if(N>0)
	{
		int n=ceil(log2(--N));
		left=v[0];
		right=v[N];
		for(int i=0;i<n;i++)
		{
			int m=pow(2,n-i);
			auto jj=begin();
			for(int j=m/2;j<N;j+=m)
			{
				jj->add_node(v[j]);
				jj=jj->next();
			}
		}
	}
}
L_interval_binary_tree::~L_interval_binary_tree()
{
	if(child!=NULL)
		delete[] child;
}
L_interval_binary_tree &L_interval_binary_tree::operator=(L_interval_binary_tree &l)
{
	if(child!=NULL)
		delete[] child;
	left=l.left;
	right=l.right;
	pos=l.pos;
	parent=l.parent;
	data=l.data;
	value=l.value;
	dvalue=l.dvalue;
	svalue=l.svalue;
	if(l.child==NULL)
		child=NULL;
	else
	{
		child=new L_interval_binary_tree[2];
		child[0]=l.child[0];
		child[1]=l.child[1];
		child[0].parent=this;
		child[1].parent=this;
	}
	return *this;
}
L_interval_binary_tree *L_interval_binary_tree::begin()
{
	if(child==NULL)
		return this;
	return child->begin();
}
L_interval_binary_tree *L_interval_binary_tree::next()
{
	if(parent==NULL)
		return NULL;
	if(pos[pos.size()-1])
		return parent->next();
	return (this+1)->begin();
}
double L_interval_binary_tree::operator()(double x)
{
	double center=(right+left)/2,radius=(right-left)/2;
	if(isLeaf())
		return value((x-center)/radius);
	if(x<=child[0].right)
		return child[0](x);
	if(x>=child[1].left)
		return child[1](x);
	return NAN;
}
double L_interval_binary_tree::du(double x)
{
	double center=(right+left)/2,radius=(right-left)/2;
	if(isLeaf())
		return dvalue((x-center)/radius);
	if(x<=child[0].right)
		return child[0].du(x);
	if(x>=child[1].left)
		return child[1].du(x);
	return NAN;
}
double L_interval_binary_tree::sigma(double x)
{
	double center=(right+left)/2,radius=(right-left)/2;
	if(isLeaf())
		return svalue((x-center)/radius);
	if(x<=child[0].right)
		return child[0].sigma(x);
	if(x>=child[1].left)
		return child[1].sigma(x);
	return NAN;
}
L_interval_binary_tree *L_interval_binary_tree::previous()
{
	if(parent==NULL)
		return NULL;
	if(pos[pos.size()-1])
		return (this-1)->end();
	return parent->previous();
}
L_interval_binary_tree *L_interval_binary_tree::end()
{
	if(child==NULL)
		return this;
	return child[1].end();
}
bool L_interval_binary_tree::isLeft()
{
	if(parent==NULL)
		return false;
	return !pos[pos.size()-1];
}
bool L_interval_binary_tree::isRight()
{
	if(parent==NULL)
		return false;
	return pos[pos.size()-1];
}
bool L_interval_binary_tree::isRoot()
{
	return parent==NULL;
}
bool L_interval_binary_tree::isLeaf()
{
	return child==NULL;
}
bool L_interval_binary_tree::add_node()
{
	return add_node((left+right)/2);
}
bool L_interval_binary_tree::add_node(double b)
{
	if(child!=NULL)
		return false;
	child=new L_interval_binary_tree[2];
	child[0].left=left;
	child[1].right=right;
	child[0].right=child[1].left=b;
	child[0].parent=child[1].parent=this;
	child[0].pos=child[1].pos=pos;
	child[0].pos.push_back(false);
	child[1].pos.push_back(true);
	return true;
}
bool L_interval_binary_tree::delete_node()
{
	if(child!=NULL)
		return false;
	delete[] child;
	child=NULL;
	return true;
}
void L_interval_binary_tree::cleanData()
{
	data.abdFlag=false;
	if(!isLeaf())
	{
		child[0].cleanData();
		child[1].cleanData();
	}
}
void L_interval_binary_tree::needUpdate()
{
	data.abdFlag=false;
	if(!isRoot())
		parent->needUpdate();
}
void L_interval_binary_tree::evaluateABD(int K,L_function p,L_function q,L_function f,double a,double c,
										 double zl0,double zl1,double zr0,double zr1,double uia,double uik,double uib)
{
	if(data.abdFlag==true)
		return;
	bool gtype=abs(zl0)<abs(zl1)&&abs(zr0)<abs(zr1);
	auto gl=[gtype,a,c,zl0,zl1,zr0,zr1](double x)
	{
		if(gtype)
			return zl1*cosh(x-a)-zl0*sinh(x-a);
		return zl0*(x-a)-zl1;
	};
	auto gr=[gtype,a,c,zl0,zl1,zr0,zr1](double x)
	{
		if(gtype)
			return zr1*cosh(x-c)-zr0*sinh(x-c);
		return zr0*(x-c)-zr1;
	};
	auto dgl=[gtype,a,c,zl0,zl1,zr0,zr1](double x)
	{
		if(gtype)
			return zl1*sinh(x-a)-zl0*cosh(x-a);
		return zl0;
	};
	auto dgr=[gtype,a,c,zl0,zl1,zr0,zr1](double x)
	{
		if(gtype)
			return zr1*sinh(x-c)-zr0*cosh(x-c);
		return zr0;
	};
	auto p_=[gtype,p](double x)
	{
		return p(x);
	};
	auto q_=[gtype,q](double x)
	{
		if(gtype)
			return q(x)+1;
		return q(x);
	};
	auto f_=[f,p,q,uia,uik,uib](double x)
	{
		return f(x)-uia-p(x)*(2*uia*x+uik)-q(x)*(uia*x*x+uik*x+uib);
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
	vec vtau(K),vp(K),vq(K),vf(K),vpsil(K),vpsir(K),vgl(K),vgr(K),vdgl(K),vdgr(K);
	double center=(right+left)/2,radius=(right-left)/2;
	for(int i=0;i<K;i++)
	{
		double x=center+radius*cos(pi*(K-i-0.5)/K);
		vtau(i)=x;
		vp(i)=p_(x);
		vq(i)=q_(x);
		vf(i)=f_(x);
		vpsil(i)=psil(x);
		vpsir(i)=psir(x);
		vgl(i)=gl(x);
		vgr(i)=gr(x);
		vdgl(i)=dgl(x);
		vdgr(i)=dgr(x);
	}
	mat P(K,K,fill::eye);
	for(int i=0;i<K;i++)
	{
		vec eta=P.col(i);
		L_Chebyshev_integral I1(conv_to<vector<double>>::from(vgl%eta));
		L_Chebyshev_integral I2(conv_to<vector<double>>::from(vgr%eta));
		vec vI1(K),vI2(K);
		for(int j=0;j<K;j++)
		{
			vI1(j)=I1.Fl((vtau(j)-center)/radius)*radius;
			vI2(j)=I2.Fr((vtau(j)-center)/radius)*radius;
		}
		P.col(i)=eta+vpsil%vI1+vpsir%vI2;
	}
	vec P_f=solve(P,vf),P_psil=solve(P,vpsil),P_psir=solve(P,vpsir);
	data.al=L_Chebyshev_integral(conv_to<vector<double>>::from(vgl%P_psil)).integral()*radius;
	data.ar=L_Chebyshev_integral(conv_to<vector<double>>::from(vgr%P_psil)).integral()*radius;
	data.bl=L_Chebyshev_integral(conv_to<vector<double>>::from(vgl%P_psir)).integral()*radius;
	data.br=L_Chebyshev_integral(conv_to<vector<double>>::from(vgr%P_psir)).integral()*radius;
	data.dl=L_Chebyshev_integral(conv_to<vector<double>>::from(vgl%P_f)).integral()*radius;
	data.dr=L_Chebyshev_integral(conv_to<vector<double>>::from(vgr%P_f)).integral()*radius;
	data.abdFlag=true;
	data.tau=vtau;
	data.P_f=P_f;
	data.P_psil=P_psil;
	data.P_psir=P_psir;
	data.gl=vgl;
	data.gr=vgr;
	data.dgl=vdgl;
	data.dgr=vdgr;
	data.s=s;
}
bool L_interval_binary_tree::evaluateABD()
{
	if(data.abdFlag)
		return true;
	if(isLeaf())
		return false;
	if(!child[0].evaluateABD())
		return false;
	if(!child[1].evaluateABD())
		return false;
	const L_interval_data &D=child[0].data;
	const L_interval_data &E=child[1].data;
	double Delta=1-E.ar*D.bl;
	data.al=(1-E.al)*(D.al-D.bl*E.ar)/Delta+E.al;
	data.ar=E.ar*(1-D.br)*(1-D.al)/Delta+D.ar;
	data.bl=D.bl*(1-E.br)*(1-E.al)/Delta+E.bl;
	data.br=(1-D.br)*(E.br-D.bl*E.ar)/Delta+D.br;
	data.dl=(1-E.al)/Delta*D.dl+E.dl+(E.al-1)*D.bl/Delta*E.dr;
	data.dr=(1-D.br)/Delta*E.dr+D.dr+(D.br-1)*E.ar/Delta*D.dl;
	data.abdFlag=true;
	return true;
}
bool L_interval_binary_tree::evaluateMu()
{
	if(isRoot())
	{
		data.m=1;
		data.ml=0;
		data.mr=0;
	}
	if(isLeaf())
	{
		data.sigma=data.ml*data.P_psil+data.mr*data.P_psir+data.m*data.P_f;
		int K=data.sigma.n_rows;
		svalue.interpolation(conv_to<vector<double>>::from(data.sigma));
		data.S=abs(svalue.f[K-2])+abs(svalue.f[K-1]-svalue.f[K-3]);
		return true;
	}
	L_interval_data &D=child[0].data;
	L_interval_data &E=child[1].data;
	D.m=data.m;
	E.m=data.m;
	D.ml=data.ml;
	E.mr=data.mr;
	vec::fixed<2> v=solve(mat::fixed<2,2>({{1,E.ar},{D.bl,1}}),
						  vec::fixed<2>{data.mr*(1-E.br)-data.m*E.dr,data.ml*(1-D.al)-data.m*D.dl});
	D.mr=v(0);
	E.ml=v(1);
	if(!child[0].evaluateMu())
		return false;
	if(!child[1].evaluateMu())
		return false;
	return true;
}
void L_interval_binary_tree::evaluateU(int K,double uia,double uik,double uib)
{
	auto ui=[uia,uik,uib](double x)
	{
		return uia*x*x+uik*x+uib;
	};
	double center=(right+left)/2,radius=(right-left)/2;
	L_Chebyshev_integral Il(conv_to<vector<double>>::from(data.gl%data.sigma));
	L_Chebyshev_integral Ir(conv_to<vector<double>>::from(data.gr%data.sigma));
	vec vIl(K),vIr(K),vui(K),vdui(K);
	for(int i=0;i<K;i++)
	{
		vIl(i)=Il.Fl((data.tau(i)-center)/radius)*radius;
		vIr(i)=Ir.Fr((data.tau(i)-center)/radius)*radius;
		vui(i)=ui(data.tau(i));
		vdui(i)=2*uia*data.tau(i)+uik;
	}
	data.u=vui+data.gr%(data.Jl+vIl)/data.s+data.gl%(vIr+data.Jr)/data.s;
	data.du=vdui+data.dgr%(data.Jl+vIl)/data.s+data.dgl%(vIr+data.Jr)/data.s;
	value.interpolation(conv_to<vector<double>>::from(data.u));
	dvalue.interpolation(conv_to<vector<double>>::from(data.du));
}

L_solve_on_mesh::L_solve_on_mesh(L_interval_binary_tree *_mesh,int _K,L_function _p,L_function _q,L_function _f,
								 double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr)
{
	mesh=_mesh;
	K=_K;
	p=_p;
	q=_q;
	f=_f;
	a=mesh->begin()->left;
	c=mesh->end()->right;
	zl0=_zl0;
	zl1=_zl1;
	zr0=_zr0;
	zr1=_zr1;
	Gl=_Gl;
	Gr=_Gr;
	if(zl0!=0||zr0!=0)
	{
		vec::fixed<2> v=solve(mat::fixed<2,2>({{zl0*a+zl1,zl0},{zr0*c+zr1,zr0}}),vec::fixed<2>({_Gl,_Gr}));
		uia=0;
		uik=v(0);
		uib=v(1);
	}
	else
	{
		vec::fixed<2> v=solve(mat::fixed<2,2>({{2*zl1*a,zl1},{2*zr1*c,zr1}}),vec::fixed<2>({_Gl,_Gr}));
		uia=v(0);
		uik=v(1);
		uib=0;
	}
	for(auto i=mesh->begin();i!=NULL;i=i->next())
		i->evaluateABD(K,p,q,f,a,c,zl0,zl1,zr0,zr1,uia,uik,uib);
	mesh->evaluateABD();
	mesh->evaluateMu();
	{
		double Jl=0;
		for(auto i=mesh->begin();i!=NULL;i=i->next())
		{
			i->data.Jl=Jl;
			Jl+=i->data.dl+i->data.ml*i->data.al+i->data.mr*i->data.bl;
		}
	}
	{
		double Jr=0;
		for(auto i=mesh->end();i!=NULL;i=i->previous())
		{
			i->data.Jr=Jr;
			Jr+=i->data.dr+i->data.ml*i->data.ar+i->data.mr*i->data.br;
		}
	}
	for(auto i=mesh->begin();i!=NULL;i=i->next())
	{
		i->evaluateU(K,uia,uik,uib);
	}
}
vector<vector<double>> L_solve_on_mesh::tau()
{
	vector<vector<double>> t;
	for(auto i=mesh->begin();i!=NULL;i=i->next())
		t.push_back(conv_to<vector<double>>::from(i->data.tau));
	return t;
}
vector<vector<double>> L_solve_on_mesh::u()
{
	vector<vector<double>> u;
	for(auto i=mesh->begin();i!=NULL;i=i->next())
		u.push_back(conv_to<vector<double>>::from(i->data.u));
	return u;
}
vector<vector<double>> L_solve_on_mesh::sigma()
{
	vector<vector<double>> s;
	for(auto i=mesh->begin();i!=NULL;i=i->next())
		s.push_back(conv_to<vector<double>>::from(i->data.sigma));
	return s;
}

L_adaptive_solver::L_adaptive_solver(L_interval_binary_tree _mesh,int _K,double _C,double _TOL,
									 L_function _p,L_function _q,L_function _f,
									 double _zl0,double _zl1,double _zr0,double _zr1,double _Gl,double _Gr)
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
	C=_C;
	TOL=_TOL;
	mesh=_mesh;
	bool DblFlag=false;
	double TEST=2*TOL;
	vector<vector<double>> Xs,Ys;
	L_solve_on_mesh solver(&mesh,K,p,q,f,zl0,zl1,zr0,zr1,Gl,Gr);
	int t=0;
	for(;TEST>TOL||!DblFlag;)
	{
		cout<<"\t自适应求解器: "<<t++<<endl;
		Xs=solver.tau();
		Ys=solver.u();
		if(TEST>TOL)
		{
			double Sdiv=0;
			for(auto i=mesh.begin();i!=NULL;i=i->next())
				if(i->data.S>Sdiv)
					Sdiv=i->data.S;
			Sdiv/=pow(2,C);
			for(auto i=mesh.begin();i!=NULL;i=i->next())
			{
				if(i->data.S>=Sdiv)
				{
					i->add_node();
					i->needUpdate();
				}
				else if(i->isLeft()&&(i->data.S+(i+1)->data.S)<Sdiv/pow(2,K))
				{
					i=i->parent;
					i->delete_node();
					i->needUpdate();
				}
			}
			DblFlag=false;
		}
		else
		{
			for(auto i=mesh.begin();i!=NULL;i=i->next())
			{
				i->add_node();
				i->needUpdate();
			}
			DblFlag=true;
		}
		solver=L_solve_on_mesh(&mesh,K,p,q,f,zl0,zl1,zr0,zr1,Gl,Gr);
		double sum=0,difference=0;
		for(int i=0;i<Xs.size();i++)
			for(int j=0;j<Xs[i].size();j++)
			{
				double u=mesh(Xs[i][j]);
				double s=abs(u+Ys[i][j]);
				double d=abs(u-Ys[i][j]);
				if(s>sum)
					sum=s;
				if(d>difference)
					difference=d;
			}
		TEST=difference/sum;
	}
}
vector<vector<double>> L_adaptive_solver::tau()
{
	vector<vector<double>> t;
	for(auto i=mesh.begin();i!=NULL;i=i->next())
		t.push_back(conv_to<vector<double>>::from(i->data.tau));
	return t;
}
vector<vector<double>> L_adaptive_solver::u()
{
	vector<vector<double>> u;
	for(auto i=mesh.begin();i!=NULL;i=i->next())
		u.push_back(conv_to<vector<double>>::from(i->data.u));
	return u;
}
vector<vector<double>> L_adaptive_solver::sigma()
{
	vector<vector<double>> s;
	for(auto i=mesh.begin();i!=NULL;i=i->next())
		s.push_back(conv_to<vector<double>>::from(i->data.sigma));
	return s;
}

