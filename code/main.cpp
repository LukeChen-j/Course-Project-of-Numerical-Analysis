#include<cmath>
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<vector>
#include<ctime>
#include"2.h"
#include"1.h"
#include"3.h"

using namespace std;
double derivative(L_function f,double x,double h=1e-5)
{
	auto F=[f,x,h](double d)
	{
		return f(x+d*h);
	};
	return (F(-2)-8*F(-1)+8*F(1)-F(2))/(12*h);

}
void plot(vector<double> X,vector<double> Y,string file,string form)
{
	if(X.size()==0)
	{
		cout<<"size error: X.size()==0\n";
		return;
	}
	if(X.size()!=Y.size())
	{
		cout<<"size error: X.size()!=Y.size()\n";
		return;
	}
	ofstream figure(file,ios::trunc);
	figure<<setprecision(18)<<"x=["<<X[0];
	for(int i=1;i<X.size();i++)
		figure<<','<<X[i];
	figure<<"];\ny=["<<Y[0];
	for(int i=1;i<Y.size();i++)
		figure<<','<<Y[i];
	figure<<"];\n"<<form;
	figure.close();
}
void plot(vector<vector<double>> Xs,vector<vector<double>> Ys,string file,string form)
{
	for(int i=0;i<Xs.size();i++)
		if(Xs[i].size()==0)
		{
			cout<<"size error: Xs["<<i<<"].size()==0\n";
			return;
		}
	for(int i=0;i<Ys.size();i++)
		if(Ys[i].size()==0)
		{
			cout<<"size error: Ys["<<i<<"].size()==0\n";
			return;
		}
	ofstream figure(file,ios::trunc);
	figure<<setprecision(18);
	for(int j=0;j<Xs.size();j++)
	{
		figure<<"x"<<j<<"=["<<Xs[j][0];
		for(int i=1;i<Xs[j].size();i++)
			figure<<','<<Xs[j][i];
		figure<<"];\n";
	}
	for(int j=0;j<Ys.size();j++)
	{
		figure<<"y"<<j<<"=["<<Ys[j][0];
		for(int i=1;i<Ys[j].size();i++)
			figure<<','<<Ys[j][i];
		figure<<"];\n";
	}
	figure<<form;
	figure.close();
}
clock_t TIME;
string timer;
void tic(string s)
{
	timer=s;
	cout<<"开始计时: "<<s<<endl;
	TIME=clock();
}
void toc()
{
	double t=clock();
	cout<<"结束计时: "<<timer<<" 用时 "<<t-TIME<<" ms\n";
}

double p1(double x)
{
	return 1e10*x;
}
double q1(double x)
{
	return -0.5*1e10;
}
double f1(double x)
{
	return 0;
}
double f1(double x,double u,double u_)
{
	return -1000*x*u_;
}
double f2(double x,double u,double u_)
{
	return -10000*x*u_;
}
double f3(double x,double u,double u_)
{
	return x*x*u*u*u-10*u_;
}
double PPP(double x)
{
	return 2e3*x;
}
double QQQ(double x)
{
	return 0;
}
double FFF(double x)
{
	return -2e3*x;
}

int main()
{
	if(1)
	{
		tic("测试 1: Cusp 自适应求解器 ε=1e-10 1e-10精度");
		L_adaptive_solver S1(L_interval_binary_tree(-1,1),16,4,1e-10,p1,q1,f1,1,0,1,0,1,2);
		toc();
		ostringstream text;
		int N1=S1.tau().size();
		cout<<"总区间数(N_1): "<<N1<<endl;
		text<<"plot(x0,y0,'.-'";
		for(int i=1;i<N1;i++)
			text<<",x"<<i<<",y"<<i<<",'.-'";
		text<<");\ngrid on;\n";
		plot(S1.tau(),S1.u(),"g1_1.m",text.str());

		tic("测试 2: Cusp 自适应求解器 ε=1e-10 1e-5精度");
		L_adaptive_solver S2(L_interval_binary_tree(-1,1),16,4,1e-5,p1,q1,f1,1,0,1,0,1,2);
		toc();
		int N2=S2.tau().size();
		cout<<"总区间数(N_2): "<<N2<<endl;

		tic("测试 3: Cusp 直接求解器 ε=1e-10 使用测试 2 中的网格");
		L_direct_solver S3(S2.mesh,16,p1,q1,f1,1,0,1,0,1,2);
		toc();

		vector<double> S4m;
		for(double i=0;i<=2*N2;i++)
			S4m.push_back(-1+i/N2);
		tic("测试 4: Cusp 直接求解器 ε=1e-10 使用 2*N_2 等距节点网格");
		L_direct_solver S4(L_interval_binary_tree(S4m),16,p1,q1,f1,1,0,1,0,1,2);
		toc();

		cout<<"以测试 1 的结果为参考值, 计算测试2、3、4的绝对误差\n";
		double e2=0;
		double e3=0;
		double e4=0;
		vector<vector<double>> S2tau=S2.tau();
		vector<vector<double>> S2u=S2.u();
		for(int i=0;i<S2tau.size();i++)
			for(int j=0;j<S2tau[i].size();j++)
			{
				double e=abs(S2u[i][j]-S1.mesh(S2tau[i][j]));
				if(e>e2)
					e2=e;
			}
		vector<double> S3tau=conv_to<vector<double>>::from(S3.tau);
		vector<double> S3u=conv_to<vector<double>>::from(S3.u);
		for(int i=0;i<S3tau.size();i++)
		{
			double e=abs(S3u[i]-S1.mesh(S3tau[i]));
			if(e>e3)
				e3=e;
		}
		vector<double> S4tau=conv_to<vector<double>>::from(S4.tau);
		vector<double> S4u=conv_to<vector<double>>::from(S4.u);
		for(int i=0;i<S4tau.size();i++)
		{
			double e=abs(S4u[i]-S1.mesh(S4tau[i]));
			if(e>e4)
				e4=e;
		}
		cout<<"e2 = "<<e2<<endl;
		cout<<"e3 = "<<e3<<endl;
		cout<<"e4 = "<<e4<<endl;
	}
	if(1)
	{
		tic("非线性1");
		L_nonlinear_solver lns1=L_nonlinear_solver(L_interval_binary_tree(-1,1),16,4,1e-2,f1,1,0,1,0,{0,1});
//		L_nonlinear_solver lns2=L_nonlinear_solver(L_interval_binary_tree(-1,1),16,4,1e-2,f2,1,0,1,0,{0,1});
		toc();
		vector<double> X,Y;
		for(int i=0;i<1001;i++)
		{
			double x=-1+i*0.002;
			X.push_back(x);
			Y.push_back(lns1.u(x));
		}
		plot(X,Y,"g2_1.m","plot(x,y);\ngrid on;\n");
	}
	if(1)
	{
		tic("非线性3");
		L_nonlinear_solver lns3=L_nonlinear_solver(L_interval_binary_tree(-1,1),16,4,1e-2,f3,0,1,0,1,{0,-0.8,0.5});
		toc();
		vector<double> X,Y;
		for(int i=0;i<1001;i++)
		{
			double x=-1+i*0.002;
			X.push_back(x);
			Y.push_back(lns3.u(x));
		}
		plot(X,Y,"g2_2.m","plot(x,y);\ngrid on;\n");
	}




	return 0;
}
