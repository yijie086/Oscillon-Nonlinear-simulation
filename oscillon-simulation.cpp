#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>

using namespace std;

#define num_threads 4

//=================To solve this pde equation: \frac{\partial^2\phi}{\partial t^2}-\frac{\partial^2\phi}{\partial r^2}-\frac{D-1}{r}\frac{\partial\phi}{\partial r}+\frac{\partial\phi}{\partial r}=0

//=================Input parameters(global variable)
const double R=128;					//calculate space
const double calR=20;
const int zoombeta=2;					//
const int sz=256*zoombeta;				//number of cells
const double dr=R/sz;				//length of cells
const double cfl=0.01;				//value of dt/dx
const double dt=dr*cfl;				//length of time step
const double tnum=8.0*pow(10,10)*zoombeta;	//number of time step
const double tlength=tnum*dt;		//length of time
const double Dim=2;					//dimention
const double g=0.60;					//phi6 coefficient
const double A=1.0;					//intial Configuration A
const double a=81;					//intial Configuration a
const double sigma=90;				//intial Configuration sigma
const double peroid=10000*zoombeta;		//output peroid
const double pi=3.1415926535898;	//const pi
const int numofrT=16;				//number of peroid beacon
double potenial(double phi)
{
	return(0.5*(phi*phi-phi*phi*phi*phi+g*phi*phi*phi*phi*phi*phi));
}
double pdpotenial(double phi)
{
	return(1.0*phi-2.0*phi*phi*phi+3.0*g*phi*phi*phi*phi*phi);
}
int myabs(int i)
{
	if(i>=0)
	{
		return(i);
	}
	else
	{
		return(-i-1);
	}
}
double pdr(double *phi,int i)					//\frac{\partial}{\partial r}
{
	if((sz-i)>=2)
	{
		return((phi[myabs(i-2)]-8.0*phi[myabs(i-1)]+8.0*phi[myabs(i+1)]-phi[myabs(i+2)])/(12.0*dr));
	}
	else
	{
		return((3.0*phi[i-4]-16.0*phi[i-3]+36.0*phi[i-2]-48.0*phi[i-1]+25.0*phi[i])/(12.0*dr));
	}
}
double pd2r(double *phi,int i)					//\frac{\partial^2}{\partial r^2}
{
	if((sz-i)>=2)
	{
		return((-phi[myabs(i-2)]+16.0*phi[myabs(i-1)]-30.0*phi[myabs(i)]+16.0*phi[myabs(i+1)]-phi[myabs(i+2)])/(12.0*dr*dr));
	}
	else
	{
		return((-10*phi[i-5]+61.0*phi[i-4]-156.0*phi[i-3]+214.0*phi[i-2]-154.0*phi[i-1]+45.0*phi[i])/(12.0*dr*dr));
	}
}
void outputphi(double j,double *phi,double peroid)				//output phi
{
	if(ceil(j/peroid)==(j/peroid))
	{
		ofstream phifile;
		phifile.open("phi.dat",ios::app);
		
		phifile<<j*dt<<"\t";						//output physics time
		for(int i=0;i<=sz;i++)
		{
			phifile<<phi[i]<<"\t";
		}
		phifile<<"\n";
	
		phifile.close();
		
		ofstream out("phiflag.dat");				//flag
		for(int i=0;i<=sz;i++)
		{
			out<<phi[i]<<" ";
		}
		out.close();
	}
}
void outputpdtphi(double j,double *pdtphi,double peroid)		//output pdtphi
{
	if(ceil(j/peroid)==(j/peroid))
	{
		ofstream phifile;
		phifile.open("pdtphi.dat",ios::app);
	
		phifile<<j*dt<<"\t";						//output physics time
		for(int i=0;i<=sz;i++)
		{
			phifile<<pdtphi[i]<<"\t";
		}
		phifile<<"\n";
	
		phifile.close();
		
		ofstream out("pdtphiflag.dat");				//flag
		for(int i=0;i<=sz;i++)
		{
			out<<pdtphi[i]<<" ";
		}
		out.close();
	}
}
void outputpdrphi(double j,double *pdrphi,double peroid)		//output pdrphi
{
	if(ceil(j/peroid)==(j/peroid))
	{
		ofstream phifile;
		phifile.open("pdrphi.dat",ios::app);
		
		phifile<<j*dt<<"\t";						//output physics time
		for(int i=0;i<=sz;i++)
		{
			phifile<<pdrphi[i]<<"\t";
		}
		phifile<<"\n";
			
		phifile.close();
		}
}
void outputr(double *r)								//output r
{
	ofstream phifile;
	phifile.open("r.dat",ios::out);
	
	for(int i=0;i<=sz;i++)
	{
		phifile<<r[i]<<"\n";
	}
	
	phifile.close();	
}
void outputrT(double *rT)
{
	ofstream phifile;
	phifile.open("rT.dat",ios::out);
	
	for(int i=0;i<=numofrT;i++)
	{
		phifile<<rT[i]<<"\n";
	}
	
	phifile.close();
}
double formular(double *r,double *phi,double *pdtphi,double *pdrphi,double *pd2rphi,double *pdrtphi,int i)
{
	if((sz-i)<=2)
	{
		return(-0.5*phi[i]-pdrtphi[i]-((Dim-1.0)/(2.0*r[i]))*pdtphi[i]);		//rewrite the pd2tphi on boundry through boundry condition
	}
	else
	{
		return(pd2rphi[i]+((Dim-1.0)/(r[i]))*pdrphi[i]-pdpotenial(phi[i]));
	}

}
double formular2(double r,double phi,double pdtphi,double pdrphi,double pd2rphi,double pdrtphi,double i)
{
	if((sz-i)<=2)
	{
		return(-0.5*phi-pdrtphi-((Dim-1.0)/(2.0*r))*pdtphi);		//rewrite the pd2tphi on boundry through boundry condition
	}
	else
	{
		return(pd2rphi+((Dim-1.0)/(r))*pdrphi-pdpotenial(phi));
	}
	
}
double energycal(double *r,double dr,double *phi,double *pdtphi,double *pdrphi)
{
	double Energy=0;
	for(int i=0;i<=sz;i++)
	{
		Energy=Energy+(0.5*pdrphi[i]*pdrphi[i]+0.5*pdtphi[i]*pdtphi[i]+0.5*phi[i]*phi[i]-0.5*phi[i]*phi[i]*phi[i]*phi[i]+0.5*g*phi[i]*phi[i]*phi[i]*phi[i]*phi[i]*phi[i])*2*pi*r[i]*dr;
	}
	return(Energy);
}
double energycallim(double *r,double dr,double *phi,double *pdtphi,double *pdrphi)
{
	double Energy=0;
	for(int i=0;i<=ceil(calR/dr);i++)
		{
			Energy=Energy+(0.5*pdrphi[i]*pdrphi[i]+0.5*pdtphi[i]*pdtphi[i]+0.5*phi[i]*phi[i]-0.5*phi[i]*phi[i]*phi[i]*phi[i]+0.5*g*phi[i]*phi[i]*phi[i]*phi[i]*phi[i]*phi[i])*2*pi*r[i]*dr;
		}
	return(Energy);
}
void outputenergy(double j,double energy,double energylim,double peroid)
{
	if(ceil(j/peroid)==(j/peroid))
	{
		ofstream phifile;
		phifile.open("energy.dat",ios::app);
		
		phifile<<j*dt<<"\t";						//output physics time
		phifile<<energy<<"\t";
		phifile<<energylim<<"\t";
		phifile<<"\n";
	
		phifile.close();
	}
}
void initialize(double *phi,double *pdtphi,double *r)
{
	for(int i=0;i<=sz;i++)
	{
		phi[i]=A*exp(-pow((r[i])*(r[i])-a,2.0)/(sigma*sigma));		//initial phi setup
		pdtphi[i]=0;												//initial pdtphi setup
	}
}
void continue_initialize(double *phi,double *pdtphi)
{
	fstream phifile;												//continue phi setup
	phifile.open("phiflag.dat",ios::in);
	for(int i=0;i<=sz;i++)
	{
		phifile>>phi[i];
	}
	phifile.close();
	
	fstream pdtphifile;												//continue pdtphi setup
	pdtphifile.open("pdtphiflag.dat",ios::in);
	for(int i=0;i<=sz;i++)
	{
		pdtphifile>>pdtphi[i];
	}
	pdtphifile.close();	
}
void peroidcal(double *phi,double *tempTphi,int j,double *tempT1,double *tempT2,double *T)
{
	for(int i=0;i<=numofrT;i++)
	{
		if((phi[i*(sz/numofrT)]<=0)&&(tempTphi[i*(sz/numofrT)]>=0))
		{
			tempT2[i]=j*dt;
			T[i]=tempT2[i]-tempT1[i];
			tempT1[i]=tempT2[i];
		}
	}
}
void outputperoid(double j,double *T,double peroid)
{
	if(ceil(j/peroid)==(j/peroid))
	{
		ofstream phifile;
		phifile.open("T.dat",ios::app);

		phifile<<j*dt<<"\t";						//output physics time
		for(int i=0;i<=numofrT;i++)
		{
			phifile<<T[i]<<"\t";
		}
		phifile<<"\n";
			
		phifile.close();
	}
}
int main(int argc, char *argv[])
{
	//---------------------initial setup
	double phi[sz+1]={0};
	double r[sz+1]={0};
	double pdrphi[sz+1]={0};
	double pdtphi[sz+1]={0};
	double pd2tphi[sz+1]={0};
	double pd2rphi[sz+1]={0};
	double pdrtphi[sz+1]={0};
	double temphi1[sz+1]={0};
	double tempdtphi[sz+1]={0};
	double tempdrphi[sz+1]={0};
	double tempd2rphi[sz+1]={0};
	double tempdrtphi[sz+1]={0};
	double temphi2[sz+1]={0};
	double k1[sz+1]={0};
	double k2[sz+1]={0};
	double k3[sz+1]={0};
	double k4[sz+1]={0};
	double k1t[sz+1]={0};
	double k2t[sz+1]={0};
	double k3t[sz+1]={0};
	double k4t[sz+1]={0};
	double y1[sz+1]={0};
	double y2[sz+1]={0};
	double y3[sz+1]={0};
	double y1t[sz+1]={0};
	double y2t[sz+1]={0};
	double y3t[sz+1]={0};
	
	double rT[numofrT+1]={0};
	double tempT1[numofrT+1]={0};
	double tempT2[numofrT+1]={0};
	double T[numofrT+1]={0};
	for(int i=0;i<=sz;i++)						//set up squre
	{
		r[i]=(i+0.5)*dr;
	}
	
	//---------------------peroid beacon
	for(int i=0;i<=numofrT;i++)
	{
		rT[i]=r[i*(sz/numofrT)];
	}
	//---------------------
	
	//---------------------initialize
	
	initialize(phi,pdtphi,r);
	//continue_initialize(phi,pdtphi);
												//two initialize can not be used at the same time!
	//---------------------
	outputr(r);
	outputrT(rT);
	//----------------------calculate
	
	for(double j=0;j<=tnum;j++)				//time step calculating
	{
		
		//cout<<j*dt<<"\t"<<phi[1]<<endl;	//test view
		omp_set_num_threads(num_threads);
		#pragma omp parallel
		{
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					pdrphi[i]=pdr(phi,i);
					pd2rphi[i]=pd2r(phi,i);
					pdrtphi[i]=pdr(pdtphi,i);
					pd2tphi[i]=formular(r,phi,pdtphi,pdrphi,pd2rphi,pdrtphi,i);
					temphi1[i]=phi[i];					//save elements for rk4
					tempdtphi[i]=pdtphi[i];
				}											//calculate the value of each parameters in the equation
		}
		outputphi(j,phi,peroid);										//output phi
		outputpdtphi(j,pdtphi,peroid);									//output pdtphi
		outputpdrphi(j,pdrphi,peroid);									//output pdrphi
		outputenergy(j,energycal(r,dr,phi,pdtphi,pdrphi),energycallim(r,dr,phi,pdtphi,pdrphi),peroid);		//output energy
		outputperoid(j,T,peroid);										//output peroid
//-------------------Euler step
/*		#pragma omp parallel
		{
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					phi[i]=phi[i]+pdtphi[i]*dt+0.5*pd2tphi[i]*dt*dt;
					pdtphi[i]=pdtphi[i]+pd2tphi[i]*dt;
				}								//one time step
		}
*/
//-------------------rk4
		#pragma omp parallel
		{
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					k1[i]=dt*pdtphi[i];
					k1t[i]=dt*formular(r,phi,pdtphi,pdrphi,pd2rphi,pdrtphi,i);
					y1[i]=phi[i]+0.5*k1[i];
					y1t[i]=pdtphi[i]+0.5*k1t[i];
				}
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					tempdrphi[i]=pdr(y1,i);
					tempd2rphi[i]=pd2r(y1,i);
					tempdrtphi[i]=pdr(y1t,i);
					k2[i]=dt*y1t[i];
					k2t[i]=dt*formular(r,y1,y1t,tempdrphi,tempd2rphi,tempdrtphi,i);
					y2[i]=phi[i]+0.5*k2[i];
					y2t[i]=pdtphi[i]+0.5*k2t[i];
				}
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					tempdrphi[i]=pdr(y2,i);
					tempd2rphi[i]=pd2r(y2,i);
					tempdrtphi[i]=pdr(y2t,i);
					k3[i]=dt*y2t[i];
					k3t[i]=dt*formular(r,y2,y2t,tempdrphi,tempd2rphi,tempdrtphi,i);
					y3[i]=phi[i]+k3[i];
					y3t[i]=pdtphi[i]+k3t[i];
				}
			#pragma omp for
				for(int i=0;i<=sz;i++)
				{
					tempdrphi[i]=pdr(y3,i);
					tempd2rphi[i]=pd2r(y3,i);
					tempdrtphi[i]=pdr(y3t,i);
					k4[i]=dt*y3t[i];
					k4t[i]=dt*formular(r,y3,y3t,tempdrphi,tempd2rphi,tempdrtphi,i);
					phi[i]=phi[i]+(k1[i]+2.0*k2[i]+2.0*k3[i]+k4[i])/6.0;
					pdtphi[i]=pdtphi[i]+(k1t[i]+2.0*k2t[i]+2.0*k3t[i]+k4t[i])/6.0;
				}
		}
		peroidcal(phi,temphi1,j,tempT1,tempT2,T);
	}
}
