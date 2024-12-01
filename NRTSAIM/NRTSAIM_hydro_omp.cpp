/*C++ CODE - MANGEAT MATTHIEU - 2024*/
/*NON-RECIPROCAL TWO SPECIES ACTIVE ISING MODEL - CONSTANT SPECIES POPULATION - HYDRODYNAMIC THEORY (FINITE DIFFERENCE METHOD)*/

//////////////////////
///// LIBRAIRIES /////
//////////////////////

//Public librairies.
#include <cstdlib>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <omp.h>

using namespace std;

//Personal libraries.
#include "lib/special_functions.cpp"

double modulo(const double &x, const int &LX);

//////////////////////////////
///// MESHGRID FUNCTIONS /////
//////////////////////////////

void meshInit(vector<double> &RHOA, vector<double> &RHOB, vector<double> &MAGA, vector<double> &MAGB, const double &rho0, const double &mag0, const int &LX, const int &init)
{
	const double alpha=0.25, kappa=0.8;
	const double Nliq=kappa/alpha, Ngas=(1-kappa)/(1-alpha), rhoA=0.5*rho0*(1+mag0), rhoB=0.5*rho0*(1-mag0);
	const int phiB=1-2*init;
	
	for (int x=0; x<LX; x++)
	{
		if (x<alpha*LX)
		{
			RHOA[x]=Nliq*rhoA;
			RHOB[x]=Ngas*rhoB;
			MAGA[x]=Nliq*rhoA;
			MAGB[x]=0.;
		}
		else if (x-0.5*LX>0 and x-0.5*LX<alpha*LX)
		{
			RHOA[x]=Ngas*rhoA;
			RHOB[x]=Nliq*rhoB;
			MAGA[x]=0.;
			MAGB[x]=phiB*Nliq*rhoB;
		}
		else
		{
			RHOA[x]=Ngas*rhoA;
			RHOB[x]=Ngas*rhoB;
			MAGA[x]=0.;
			MAGB[x]=0.;
		}
	}
}

void finiteDiff(vector<double> &RHOA, vector<double> &RHOB, vector<double> &MAGA, vector<double> &MAGB, const double &D0, const double &v0, const double &gamma0, const double &beta, const double &JAB, const double &JBA, const double &r, const int &LX)
{
	static const double rAB=square(JAB)*r, rBA=square(JBA)*r;
	vector<double> DRHOA(LX,0), DRHOB(LX,0), DMAGA(LX,0), DMAGB(LX,0);
		
	#pragma omp parallel for default(shared)
	for (int x=0; x<LX; x++)
	{
			const double rho=RHOA[x]+RHOB[x];
			const double mAeff=2*beta*(MAGA[x]+JAB*MAGB[x])/rho, mBeff=2*beta*(MAGB[x]+JBA*MAGA[x])/rho;
			const int xm=modulo(x-1,LX), xp=modulo(x+1,LX);
			
			DRHOA[x]=D0*(RHOA[xm]+RHOA[xp]-2*RHOA[x]) - v0*(MAGA[xp]-MAGA[xm]);
			DRHOB[x]=D0*(RHOB[xm]+RHOB[xp]-2*RHOB[x]) - v0*(MAGB[xp]-MAGB[xm]);
			DMAGA[x]=D0*(MAGA[xm]+MAGA[xp]-2*MAGA[x]) - v0*(RHOA[xp]-RHOA[xm]) + 2*gamma0*exp(0.5*(r+rAB)/rho)*((RHOA[x]-r/(2*beta))*sinh(mAeff)-MAGA[x]*cosh(mAeff));
			DMAGB[x]=D0*(MAGB[xm]+MAGB[xp]-2*MAGB[x]) - v0*(RHOB[xp]-RHOB[xm]) + 2*gamma0*exp(0.5*(r+rBA)/rho)*((RHOB[x]-r/(2*beta))*sinh(mBeff)-MAGB[x]*cosh(mBeff));
	}
	
	#pragma omp parallel for default(shared)
	for (int x=0; x<LX; x++)
	{
		RHOA[x]+=DRHOA[x];
		RHOB[x]+=DRHOB[x];
		MAGA[x]+=DMAGA[x];
		MAGB[x]+=DMAGB[x];
	}
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
double modulo(const double &x, const int &LX)
{
	if (x<0)
	{
		return x+LX;
	}
	else if (x>=LX)
	{
		return x-LX;
	}
	else
	{
		return x;
	}
}

//Total average on the all space.
double average(const vector<double> &RHO, const int &LX)
{
	double rhoAv=0;
	for (int x=0; x<LX; x++)
	{
		rhoAv+=RHO[x];
	}
	return double(rhoAv)/LX;
}

//Export density.
void exportProfiles(const vector<double> &RHOA, const vector<double> &RHOB, const vector<double> &MAGA, const vector<double> &MAGB, const double &beta, const double &epsilon, const double &rho0, const double &mag0, const double &JAB, const double &JBA, const int &LX, const int &NX, const int &init, const double &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_NRTSAIM_dynamics1d/");
	stringstream ssDensity;
	ssDensity << "./data_NRTSAIM_dynamics1d/NRTSAIM_profile_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_mag0=" << mag0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_init=" << init << "_t=" << t << ".txt";
	string nameDensity = ssDensity.str();
	
	ofstream fileDensity(nameDensity.c_str(),ios::trunc);
	fileDensity.precision(6);
	
	//Write in the file.
	static double dx=double(LX)/NX;
	for (int x=0; x<NX; x++)
	{
		fileDensity << x*dx << " " << RHOA[x] << " " << RHOB[x] << " " << MAGA[x] << " " << MAGB[x] << endl;
	}
	fileDensity.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &epsilon, double &rho0, double &mag0, double &JAB, double &JBA, int &LX, double &dt, double &tmax, int &init, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta=" ))
		{
			beta=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-epsilon=" ))
		{
			epsilon=atof(argv[i]+9);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-mag0=" ))
		{
			mag0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-JAB=" ))
		{
			JAB=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-JBA=" ))
		{
			JBA=atof(argv[i]+5);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-dt=" ))
		{
			dt=atof(argv[i]+4);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-threads=" ))
		{
			THREAD_NUM=atoi(argv[i]+9);
		}
		else
		{
			cerr << "BAD ARGUMENT : " << argv[i] << endl;
			abort();
		}
	}
}

/////////////////////
///// MAIN CODE /////
/////////////////////

int main(int argc, char *argv[])
{
	//UNIT PARAMETERS.
	const double D=1.,r=1.;
	
	//Physical parameters: beta=inverse temperature, epsilon=self-propulsion parameter, rho0=average density, mag0=average species magnetization (divided by rho0), (JAB,JBA)=non-reciprocal couplings, LX=size of the box.
	double beta=1.25, epsilon=0.9, rho0=2.5, mag0=0, JAB=1.,JBA=-1.;
	int LX=1024;
	
	//Numerical parameters: dt=discrete time interval, tmax=maximal time, dx=discrete space interval, init=initial condition, THREAD_NUM=number of threads used in OpenMP.
	double dt=0.001, tmax=10000., dx=0.2;
	int init=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,epsilon,rho0,mag0,JAB,JBA,LX,dt,tmax,init,THREAD_NUM);
	
	//OpenMP threads creation.
	const int OMP_MAX_THREADS=omp_get_max_threads();
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;
	
	//Global variables.
	const int NX=LX/dx;
	const double D0=D*dt/square(dx), v0=D*epsilon*dt/dx, gamma0=dt;

	//Density and magnetization on each node.
	vector<double> RHOA(NX,0), RHOB(NX,0), MAGA(NX,0), MAGB(NX,0);
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_NRTSAIM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_NRTSAIM_averages/NRTSAIM_averages_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_mag0=" << mag0 << "_JAB=" << JAB << "_JBA=" << JBA << "_LX=" << LX << "_init=" << init << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	const int Nsteps=int(tmax/dt), DeltaT=int(5./dt);
	meshInit(RHOA,RHOB,MAGA,MAGB,rho0,mag0,NX,init);
	
	//Time evolution.
	for(int t=0;t<=Nsteps;t++)
	{
		//Export data.
		if (t%DeltaT==0 or t==Nsteps)
		{	
			//Export 1d profiles.
			exportProfiles(RHOA,RHOB,MAGA,MAGB,beta,epsilon,rho0,mag0,JAB,JBA,LX,NX,init,t*dt);
			
			//Export the averaged densities and magnetizations.
			const double rhoA=average(RHOA,NX), rhoB=average(RHOB,NX), magA=average(MAGA,NX), magB=average(MAGB,NX);
			cout << "time=" << t*dt << " -rhoA=" << rhoA << " -rhoB=" << rhoB << " -magA=" << magA << " -magB=" << magB << running_time.TimeRun(" ") << endl;
			fileAverages <<  t*dt << " " << rhoA << " " << rhoB << " " << magA << " " << magB << endl;
		}
		
		//At each time-step update densities and magnetizations.
		finiteDiff(RHOA,RHOB,MAGA,MAGB,D0,v0,gamma0,beta,JAB,JBA,r,NX);
	}
	return 0;
}
