/*C++ CODE - MANGEAT MATTHIEU - 2024*/
/*TWO SPECIES ACTIVE ISING MODEL - CONSTANT SPECIES POPULATION - HYDRODYNAMIC THEORY (FINITE DIFFERENCE METHOD)*/

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

double modulo(const double &x, const int &NX);

//////////////////////////////
///// MESHGRID FUNCTIONS /////
//////////////////////////////

void meshInit(vector<double> &RHO, vector<double> &MAG, vector<double> &VS, vector<double> &VA, const double &rho0, const double &mag0, const int &LX, const int &init)
{
	const double alpha=0.25, kappa=0.8;
	const double Nliq=kappa/alpha, Ngas=(1-kappa)/(1-alpha), rhoA=0.5*rho0*(1+mag0), rhoB=0.5*rho0*(1-mag0);
	const int phiB=1-2*init;

	for (int x=0; x<LX; x++)
	{
		if (x<alpha*LX)
		{
			RHO[x]=Nliq*rhoA+Ngas*rhoB;
			MAG[x]=Nliq*rhoA-Ngas*rhoB;
			VS[x]=Nliq*rhoA;
			VA[x]=Nliq*rhoA;
		}
		else if (x-0.5*LX>0 and x-0.5*LX<alpha*LX)
		{
			RHO[x]=Ngas*rhoA+Nliq*rhoB;
			MAG[x]=Ngas*rhoA-Nliq*rhoB;
			VS[x]=phiB*Nliq*rhoB;
			VA[x]=-phiB*Nliq*rhoB;
		}
		else
		{
			RHO[x]=Ngas*rho0;
			MAG[x]=Ngas*rho0*mag0;
			VS[x]=0.;
			VA[x]=0.;
		}
	}
}

void finiteDiff(vector<double> &RHO, vector<double> &MAG, vector<double> &VS, vector<double> &VA, const double &D0, const double &v0, const double &gamma0, const double &beta, const double &r, const int &LX)
{
	vector<double> DRHO(LX,0), DMAG(LX,0), DVS(LX,0), DVA(LX,0);
	
	#pragma omp parallel for default(shared)
	for (int x=0; x<LX; x++)
	{
			const double va=2*beta*VA[x]/RHO[x];
			const int xm=modulo(x-1,LX), xp=modulo(x+1,LX);
			
			DRHO[x]=D0*(RHO[xm]+RHO[xp]-2*RHO[x]) - v0*(VS[xp]-VS[xm]);
			DMAG[x]=D0*(MAG[xm]+MAG[xp]-2*MAG[x]) - v0*(VA[xp]-VA[xm]);
			DVS[x]=D0*(VS[xm]+VS[xp]-2*VS[x]) - v0*(RHO[xp]-RHO[xm]) + 2*gamma0*exp(0.5*r/RHO[x])*(MAG[x]*sinh(va)-VS[x]*cosh(va));
			DVA[x]=D0*(VA[xm]+VA[xp]-2*VA[x]) - v0*(MAG[xp]-MAG[xm]) + 2*gamma0*exp(0.5*r/RHO[x])*((RHO[x]-r/(2*beta))*sinh(va)-VA[x]*cosh(va));
	}
	
	#pragma omp parallel for default(shared)
	for (int x=0; x<LX; x++)
	{
		RHO[x]+=DRHO[x];
		MAG[x]+=DMAG[x];
		VS[x]+=DVS[x];
		VA[x]+=DVA[x];
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
void exportProfiles(const vector<double> &RHO, const vector<double> &MAG, const vector<double> &VS, const vector<double> &VA, const double &beta, const double &epsilon, const double &rho0, const double &mag0, const int &LX, const int &NX, const int &init, const double &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_TSAIM_dynamics1d/");
	stringstream ssDensity;
	ssDensity << "./data_TSAIM_dynamics1d/TSAIM_profile_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_mag0=" << mag0 << "_LX=" << LX << "_init=" << init << "_t=" << t << ".txt";
	string nameDensity = ssDensity.str();
	
	ofstream fileDensity(nameDensity.c_str(),ios::trunc);
	fileDensity.precision(6);
		
	//Write in the file.
	static double dx=double(LX)/NX;
	for (int x=0; x<NX; x++)
	{
		fileDensity << x*dx << " " << 0.5*(RHO[x]+MAG[x]) << " " << 0.5*(RHO[x]-MAG[x]) << " " << 0.5*(VS[x]+VA[x]) << " " << 0.5*(VS[x]-VA[x]) << endl;
	}
	fileDensity.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta, double &epsilon, double &rho0, double &mag0, int &LX, double &dx, double &dt, double &tmax, int &init, int &THREAD_NUM)
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
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-dx=" ))
		{
			dx=atof(argv[i]+4);
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
	const double D=1., r=1.;
	
	//Physical parameters: beta=inverse temperature, epsilon=self-propulsion parameter, rho0=average density, mag0=average species magnetization (divided by rho0), LX=size of the box.
	double beta=1., epsilon=0.9, rho0=1.01, mag0=0;
	int LX=1024;
	
	//Numerical parameters: dt=discrete time interval, tmax=maximal time, dx=discrete space interval, init=initial condition, THREAD_NUM=number of threads used in OpenMP.
	double dt=0.001, tmax=10000., dx=0.2;
	int init=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta,epsilon,rho0,mag0,LX,dx,dt,tmax,init,THREAD_NUM);
	
	//OpenMP threads creation.
	const int OMP_MAX_THREADS=omp_get_max_threads();
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;
	
	//Global variables.
	const int NX=LX/dx;
	const double D0=D*dt/square(dx), v0=D*epsilon*dt/dx, gamma0=dt;

	//Density and magnetization on each sites.
	vector<double> RHO(NX,0), MAG(NX,0), VS(NX,0), VA(NX,0);
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_TSAIM_averages/");
	stringstream ssAverages;
	ssAverages << "./data_TSAIM_averages/TSAIM_averages_beta=" << beta << "_epsilon=" << epsilon << "_rho0=" << rho0 << "_mag0=" << mag0 << "_LX=" << LX << "_init=" << init << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	const int Nsteps=int(tmax/dt), DeltaT=int(5./dt);
	meshInit(RHO,MAG,VS,VA,rho0,mag0,NX,init);
	
	//Time evolution.
	for(int t=0;t<=Nsteps;t++)
	{
		//Export data.
		if (t%DeltaT==0 or t==Nsteps)
		{
			//Export 1d profiles.
			exportProfiles(RHO,MAG,VS,VA,beta,epsilon,rho0,mag0,LX,NX,init,t*dt);
			
			//Export the averaged densities and magnetizations.
			const double n0=average(RHO,NX), m0=average(MAG,NX), vs=average(VS,NX), va=average(VA,NX);
			cout << "time=" << t*dt << " -rho=" << n0 << " -mag=" << m0 << " -vs=" << vs << " -va=" << va << running_time.TimeRun(" ") << endl;
			fileAverages <<  t*dt << " " << n0 << " " << m0 << " " << vs << " " << va << endl;
			
		}
		
		//At each time-step update densities and magnetizations.
		finiteDiff(RHO,MAG,VS,VA,D0,v0,gamma0,beta,r,NX);
	}
	return 0;
}
