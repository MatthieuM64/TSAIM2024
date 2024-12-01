/*C++ CODE - MANGEAT MATTHIEU - 2024*/
/*TWO SPECIES ACTIVE ISING MODEL - WITH SPECIES INTERACTION*/

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
#include "lib/random_OMP.cpp"
#include "lib/special_functions.cpp"

double modulo(const double &x, const double &LX);
void modulo(int &x, const int &LX);

//////////////////////////
///// PARTICLE CLASS /////
//////////////////////////

class particle
{
	public:
	
	int x, y; //position.
	int spin; //spin = +/-1.
	int species; //species = 0(A) / 1(B).
	
	particle(const int &LX, const int &LY, const int &species0, const int &init);
	void band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX);
	void move(const double &r, const double &D, const double &v, const double &theta, const int &LX, const int &LY);
	void spinFlip();
	void speciesFlip();
};

//Creation of bands
void particle::band_formation(const double &x0, const double &alpha, const double &kappa, const int &LX)
{
	const double r=ran();
	if (r<kappa)
	{
		double XX=modulo(x0+alpha*r/kappa,1);
		x=int(XX*LX);
	}
	else
	{
		double XX=modulo(x0+(1-alpha)*(r-kappa)/(1-kappa)+alpha,1);
		x=int(XX*LX);
		spin=2*int(2*ran())-1; //Gaseous phase -> change the spin to random.
	}
}

//Creation of the particle.
particle::particle(const int &LX, const int &LY, const int &species0, const int &init)
{
	species=species0;
	static const double alpha=0.25, kappa=0.8;
	if (init==0) //Disordered state.
	{
		spin=2*int(2*ran())-1;
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else if (init==1) //APF band state.
	{
		spin=1-2*species;
		//x=int(LX/2*species+ran()*LX/4);
		band_formation(0.5*species,alpha,kappa,LX);
		y=int(ran()*LY);
	}
	else if (init==2) //PF band state.
	{
		spin=+1;
		//x=int(LX/2*species+ran()*LX/4);
		band_formation(0.5*species,alpha,kappa,LX);
		y=int(ran()*LY);
	}
	else if (init==3) //High-density APF state.
	{
		spin=1-2*species;
		if (species==0)
		{
			x=int(0.5*ran()*LX);
		}
		else
		{
			x=int(0.5*(1+ran())*LX);
		}
		y=int(ran()*LY);
	}
	else if (init==4) //High-density PF state.
	{
		spin=+1;
		if (species==0)
		{
			x=int(0.5*ran()*LX);
		}
		else
		{
			x=int(0.5*(1+ran())*LX);
		}
		y=int(ran()*LY);
	}
	else if (init==5) //APF liquid state.
	{
		spin=1-2*species;
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else if (init==6) //PF liquid state.
	{
		spin=+1;
		x=int(ran()*LX);
		y=int(ran()*LY);
	}
	else
	{
		cerr << "BAD INIT VALUE: " << init << endl;
		abort();
	}
}

//Update the particle position.
void particle::move(const double &r, const double &D, const double &v, const double &theta, const int &LX, const int &LY)
{
	static double Whop=4*D+theta*v;
	if (r<(D+0.5*(theta+1)*v)/Whop) //Move in favoured direction (spin).
	{
		x+=spin;
	}
	else if (r<(2*D+theta*v)/Whop) //Move in unfavoured direction (-spin).
	{
		x-=spin;
	}
	else if (r<(3*D+theta*v)/Whop) //Move up.
	{
		y++;
	}
	else //Move down.
	{
		y--;
	}
	modulo(x,LX);
	modulo(y,LY);
}

//Flip the spin.
void particle::spinFlip()
{
	spin*=-1;
}

//Flip the species.
void particle::speciesFlip()
{
	species=1-species;
}

///////////////////////////
///// BASIC FUNCTIONS /////
///////////////////////////

//Modulo function.
double modulo(const double &x, const double &LX)
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

//Modulo function.
void modulo(int &x, const int &LX)
{
	if (x<0)
	{
		x+=LX;
	}
	else if (x>=LX)
	{
		x-=LX;
	}
}

//Total average on the whole space.
double average(const vector< vector<int> > &RHO, const int &LX, const int &LY)
{
	int rhoAv=0;
	for (int X0=0; X0<LX; X0++)
	{
		for (int Y0=0; Y0<LY; Y0++)
		{
			rhoAv+=RHO[X0][Y0];
		}
	}
	return double(rhoAv)/(LX*LY);
}

//Average along the y-axis.
vector<double> average1d(const vector< vector<int> > &RHO, const int &LX, const int &LY)
{
	vector<double> rhoAv(LX,0.);
	for (int X0=0; X0<LX; X0++)
	{
		for (int Y0=0; Y0<LY; Y0++)
		{
			rhoAv[X0]+=double(RHO[X0][Y0])/LY;
		}
	}
	return rhoAv;
}

//Export averaged densities and magnetizations along y-axis.
void exportProfiles(const vector<double> &RHOA, const vector<double> &RHOB, const vector<double> &MAGA, const vector<double> &MAGB, const double &beta1, const double &beta2, const double &D, const double &v, const double &theta, const double &rho0, const double &mag0, const double &gamma, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_TSAIM_species_dynamics1d/");
	stringstream ss;
	ss << "./data_TSAIM_species_dynamics1d/TSAIM_species_profile_beta1=" << beta1 << "_beta2=" << beta2 << "_D=" << D << "_v=" << v << "_theta=" << theta << "_rho0=" << rho0 << "_mag0=" << mag0 << "_gamma=" << gamma << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameProfile = ss.str();
	
	ofstream fileProfile(nameProfile.c_str(),ios::trunc);
	fileProfile.precision(6);
	
	//Write in the file.
	fileProfile << "#site\tRHOA\tRHOB\tMAGA\tMAGB" << endl;
	for (int X0=0;X0<LX;X0++)
	{
		fileProfile << X0 << "\t" << RHOA[X0] << "\t" << RHOB[X0] << "\t" << MAGA[X0] << "\t" << MAGB[X0] << endl;	
	}
	fileProfile.close();
}

//Export densities.
void exportDensity(const vector< vector<int> > &RHOA, const vector< vector<int> > &RHOB, const double &beta1, const double &beta2, const double &D, const double &v, const double &theta, const double &rho0, const double &mag0, const double &gamma, const int &LX, const int &LY, const int &init, const int &RAN, const int &t)
{
	//Creation of the file.
	int returnSystem=system("mkdir -p data_TSAIM_species_dynamics2d/");
	stringstream ssDensity;
	ssDensity << "./data_TSAIM_species_dynamics2d/TSAIM_species_density_beta1=" << beta1 << "_beta2=" << beta2 << "_D=" << D << "_v=" << v << "_theta=" << theta << "_rho0=" << rho0 << "_mag0=" << mag0 << "_gamma=" << gamma << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << "_t=" << t << ".txt";
	string nameDensity = ssDensity.str();
	
	ofstream fileDensity(nameDensity.c_str(),ios::trunc);
	fileDensity.precision(6);
		
	//Write in the file.
	for (int Y0=0;Y0<LY;Y0++)
	{
		for (int X0=0; X0<LX; X0++)
		{
			const double rhoA=RHOA[X0][Y0], rhoB=RHOB[X0][Y0];
			if (rhoA>rhoB or (rhoA==rhoB and ran()<0.5))
			{
				fileDensity << rhoA << "\t";
			}
			else
			{
				fileDensity << -rhoB << "\t";
			}
		}
		fileDensity << endl;
	}
	fileDensity.close();
}

//Read parameters from command line.
void ReadCommandLine(int argc, char** argv, double &beta1, double &beta2, double &D, double &v, double &theta, double &rho0, double &mag0, double &gamma, int &LX, int &LY, int &init, int &tmax, int &RAN, int &THREAD_NUM)
{
 	for( int i = 1; i<argc; i++ )
	{
		if (strstr(argv[i], "-beta1=" ))
		{
			beta1=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-beta2=" ))
		{
			beta2=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-D=" ))
		{
			D=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-v=" ))
		{
			v=atof(argv[i]+3);
		}
		else if (strstr(argv[i], "-theta=" ))
		{
			theta=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-rho0=" ))
		{
			rho0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-mag0=" ))
		{
			mag0=atof(argv[i]+6);
		}
		else if (strstr(argv[i], "-gamma=" ))
		{
			gamma=atof(argv[i]+7);
		}
		else if (strstr(argv[i], "-LX=" ))
		{
			LX=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-LY=" ))
		{
			LY=atoi(argv[i]+4);
		}
		else if (strstr(argv[i], "-init=" ))
		{
			init=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-tmax=" ))
		{
			tmax=atoi(argv[i]+6);
		}
		else if (strstr(argv[i], "-ran=" ))
		{
			RAN=atoi(argv[i]+5);
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
	//Physical parameters: beta1=inverse temperature of spin-interaction, beta2=inverse temperature of species-interaction, D=diffusion, v=self-propulsion velocity, theta=model parameter, rho0=average density, mag0=average species magnetization (divided by rho0), gamma=gamma2/gamma1, LX*LY=size of the box.
	double beta1=0.75, beta2=1.5, D=1., v=1.8, theta=0., rho0=10., mag0=0., gamma=0.5;
	int LX=1024, LY=128;
	
	//Numerical parameters: init=initial condition, tmax=maximal time, RAN=index of RNG, THREAD_NUM=number of threads used in OpenMP.
	int init=2, tmax=600000, RAN=0, THREAD_NUM=4;

	//Read imported parameters in command line.
	ReadCommandLine(argc,argv,beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,tmax,RAN,THREAD_NUM);
	
	//OpenMP threads creation.
	omp_set_dynamic(0);
	omp_set_num_threads(THREAD_NUM);
	cout << OMP_MAX_THREADS << " maximum threads on this node. " << THREAD_NUM << " threads will be used." << endl;

	//Start the random number generator.
	init_gsl_ran();
	for (int k=0; k<THREAD_NUM; k++)
	{
		gsl_rng_set(GSL_r[k],THREAD_NUM*RAN+k);
	}
	
	//Number of particles.
	int Npart=int(rho0*LX*LY);
	const int NA=int(0.5*Npart*(1+mag0)), NB=Npart-NA;
	
	//Densities and magnetizations on each sites.
	vector< vector<int> > RHOA(LX,vector<int>(LY,0)), RHOB(LX,vector<int>(LY,0)), MAGA(LX,vector<int>(LY,0)), MAGB(LX,vector<int>(LY,0));
	
	//Creation of the file to export global averages.
	const int dossier=system("mkdir -p ./data_TSAIM_species_averages/");
	stringstream ssAverages;
	ssAverages << "./data_TSAIM_species_averages/TSAIM_species_averages_beta1=" << beta1 << "_beta2=" << beta2 << "_D=" << D << "_v=" << v << "_theta=" << theta << "_rho0=" << rho0 << "_mag0=" << mag0 << "_gamma=" << gamma << "_LX=" << LX << "_LY=" << LY << "_init=" << init << "_ran=" << RAN << ".txt";
	string nameAverages = ssAverages.str();
	ofstream fileAverages(nameAverages.c_str(),ios::trunc);
	fileAverages.precision(6);
	
	//Creation of particles.
	vector<particle> ISING;
	for (int i=0;i<NA;i++)
	{
		particle IsingA(LX,LY,0,init);
		RHOA[IsingA.x][IsingA.y]++;
		MAGA[IsingA.x][IsingA.y]+=IsingA.spin;
		ISING.push_back(IsingA);
	}
	
	for (int i=0;i<NB;i++)
	{
		particle IsingB(LX,LY,1,init);
		RHOB[IsingB.x][IsingB.y]++;
		MAGB[IsingB.x][IsingB.y]+=IsingB.spin;
		ISING.push_back(IsingB);
	}
	
	//Time increment.
	const double dt=1./(4*D+theta*v+exp(2*beta1)+gamma*exp(2*beta2));
	
	//Probability to hop.
	const double proba_hop=(4*D+theta*v)*dt;
	
	//Number of particles per threads
	const int Nth=Npart/THREAD_NUM;
	
	//Lock for the update of magnetization
	vector< vector<omp_lock_t> > lock_site(LX,vector<omp_lock_t>(LY));
	for (int x=0; x<LX; x++)
	{
		for (int y=0; y<LY; y++)
		{
			omp_init_lock(&lock_site[x][y]);
		}
	}
	
	//Time evolution (Monte-Carlo steps).
	for(int t=0;t<=tmax;t++)
	{
		//Export data.
		if (t%150==0 or t==tmax)
		{
			//Export 2d snapshots.
			exportDensity(RHOA,RHOB,beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,RAN,t);
			
			//Export 1d profiles, integrated over y.
			const vector<double> rhoA=average1d(RHOA,LX,LY), rhoB=average1d(RHOB,LX,LY), magA=average1d(MAGA,LX,LY), magB=average1d(MAGB,LX,LY);
			exportProfiles(rhoA,rhoB,magA,magB,beta1,beta2,D,v,theta,rho0,mag0,gamma,LX,LY,init,RAN,t);
			
			//Export the averaged densities and magnetizations.
			const double nA=average(RHOA,LX,LY), nB=average(RHOB,LX,LY), mA=average(MAGA,LX,LY), mB=average(MAGB,LX,LY);
			const double rho=nA+nB, mag=(nA-nB)/(nA+nB), vs=(mA+mB)/(nA+nB), va=(mA-mB)/(nA+nB);
			cout << "time=" << t << " -rho=" << rho << " -mag=" << mag << " -vs=" << vs << " -va=" << va << running_time.TimeRun(" ") << endl;
			fileAverages <<  t << " " << rho << " " << mag << " " << vs << " " << va << endl;
		}
		
		//At each time-step move all particles randomly.
		#pragma omp parallel
		{
			const int actual_thread=omp_get_thread_num();
			for (int i=0; i<Nth; i++)
			{	
				//Choose a particle randomly (j) at the site (X0,Y0).
				const int j=Nth*(actual_thread+ran());
				const int X0=ISING[j].x, Y0=ISING[j].y, species=ISING[j].species, spin=ISING[j].spin;
				const int speciesSpin=1-2*species; // +/- value of the species-spin (sA=+1,sB=-1).
				
				//Probabilities to flip the spin and the species.
				omp_set_lock(&lock_site[X0][Y0]);
				const int mag=RHOA[X0][Y0]-RHOB[X0][Y0], va=MAGA[X0][Y0]-MAGB[X0][Y0], rho=RHOA[X0][Y0]+RHOB[X0][Y0];
				omp_unset_lock(&lock_site[X0][Y0]);
				const double proba_spin_flip=exp(-2*beta1*spin*speciesSpin*double(va)/double(rho))*dt;
				const double proba_species_flip=gamma*exp(-2*beta2*speciesSpin*double(mag)/double(rho))*dt;
				const double ptot=proba_hop+proba_spin_flip+proba_species_flip;
				
				//Verify that the probability to wait is positive.
				if (ptot>1)
				{
					cerr << "THE PROBABILITY TO WAIT IS NEGATIVE: " << 1-ptot << " " << RHOA[X0][Y0] << " " << RHOB[X0][Y0] << " " << MAGA[X0][Y0] << " " << MAGB[X0][Y0] << endl;
					abort();
				}
				
				const double random_number=ran();
				//The particle hops: perform the hopping on the particle and modify the densities (const spin).
				if (random_number<proba_hop)
				{
					ISING[j].move(random_number/proba_hop,D,v,theta,LX,LY);
					const int XN=ISING[j].x, YN=ISING[j].y;
					
					//Update densities and magnetizations on old site.
					omp_set_lock(&lock_site[X0][Y0]);
					RHOA[X0][Y0]-=1-species;
					RHOB[X0][Y0]-=species;
					MAGA[X0][Y0]-=spin*(1-species);
					MAGB[X0][Y0]-=spin*species;
					omp_unset_lock(&lock_site[X0][Y0]);
					
					//Update densities and magnetizations on new site.
					omp_set_lock(&lock_site[XN][YN]);
					RHOA[XN][YN]+=1-species;
					RHOB[XN][YN]+=species;
					MAGA[XN][YN]+=spin*(1-species);
					MAGB[XN][YN]+=spin*species;
					omp_unset_lock(&lock_site[XN][YN]);
				}
				//The spin-state flips: perform the flipping on the particle and modify the densities (on-site, ind. of D,v,theta).
				else if (random_number<proba_hop+proba_spin_flip)
				{
					ISING[j].spinFlip();
					//Update the magnetizations (densities are not modified).
					omp_set_lock(&lock_site[X0][Y0]);
					MAGA[X0][Y0]-=2*spin*(1-species);
					MAGB[X0][Y0]-=2*spin*species;
					omp_unset_lock(&lock_site[X0][Y0]);
				}
				//The species-state flips: perform the flipping on the particle and modify the densities (on-site, ind. of D,v,theta).
				else if (random_number<ptot)
				{
					ISING[j].speciesFlip();
					//Update the magnetizations (densities are modified).
					omp_set_lock(&lock_site[X0][Y0]);
					RHOA[X0][Y0]-=speciesSpin;
					RHOB[X0][Y0]+=speciesSpin;
					MAGA[X0][Y0]-=spin*speciesSpin;
					MAGB[X0][Y0]+=spin*speciesSpin;
					omp_unset_lock(&lock_site[X0][Y0]);
				}
				//Else do nothing (proba_wait)
			}		
		}
	}
	return 0;
}
