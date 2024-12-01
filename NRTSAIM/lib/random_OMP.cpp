/*LIBRARY FOR RANDOM NUMBER GENERATOR*/

#ifndef DEF_RANDOM_MANGEAT_CPP
#define DEF_RANDOM_MANGEAT_CPP

#include <gsl/gsl_randist.h>

const int OMP_MAX_THREADS=omp_get_max_threads();
vector<gsl_rng*> GSL_r(OMP_MAX_THREADS);

void init_gsl_ran()
{
	for (int i=0; i<OMP_MAX_THREADS; i++)
	{
		GSL_r[i]=gsl_rng_alloc(gsl_rng_mt19937);
	}
}

double ran()
{
	double r=gsl_rng_uniform(GSL_r[omp_get_thread_num()]);
	if (r!=0)
	{
		return r;
	}
	else
	{
		return ran();
	}
}

double gaussian()
{
	static bool eval=false;
	static double next;
	if (not eval)
	{
		double phi=2*M_PI*ran();
		double psi=ran();
		double rad=sqrt(-2*log(psi));
		eval=true;
		next=rad*sin(phi);
		return rad*cos(phi);
	}
	else
	{
		eval=false;
		return next;
	}
}

#endif


