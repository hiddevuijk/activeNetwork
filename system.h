#ifndef GUARD_SYSTEM_H
#define GUARD_SYSTEM_H



#include "ConfigFile.h"

#include <iostream>
#include <vector>
#include <boost/random.hpp>

;

struct System {
public:
	// initialize from ConfigFile object
	System(ConfigFile config);

	// random number generator
	const boost::normal_distribution<double> ndist;
	const boost::uniform_real<double> udist;

	int seed;
	boost::mt19937 rng;		
	boost::variate_generator<boost::mt19937&,
		boost::normal_distribution<double> > rndist;
	boost::variate_generator<boost::mt19937&,
		boost::uniform_real<double> > rudist;

	// fixed system parameters
	unsigned int N;
	double d;
	double k;
	double v0;
	double tau;
	double gamma;
	double m;

	double L;
	double dt;

	double sqrt_2dt;

	// state of the system
	double t;
	std::vector<double> r;
	std::vector<double> v;
	std::vector<double> eta;
	std::vector<double> forces;
	// increment time
	void step();

	void force();
	
};

void System::force()
{
	
	for(unsigned int i=0; i<N; ++i) 
		forces[i] = 0.;

	double f;
	for(unsigned int i=0; i<(N-1); ++i) {
		f = k*(r[i+1] - r[i] - d);
		forces[i] += f;
		forces[i+1] -= f;
	}
	f = k*(r[0] + L - d - r[N-1] );
	forces[0] -= f;
	forces[N-1] += f;
}


void System::step()
{
	force();

	for(unsigned int i=0;i<N;++i) {
		// update eta
		eta[i] += -eta[i]/tau + sqrt_2dt*rndist();
		// update velocities
		v[i] += (-gamma*v[i]*dt + forces[i]*dt + v0*eta[i])/m;
		// update positions
		r[i] += v[i]*dt; 
	}
	t += dt;
}



System::System(ConfigFile config)
:
	ndist(0.,1.),udist(0,1),
	seed(config.read<unsigned int>("seed")),
	rng(seed), rndist(rng,ndist), rudist(rng,udist)
{
	XYZ rr;

	// init parameters
	N = config.read<unsigned int>("N");
	d = config.read<double>("d");
	k = config.read<double>("k");
	v0 = config.read<double>("v0");
	tau = config.read<double>("tau");
	gamma = config.read<double>("gamma");
	m = config.read<double>("m");
	L = N*d;

	dt = config.read<double>("dt");
	sqrt_2dt = std::sqrt(2*dt);

	// init state
	t = 0.0;
	r = std::vector<double>(N);
	v = std::vector<double>(N,0.);
	eta = std::vector<double>(N,0.);
	forces = std::vector<double>(N,0.);

	for(unsigned int i=0; i<N; ++i)
		r[i] = i*d;


}



class Integration {
public:
	Integration( ConfigFile config)
	{
		Nt_init = config.read<unsigned int>("Nt_init");
		Nt = config.read<unsigned int>("Nt");
		sample_freq = config.read<unsigned int>("sample_freq");
		print_freq = config.read<unsigned int>("print_freq");
		t_unit = config.read<unsigned int>("t_unit");
		dt = config.read<double>("dt");
		bs = config.read<double>("bs");
	}

	unsigned int Nt_init;
	unsigned int Nt;
	unsigned int sample_freq;
	unsigned int print_freq;
	unsigned int t_unit;
	double dt;
	double bs;
};


#endif
