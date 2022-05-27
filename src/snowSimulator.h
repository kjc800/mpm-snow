#ifndef MPM_SNOW_SIMULATOR_H
#define MPM_SNOW_SIMULATOR_H

#pragma once

#include "grid.h"
#include "particle.h"

#include <vector>
#include <Eigen/Dense>

using namespace std;

class SnowSimulator
{
public:
	SnowSimulator(size_t seconds, double delta, size_t frames = 25);
	void init(Grid& g, vector<Particle*>& p);

	void simulate();

private:

	void simulateStep(size_t frame_num);

	// Outputs a Mitsuba 2 formatted XML file
	void outputXML(size_t frame_count);
	void drawSnowParticles(ofstream& frame);


	Grid grid;
	vector<Particle*> particles;

	size_t time_steps;
	size_t fps;
	double delta;
	
	// Parameters of the snow simulation that can be tweaked as desired.
	const double COMPRESSION = 0.025;
	const double STRETCH = 0.0075;
	const double ALPHA = 0.95;
	const double YOUNGS_MODULUS = 140000.0;
	const double HARDENING = 10;
	const double POISSON_RATIO = 0.2;
	const double FRICTION = 0.5;

	const double MU = YOUNGS_MODULUS / (2 * (1 + POISSON_RATIO));
	const double LAMBDA = (YOUNGS_MODULUS * POISSON_RATIO) / ((1 + POISSON_RATIO) * (1 - 2 * POISSON_RATIO));
};


#endif // !MPM_SNOW_SIMULATOR_H