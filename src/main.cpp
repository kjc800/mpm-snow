#pragma once

#include "snowSimulator.h"
#include "grid.h"
#include "particle.h"

#include <iostream>
#include <random>
#include <Eigen/Dense>


using namespace std;
using namespace Eigen;

#define SECONDS 2
#define FPS 24
#define RESOLUTION 10
#define SIDE_LENGTH 1
#define DELTA 0.00001

vector<Particle*> buildSnowball()
{
	vector<Particle*> particles;
	double radius = 0.1;

	for (int i = 0; i < 1000; i++)
	{
		double x;
		double y;
		double z;

		do 
		{
			x = ((rand() % 101) - 50) / 100.0;
			y = ((rand() % 101) - 50) / 100.0;
			z = ((rand() % 101) - 50) / 100.0;

		} while (x * x + y * y + z * z > radius * radius);

		//x = ((rand() % 101) - 50) / 100.0;
		//y = ((rand() % 101) - 50) / 100.0;
		//z = ((rand() % 101) - 50) / 100.0;

		Vector3d position, vel;
		position << x, y, z + 0.5;
		vel << 0.0, 0.0, -1;

		Particle* p = new Particle(position, 0.0017, vel, i);
		particles.push_back(p);
	}

	return particles;
}

vector<Particle*> buildSnowblock()
{
	vector<Particle*> particles;

	for (double x = -0.2; x < 0.2; x += 0.03)
	{
		for (double y = 0.4; y < 0.6; y += 0.03)
		{
			for (double z = -0.2; z < 0.0; z += 0.03)
			{
				Vector3d position, vel;
				position << x, y, z + 0.5;
				vel << 0.0, 0.0, -1;

				Particle* p = new Particle(position, 0.0017, vel, 0);
				particles.push_back(p);
			}
		}
	}

	return particles;
}

vector<Particle*> buildSimpleSnow()
{
	vector<Particle*> particles;

	for (double x = -0.05; x < 0.05; x += 0.03)
	{
		double y = 0.403;
		double z = 0.0; 

		Vector3d position, vel;
		position << x, y, z + 0.5;
		vel << 0.0, 0.0, -1;

		Particle* p = new Particle(position, 0.0017, vel, 0);
		particles.push_back(p);

	}

	return particles;
}

vector<Particle*> buildSnowflake()
{
	vector<Particle*> particles;

	Vector3d position, vel;
	position << 0.503, 0.503, 0.503;
	vel << 0.0, 0.0, -1;

	Particle* p = new Particle(position, 0.0017, vel, 0);
	particles.push_back(p);

	

	return particles;
}


 


int main(int argc, char **argv) {

	SnowSimulator snow_simulator = SnowSimulator(SECONDS, 0.001, 24);
	Grid grid = Grid(RESOLUTION, SIDE_LENGTH);

	vector<Particle*> particles = buildSnowball();


	grid.setParticles(particles);

	snow_simulator.init(grid, particles);
	snow_simulator.simulate();

	return 0;
}
