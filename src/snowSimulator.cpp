#pragma once

#include "snowSimulator.h"
#include "particle.h"

#include <iostream>

using namespace std;

/*######################################################*/
/*                    PUBLIC METHODS                    */
/*######################################################*/

SnowSimulator::SnowSimulator(size_t seconds, double delta, size_t frames)
{
	time_steps = seconds * (1.0 / delta);
	fps = frames;
	this->delta = delta;
}

void SnowSimulator::init(Grid& g, vector<Particle*>& p)
{
	grid = g;
	particles = p;
}

void SnowSimulator::simulate()
{
	// Simulates FPS frames per second, until it has simulated
	// TIME_STEPS seconds. Should simulate a total of 
	// FPS * TIME_STEPS frames total.
	size_t frame_count = 0;

	for (int time = 0; time < time_steps; time++)
	{
		if (time % 20 == 0)
		{
			outputXML(frame_count);
		}
		//outputXML(frame_count);
		simulateStep(frame_count);
		grid.clear();
		
		std::cout << "rendering step: " << to_string(time) << endl;
		if (time % 20 == 0)
		{
			frame_count++;
		}

	/*	frame_count++;*/
	}
}

void SnowSimulator::simulateStep(size_t frame_num)
{
	// 1. Rasterize Particles
	grid.rasterize();

	// 2. Compute Initial Volume (will be constant)
	if (frame_num == 0)
	{
		grid.setVolumes();
	}

	// 3. Compute and Set Grid Forces
	Vector3d gravity;
	gravity << 0.0, 0.0, -200;
	grid.setForces(MU, LAMBDA, HARDENING, gravity);

	// 4. Update Grid Velocities
	// 5. Grid Based Collisions
	// 6. Integrate Velocities
	grid.setGridVelocitiesEuler(delta, FRICTION);

	// 7. Update Deformation Gradients
	grid.setDeformation(delta, COMPRESSION, STRETCH);

	// 8. Update Particle Velocities
	// 9. Particle Based Collisions
	grid.setParticleVelocities(ALPHA, FRICTION);

	// 10. Update Particle Positions
	grid.setParticlePositions(delta);
}

