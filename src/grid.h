#ifndef GRID_H
#define GRID_H

#pragma once

#include <vector>
#include "CGL/CGL.h"
#include "CGL/misc.h"

#include <unordered_map>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class Particle;

// Internally grid will be represented as a 3D array.
// All grids are cubes. Assume physical dimensions are
// supplied in meters. By convention, the Z axis points upwards.
struct GridEntry 
{
	GridEntry()
	{
		mass = 0.0;
		velocity = Vector3d::Zero();
		previousVelocity = Vector3d::Zero();
		force = Vector3d::Zero();
	}

	~GridEntry() {};

	double mass;
	Vector3d velocity;
	Vector3d previousVelocity;
	Vector3d force;

	bool isCollision = false;
};

class Grid
{
public:
	Grid();
	Grid(double resolution, double physical_length);

	// We can reset the grid for the next time step by simply deleting
	// and replacing all existing grid entries with new blank ones.
	void clear();

	double getResolution();
	double getLength();

	void setParticles(vector<Particle*>& p);
	GridEntry* entry(int x, int y, int z);

	// Core functions called from the base simulator
	void rasterize();
	void setVolumes();
	void setForces(double mu, double lambda, double hardening, Vector3d& external_acceleration);
	void setGridVelocitiesEuler(double timestep, double friction_coef);
	void setDeformation(double timestep, double compression, double stretch);
	void setParticleVelocities(double alpha, double friction_coef);
	void setParticlePositions(double timestep);

private:
	// Functions that are used only to rasterize initial grid cells
	// these are called at the beginning of the grid rasterization process
	void setGridMassAndVelocity(int x, int y, int z);
	double computeMass(int x, int y, int z);
	Vector3d computeVelocity(int x, int y, int z);
	
	// Functions that are used internally to handle all grid parameter updates
	Vector3d computeForce(int x, int y, int z, double mu_0, double lambda_0, double hardening, Vector3d& external_acceleration);
	void updateGridVelocity(int x, int y, int z, double timestep, double friction);
	void updateDeformation(Particle* p, double timestep, double critical_compression, double critical_stretch);

	// Functions that will handle collisions 
	Vector3d collideGridEntry(Vector3d& velocity, double mu);
	Vector3d collideParticle(Vector3d& velocity, double mu);

	// Functions that are used to compute and update PARTICLE (not grid) velocities
	Vector3d computePICVelocity(Particle* p);
	Vector3d computeFLIPVelocity(Particle* p);
	void updateParticleVelocity(Particle* p, double alpha, double friction);

	void updateParticlePosition(Particle* p, double timestep);

	bool isImpactful(Particle* p, int i, int j, int k);
	int getHashedKey(int x, int y, int z);

	// Dyadic B Spline function, used by computeWeight
	double N(double x);
	double dN(double x);

	// Gradient and Weight helper functions
	double computeWeight(Particle* p, int i, int j, int k);
	Vector3d computeWeightGradient(Particle* p, int i, int j, int k);
	Matrix3d computeVelocityGradient(Particle* p);

	// Private Fields that represent the internal grid
	vector<vector<vector<GridEntry*>>> grid;
	vector<Particle*> particles;

	double resolution;
	double side_length;
	double grid_spacing;
};



#endif