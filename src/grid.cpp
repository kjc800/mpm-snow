#pragma once

#include "grid.h"
#include "particle.h"

using namespace std;

/*######################################################*/
/*                    PUBLIC METHODS                    */
/*######################################################*/

Grid::Grid() 
{
	resolution = 0.0;
	side_length = 0.0;
}

Grid::Grid(double resolution, double physical_length)
{
	this->resolution = resolution;
	this->side_length = physical_length;
	this->grid_spacing = side_length / resolution;

	for (int x = 0; x < resolution; x++)
	{
		vector<vector<GridEntry*>> dim1;
		for (int y = 0; y < resolution; y++)
		{
			vector<GridEntry*> dim2;
			for (int z = 0; z < resolution; z++)
			{
				 GridEntry* grid_entry = new GridEntry();
				 if (z * grid_spacing < 0.15 && z * grid_spacing + grid_spacing > 0.15)
				 {
					 grid_entry->isCollision = true;
				 }
				dim2.push_back(grid_entry);
			}
			dim1.push_back(dim2);
		}
		grid.push_back(dim1);
	}
	
}

void Grid::clear()
{
	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				// Clears old entry while preserving the previous velocities
				grid[x][y][z]->mass = 0.0;
				grid[x][y][z]->velocity = Vector3d::Zero();
				grid[x][y][z]->force = Vector3d::Zero();
			}
		}
	}
}

double Grid::getResolution()
{
	return resolution;
}

double Grid::getLength()
{
	return side_length;
}

void Grid::setParticles(vector<Particle*>& p)
{
	particles = p;
}

GridEntry* Grid::entry(int x, int y, int z)
{
	return grid[x][y][z];
}

void Grid::rasterize()
{
	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				setGridMassAndVelocity(x, y, z);
			}
		}
	}
}

void Grid::setVolumes()
{
	// This function should only be called after rasterizing
	// particles with a call to rasterize.

	for (Particle* p : particles)
	{
		// Deprecated for now, weird things happen when grid res is low
		// use density supplied in paper for now.

		/*double density = 0.0;
		for (int x = 0; x < resolution; x++)
		{
			for (int y = 0; y < resolution; y++)
			{
				for (int z = 0; z < resolution; z++)
				{
					double grid_entry_mass = grid[x][y][z]->mass;
					double grid_weight = computeWeight(p, x, y, z);
					density += grid_entry_mass * grid_weight;
				}
			}
		}

		density = density / (double(grid_spacing * grid_spacing * grid_spacing));*/
		p->volume = p->mass / 400;
	}
}

void Grid::setForces(double mu, double lambda, double hardening, Vector3d& external_acceleration)
{
	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				grid[x][y][z]->force = computeForce(x, y, z, mu, lambda, hardening, external_acceleration);
			}
		}
	}
}

void Grid::setGridVelocitiesEuler(double timestep, double mu)
{
	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				updateGridVelocity(x, y, z, timestep, mu);
			}
		}
	}
}

void Grid::setDeformation(double timestep, double compression, double stretch)
{
	for (Particle* p : particles)
	{
		updateDeformation(p, timestep, compression, stretch);
	}
}

void Grid::setParticleVelocities(double alpha, double mu)
{
	for (Particle* p : particles)
	{
		updateParticleVelocity(p, alpha, mu);
	}
}

void Grid::setParticlePositions(double timestep)
{
	for (Particle* p : particles)
	{
		updateParticlePosition(p, timestep);
	}
}

/*######################################################*/
/*                   PRIVATE METHODS                    */
/*######################################################*/

void Grid::setGridMassAndVelocity(int x, int y, int z)
{
	double entry_mass = 0.0;
	Vector3d entry_velocity = Vector3d::Zero();

	for (Particle* p : particles)
	{
		bool hasImpact = isImpactful(p, x, y, z);

		double weight;
		if (hasImpact)
		{
			weight = computeWeight(p, x, y, z);
			p->weights[getHashedKey(x, y, z)] = weight;
			entry_mass += p->mass * weight;
			entry_velocity += p->velocity * p->mass * weight;
		}
	}

	if (entry_mass == 0.0)
	{
		grid[x][y][z]->mass = 0.0;
		grid[x][y][z]->velocity = Vector3d::Zero();
	}
	else
	{
		grid[x][y][z]->mass = entry_mass;
		grid[x][y][z]->velocity = entry_velocity;
	}
}

/*HELPER FUNCTION FOR COMPUTING GRID BASED FORCES*/
Vector3d Grid::computeForce(int x, int y, int z, double mu_0, double lambda_0, double hardening, Vector3d& external_acceleration)
{
	Vector3d force = Vector3d::Zero();
	for (Particle* p : particles)
	{
		bool hasImpact = isImpactful(p, x, y, z);

		if (!hasImpact)
		{
			continue;
		}
		Vector3d weight_grad = computeWeightGradient(p, x, y, z);

		p->weight_gradients[getHashedKey(x, y, z)] = weight_grad;

		// Compute Jacobi SVD to solve for polar 
		JacobiSVD<Matrix3d> svd(p->elastic_deformation_gradient, ComputeFullU | ComputeFullV);


		// Force calcualtion parameters. Re is derived from polar decomposition
		// and Je is the determinant of the plastic deofrmation gradient of p
		Matrix3d Re = svd.matrixU() * svd.matrixV().transpose();

		double Je = p->elastic_deformation_gradient.determinant();
		double Jp = p->plastic_deformation_gradient.determinant();
		double mu = mu_0 * exp(hardening * (1 - Jp));
		double lambda = lambda_0 * exp(hardening * (1 - Jp));


		Matrix3d sum1 = 2 * mu * (p->elastic_deformation_gradient - Re) * p->elastic_deformation_gradient.transpose();
		Matrix3d sum2 = lambda * (Je - 1) * Je * Matrix3d::Identity();
		force += p->volume * (sum1 + sum2) * weight_grad;
	}

	Vector3d final_force = -1 * force + grid[x][y][z]->mass * external_acceleration;
	return final_force;
}

/*HELPER FUNCTION FOR UPDATING GRID BASED VELOCITY*/
void Grid::updateGridVelocity(int x, int y, int z, double timestep, double friction)
{
	// current velocity and mass
	Vector3d v_n = grid[x][y][z]->velocity;
	Vector3d f = grid[x][y][z]->force;
	double m = grid[x][y][z]->mass;

	Vector3d v_next;

	if (m == 0)
	{
		v_next = Vector3d::Zero();
	}
	else
	{
		v_next = v_n + timestep * (1.0 / m) * f;

		if (grid[x][y][z]->isCollision)
		{
			v_next = collideGridEntry(v_next, friction);
		}
	}

	grid[x][y][z]->previousVelocity = v_n;
	grid[x][y][z]->velocity = v_next;
}

Vector3d Grid::collideGridEntry(Vector3d& velocity, double mu)
{
	Vector3d vel = velocity;
	Vector3d vRel = vel - Vector3d::Zero();

	Vector3d floor_n = Vector3d::Zero();
	floor_n[2] = 1.0;

	double v_n = vRel.dot(floor_n);

	if (v_n >= 0)
	{
		return vel;
	}

	Vector3d v_tan = vRel - v_n * floor_n;

	Vector3d v_rel_prime;
	if (v_tan.norm() <= -mu * v_n)
	{
		v_rel_prime = Vector3d::Zero();
	}
	else
	{
		v_rel_prime = v_tan + mu * v_n * (v_tan / v_tan.norm());
	}

	return v_rel_prime;
}


/*HELPER FUNCTION FOR UPDATING PARTICLE DEFORMATION GRADIENTS*/
void Grid::updateDeformation(Particle* p, double timestep, double critical_compression, double critical_stretch)
{
	Matrix3d velocity_gradient = computeVelocityGradient(p);

	Matrix3d next_deformation_grad = (Matrix3d::Identity() + timestep * velocity_gradient) * 
		p->elastic_deformation_gradient * p->plastic_deformation_gradient;

	Matrix3d temp_elastic = (Matrix3d::Identity() + timestep *
		velocity_gradient) * p->elastic_deformation_gradient;

	// Compute Jacobi SVD to solve for polar decomposition
	JacobiSVD<Matrix3d> svd(temp_elastic, ComputeFullU | ComputeFullV);
	Matrix3d sigma;
	Vector3d sv = svd.singularValues();
	
	// clamped singular values
	double sv0 = min(max(sv[0], 1 - critical_compression), 1 + critical_stretch);
	double sv1 = min(max(sv[1], 1 - critical_compression), 1 + critical_stretch);
	double sv2 = min(max(sv[2], 1 - critical_compression), 1 + critical_stretch);
	
	//double fil = min(max(0.0  , 1 - critical_compression), 1 + critical_stretch);

	sigma << sv0, 0.0, 0.0,
			 0.0, sv1, 0.0,
			 0.0, 0.0, sv2;
	
	Matrix3d next_elastic_deformation = svd.matrixU() * sigma * svd.matrixV().transpose();
	Matrix3d next_plastic_deformation = next_elastic_deformation.inverse() * next_deformation_grad;

	p->deformation_gradient = next_deformation_grad;
	p->elastic_deformation_gradient = next_elastic_deformation;
	p->plastic_deformation_gradient = next_plastic_deformation;
}

/*HELPER FUNCTION FOR UPDATING AND COMPUTING PARTICLE VELOCITIES*/
Vector3d Grid::computePICVelocity(Particle* p)
{
	Vector3d velocity = Vector3d::Zero();

	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				bool hasImpact = isImpactful(p, x, y, z);
				if (hasImpact)
				{
					velocity += grid[x][y][z]->velocity * p->weights[getHashedKey(x, y, z)];
				}
			}
		}
	}

	return velocity;
}

Vector3d Grid::computeFLIPVelocity(Particle* p)
{
	Vector3d velocity = p->velocity;

	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				bool hasImpact = isImpactful(p, x, y, z);
				if (hasImpact)
				{
					velocity += (grid[x][y][z]->velocity - grid[x][y][z]->previousVelocity)
						* p->weights[getHashedKey(x, y, z)];
				}
			}
		}
	}

	return velocity;
}

void Grid::updateParticleVelocity(Particle* p, double alpha, double friction)
{
	Vector3d v_next = (1 - alpha) * computePICVelocity(p) + alpha * computeFLIPVelocity(p);

	if (p->position[2] <= 0.16)
	{
		p->velocity = collideParticle(v_next, 0.5);
		p->position[2] = 0.16;
	}
	else
	{
		p->velocity = v_next;
	}
}

Vector3d Grid::collideParticle(Vector3d& velocity, double mu)
{
	Vector3d vel = velocity;
	Vector3d vRel = vel - Vector3d::Zero();

	Vector3d floor_n = Vector3d::Zero();
	floor_n[2] = 1.0;

	double v_n = vRel.dot(floor_n);

	if (v_n >= 0)
	{
		return vel;
	}

	Vector3d v_tan = vRel - v_n * floor_n;

	Vector3d v_rel_prime;
	if (v_tan.norm() <= -mu * v_n)
	{
		v_rel_prime = Vector3d::Zero();
	}
	else
	{
		v_rel_prime = v_tan + mu * v_n * (v_tan / v_tan.norm());
	}

	return v_rel_prime;
}

void Grid::updateParticlePosition(Particle* p, double timestep)
{
	p->position = p->position + timestep * p->velocity;
}


// Weighting function and weighting gradients
double Grid::N(double x)
{
	if (abs(x) < 1 && abs(x) >= 0)
	{
		return 0.5 * abs(x * x * x) - (x * x) + (2.0 / 3.0);
	}
	else if (abs(x) < 2 && abs(x) >= 1)
	{
		return (-1 * abs(x * x * x) / 6.0) + (x * x) - 2 * abs(x) + (4.0 / 3.0);
	}
	return 0.0;
}

double Grid::dN(double x)
{
	if (x >= -2 && x < -1)
	{
		return 0.5 * (x * x) - 2 * abs(x) + 2;
	}
	else if (x >= -1 && x < 0)
	{
		return -1.5 * (x * x) + 2 * abs(x);
	}
	else if (x >= 0 && x < 1)
	{
		return 1.5 * (x * x) - 2 * abs(x);
	}
	else if (x >= 1 && x < 2)
	{
		return -0.5 * (x * x) + 2 * abs(x) - 2;
	}
	return 0.0;
}

bool Grid::isImpactful(Particle* p, int i, int j, int k)
{
	double h = grid_spacing;

	double xp = p->position[0];
	double yp = p->position[1];
	double zp = p->position[2];

	if (abs(xp - i * h) >= 2 * h)
	{
		return false;
	}
	
	if (abs(yp - j * h) >= 2 * h)
	{
		return false;
	}

	if (abs(zp - k * h) >= 2 * h)
	{
		return false;
	}
	return true;
}

int Grid::getHashedKey(int x, int y, int z)
{
	return x * 31 * 31 + y * 31 + z;
}

double Grid::computeWeight(Particle* p, int i, int j, int k)
{
	double coef = 1.0 / grid_spacing;

	// intermediate products in the interpolation weights
	double prod1 = N(coef * (p->position[0] - i * grid_spacing));
	double prod2 = N(coef * (p->position[1] - j * grid_spacing));
	double prod3 = N(coef * (p->position[2] - k * grid_spacing));

	return prod1 * prod2 * prod3;
}

Vector3d Grid::computeWeightGradient(Particle* p, int i, int j, int k)
{
	Vector3d gradient;

	// Some heler values to reduce the length of equations lmao
	double h = grid_spacing;
	double q = 1.0 / h;

	double xp = p->position[0];
	double yp = p->position[1];
	double zp = p->position[2];

	double row0 = dN(q * (xp - i * h)) * N(q * (yp - j * h)) * N(q * (zp - k * h));
	double row1 = N(q * (xp - i * h)) * dN(q * (yp - j * h)) * N(q * (zp - k * h));
	double row2 = N(q * (xp - i * h)) * N(q * (yp - j * h)) *  dN(q * (zp - k * h));

	gradient[0] = row0;
	gradient[1] = row1;
	gradient[2] = row2;

	return gradient;
}

Matrix3d Grid::computeVelocityGradient(Particle* p)
{
	Matrix3d gradient = Matrix3d::Zero();

	for (int x = 0; x < resolution; x++)
	{
		for (int y = 0; y < resolution; y++)
		{
			for (int z = 0; z < resolution; z++)
			{
				bool hasImpact = isImpactful(p, x, y, z);
				if (hasImpact)
				{
					Vector3d velocity = grid[x][y][z]->velocity;
					RowVector3d weight_grad = p->weight_gradients[getHashedKey(x, y, z)].transpose();
					gradient += velocity * weight_grad;
				}
			}
		}
	}
	return gradient;
}

