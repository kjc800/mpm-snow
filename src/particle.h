#ifndef PARTICLE_H
#define PARTICLE_H

#pragma once

#include "CGL/CGL.h"
#include "CGL/misc.h"

#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class Particle
{
public:

	// Particle constructors. By default the deformation gradient is the 
	// 3x3 identity matrix
	Particle(Vector3d p, double m, Vector3d v, int id,
		Matrix3d e_dg = Matrix3d::Identity(), 
		Matrix3d p_dg = Matrix3d::Identity(),
		Matrix3d dg =   Matrix3d::Identity());
	Particle();

	Vector3d position;
	Vector3d velocity;
	
	Matrix3d elastic_deformation_gradient;
	Matrix3d plastic_deformation_gradient;
	Matrix3d deformation_gradient;

	double mass;
	double volume;
	
	int id;

	unordered_map<int, double> weights;
	unordered_map<int, Vector3d> weight_gradients;
};

#endif