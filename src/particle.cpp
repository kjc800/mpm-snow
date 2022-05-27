#pragma once

#include "particle.h"

using namespace Eigen;

Particle::Particle()
{
	position = Vector3d::Zero();
	mass = 0.0;
	velocity = Vector3d::Zero();
	elastic_deformation_gradient = Matrix3d::Identity();
	plastic_deformation_gradient = Matrix3d::Identity();
	deformation_gradient = Matrix3d::Identity();
}

Particle::Particle(Vector3d p, double m, Vector3d v, int id, 
	Matrix3d e_dg, Matrix3d p_dg, Matrix3d dg)
{
	this->position = p;
	this->mass = m;
	this->velocity = v;
	this->elastic_deformation_gradient = e_dg;
	this->plastic_deformation_gradient = p_dg;
	this->deformation_gradient = dg;
	this->id = id;
}
