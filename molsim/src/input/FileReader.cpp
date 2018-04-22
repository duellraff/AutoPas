/*
 * FileReader.cpp
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#include "FileReader.h"
#include "../ParticleType.h"
#include "particle_input.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>

using namespace autopas;
using namespace std;

log::Logger logg(log::Info, "../../molsim");

FileReader::FileReader() {
}

FileReader::~FileReader() {
}

/**
 * Parse a vector_t element into an array
 * @param src vector_t xml-element
 * @param dest output array
 **/
void parseVector(vector_t src, std::array<double,3>& dest) {
	dest[0] = src.x();
	dest[1] = src.y();
	dest[2] = src.z();
}

/**
 * Parse a int vector_t element into an array
 * @param src vector_t xml-element
 * @param dest output array
 **/
void parseIntVector(int_vector_t src, std::array<int,3>& dest) {
	dest[0] = src.x();
	dest[1] = src.y();
	dest[2] = src.z();
}

/**
 * Read types of particles from a xml node
 * @param pinput xml node containting particletype_t elements
 * @param Rcutoff the potential RtruncLJ
 **/
void readTypes(const unique_ptr<input_t>& pinput, double Rcutoff) {
	logg.debug() << "Read spheres"<<endl;

	int id = 0;
	double mass;
	double epsilon;
	double sigma;
	double RtruncLJ;
	bool fixed;

	for (input_t::types_input_const_iterator i = pinput->types_input().begin();
		         i != pinput->types_input().end();
		         ++i)
	{
		particletype_t type = *i;
		id = type.id();
		mass = type.mass();
		epsilon = type.epsilon();
		sigma = type.sigma();
		if (type.RtruncLJ().present()) {
			RtruncLJ = type.RtruncLJ().get();
		}
		else {
			RtruncLJ = Rcutoff;
		}

		if (type.fixed().present()) {
			fixed = type.fixed().get();
		} else {
			fixed = false;
		}

		ParticleType::setType(id, mass, epsilon, sigma, RtruncLJ, fixed);
	}
	logg.info() << "Read " << pinput->types_input().size() << " different types."<<endl;
}

/**
 * Read single particles from a xml node
 * @param particlesprop the container to add the particles
 * @param pinput xml node containing single_t-elements
 **/
long unsigned readSingleParticle(ContProperties* particlesprop, const unique_ptr<input_t>& pinput, long unsigned index) {
	logg.debug() << "Read single particles"<<endl;

	std::array<double,3> x = {0,0,0};
	std::array<double,3> f = {0,0,0};
	std::array<double,3> v = {1,1,1};
	int id = 0;

	for (input_t::single_input_const_iterator i = pinput->single_input().begin();
		         i != pinput->single_input().end();
		         ++i)
	{
		single_t single_particle = *i;
		parseVector(single_particle.coord(), x);
		parseVector(single_particle.velocity(),v);
		id = single_particle.type();

		ParticleMS p(x, v, id, index);

		if (single_particle.force().present()) {
			parseVector(single_particle.force().get(), f);
			p.setF(f);
		}
		particlesprop->getPointerContainer()->addParticle(p);
		++index;
	}
	return index;
}

/**
 * Read cuboids of particles from a xml node
 * @param particles the container to add the particles
 * @param pinput xml node containting cuboid_t-elements
 **/
long unsigned readCuboid(ContProperties* particlesprop, const unique_ptr<input_t>& pinput, long unsigned index) {
	logg.info() << "Read cuboids"<<endl;

	std::array<double,3> x = {0,0,0};
	std::array<double,3> v = {1,1,1};
	std::array<int,3> dim = {0,0,0};
	double h = 1;
	int id = 0;

	for (input_t::cuboid_input_const_iterator i = pinput->cuboid_input().begin();
		         i != pinput->cuboid_input().end();
		         ++i)
	{
		cuboid_t cuboid_particle = *i;
		parseVector(cuboid_particle.coord(), x);
		parseVector(cuboid_particle.velocity(),v);
		parseIntVector(cuboid_particle.dimension(), dim);
		h = cuboid_particle.mesh();
		id = cuboid_particle.type();

		return generateCuboid(particlesprop->getPointerContainer(), *particlesprop, id, x, dim, h, v, index);
	}
}

/**
 * Read membrane of particles from a xml node
 * @param particles the container to add the particles
 * @param pinput xml node containing membrane_t-elements
 **/
long unsigned readMembrane(ContProperties* particlesprop, const unique_ptr<input_t>& pinput, long unsigned index) {
	logg.debug() << "Read membrane"<<endl;

	double k = 1;
	double r0 = 0;
	double t_end = 10000;
	std::array<double,3> f = {0,0,0};
	std::list<std::array<int, 3>> coord_force;
	std::array<double,3> x = {0,0,0};
	std::array<double,3> v = {1,1,1};
	std::array<int,3> dim = {0,0,0};
	double h = 1;
	int id = 0;

	for (input_t::membrane_input_const_iterator i = pinput->membrane_input().begin();
		         i != pinput->membrane_input().end();
		         ++i)
	{
		membrane_t membrane_particle = *i;
		k = membrane_particle.stiffness();
		r0 = membrane_particle.r_zero();
		if (membrane_particle.force().present()) {
				parseVector(membrane_particle.force().get(), f);
				}
		if (membrane_particle.t_end_force().present()) {
				 t_end = membrane_particle.t_end_force().get();
				}


		for (membrane_t::coord_force_const_iterator j = membrane_particle.coord_force().begin();
				         j != membrane_particle.coord_force().end();
				         ++j)
			{
				int_vector_t coord_f = *j;
				std::array<int,3> fcoord = {0,0,0};
				parseIntVector(coord_f, fcoord);
				coord_force.push_back(fcoord);

			}
		parseVector(membrane_particle.coord(), x);
		parseVector(membrane_particle.velocity(),v);
		parseIntVector(membrane_particle.dimension(), dim);
		h = membrane_particle.mesh();
		id = membrane_particle.type();

		return generateMembrane(particlesprop->getPointerContainer(), *particlesprop, id, k, r0, t_end, f, coord_force, x, dim, h, v, index);
	}
}


/**
 * Read sphere of particles from a xml node
 * @param particles the container to add the particles
 * @param pinput xml node containting sphere_t-elements
 **/
long unsigned readSphere(ContProperties* particlesprop, const unique_ptr<input_t>& pinput, long unsigned index){
		logg.debug() << "Read spheres"<<endl;

	std::array<double,3> x = {0,0,0};
	std::array<double,3> v = {1,1,1};
	int r = 1;
	double h = 1;
	int id = 0;

	for (input_t::sphere_input_const_iterator i = pinput->sphere_input().begin();
		         i != pinput->sphere_input().end();
		         ++i)
	{
		sphere_t sphere_particle = *i;
		parseVector(sphere_particle.coord(), x);
		parseVector(sphere_particle.velocity(),v);
		r = sphere_particle.radius();
		h = sphere_particle.mesh();
		id = sphere_particle.type();

		return generateSphere(particlesprop->getPointerContainer(), *particlesprop, id, x, r, h, v, index);
	}
}


void FileReader::readFile(ContProperties* particlesprop, const char* filename, double Rcutoff) {
	unique_ptr<input_t> pinput (particle_input (filename, xml_schema::flags::dont_validate));

	long unsigned index = 0;
	readTypes(pinput, Rcutoff);
	readSingleParticle(particlesprop, pinput, index);
	readCuboid(particlesprop, pinput, index);
	readMembrane(particlesprop, pinput, index);
	readSphere(particlesprop, pinput, index);
}

