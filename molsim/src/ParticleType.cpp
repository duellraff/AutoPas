/*
 * ParticleType.cpp
 *
 *  Created on: Dec 12, 2017
 *      Author: jan
 */

#include "ParticleType.h"
#include <math.h>
#include <iostream>

std::vector<ParticleType> ParticleType::types;
std::unordered_map<int, int> ParticleType::type_map;
std::vector<int> ParticleType::type_list;
std::vector<float> ParticleType::epsilon_pairs;

ParticleType::ParticleType(double mass, double epsilon, double sigma, double RtruncLJ, bool fixed) {
	this->mass = mass;
	this->epsilon = epsilon;
	this->sigma = sigma;
	this->Hpart = 1.1225 * sigma;
	this->RtruncLJ = RtruncLJ;
	this->fixed = fixed;
}

ParticleType::~ParticleType() {
	// TODO Auto-generated destructor stub
}

void ParticleType::setType(int type, double mass, double epsilon, double sigma, double RtruncLJ, bool fixed) {
	if (type_map.count(type) == 0) {
		// calculate each epsilon-pair with already defined types
		for (int i = 0; i < type_list.size(); ++i) {
			epsilon_pairs.push_back(sqrt(epsilon * getEpsilon(i)));
		}

		type_map[type] = type_list.size();
		types.push_back(ParticleType(mass, epsilon, sigma, RtruncLJ, fixed));
		type_list.push_back(type);
	}
}
