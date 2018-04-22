/*
 * CSVWriter.cpp
 *
 *  Created on: 15.01.2018
 *      Author: nguyenj
 */

#include "CSVWriter.h"
#include "../utils/Vector.h"
#include <cstdlib>
#include <sstream>

using namespace std;

namespace outputWriter {

CSVWriter::CSVWriter(int bucketNum, float sizeX) {
	this->bucketNum = bucketNum;
	this->sizeX = sizeX;
}

CSVWriter::~CSVWriter() {
	// TODO Auto-generated destructor stub
}

void CSVWriter::plotParticles(ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>* particles, const std::string& filename, int iteration) {
	std::ofstream file;
	std::ofstream locfile;
	stringstream strstr;
	stringstream locstrstr;
	strstr << filename << ".csv";
	locstrstr << filename << "Profile.csv." << iteration;


	//open files
	file.open(strstr.str().c_str(), fstream::out | fstream::app);
	locfile.open(locstrstr.str().c_str());

	float bucketSize = sizeX / bucketNum;
	float buckets[bucketNum];
	for (int i = 0; i < bucketNum; i++) {
		buckets[i] = 0.0f;
	}

	//<ParticleMS, FullParticleCell<ParticleMS>>::iterator
	auto iterator = particles->begin();
	while (iterator.isValid()) {
		Particle& p = *iterator;
		std::array<double, 3> x = p.getR();
		std::array<double, 3> v = p.getV();

		int bucket = (int)(x[0] / bucketSize);
		buckets[bucket] += v[1];

		++iterator;
	}

	file << buckets[0];
	locfile << buckets[0];
	for (int i = 1; i < bucketNum; i++) {
		file << ", " << buckets[i];
		locfile << endl <<buckets[i];
	}
	file << endl;

	file.close();
}

} // end namespace
