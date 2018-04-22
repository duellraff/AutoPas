/*
 * CSVWriter.h
 *
 *  Created on: 15.01.2018
 *      Author: nguyenj
 */

#ifndef CSVWRITER_H_
#define CSVWRITER_H_

#include "../ParticleMS.h"
#include "../../mdutils.h"
#include <fstream>
#include <list>

namespace outputWriter {

class CSVWriter {

private:
	int bucketNum;
	float sizeX;

public:
	CSVWriter(int bucketNum, float sizeX);

	virtual ~CSVWriter();

	void plotParticles(ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>>* particles, const std::string& filename, int iteration);

};

}

#endif /* CSVWRITER_H_ */
