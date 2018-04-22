/*
 * FileReader.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef FILE_READER_H_
#define FILE_READER_H_

#include "../ParticleMS.h"
#include "../../mdutils.h"
#include "../ContainerProperties.h"
#include "../GenerateContainer.h"
#include <list>

class FileReader {

public:
	FileReader();
	virtual ~FileReader();

	/**
	 * Fill a particle container
	 * with particles from a xml-file
	 * @param particles the container to fill
	 * @param filename the path to the xml-file containing information to generate the particles
	 * @param Rcutoff distance to set the truncation distance of the LJ-calculation
	 **/
	void readFile(ContProperties* particlesprop, const char* filename, double Rcutoff);

};

#endif /* FILE_READER_H_ */
