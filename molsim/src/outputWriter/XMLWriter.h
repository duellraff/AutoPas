/*
 * XMLWriter.h
 *
 *  Created on: Dec 14, 2017
 *      Author: jan
 */

#ifndef SRC_OUTPUTWRITER_XMLWRITER_H_
#define SRC_OUTPUTWRITER_XMLWRITER_H_

#include "../ParticleMS.h"
#include "../input/particle_input.h"

namespace outputWriter {

class XMLWriter {
private:
	input_t* inputFile;

public:
	XMLWriter();
	virtual ~XMLWriter();

	/**
	 * set up internal data structures and prepare to plot a particle.
	 */
	void initializeOutput();

	/**
	 * plot type, mass, position, velocity and force of a particle.
	 *
	 * @note: initializeOutput() must have been called before.
	 */
	void plotParticle(ParticleMS& p);

	/**
	 * writes the final output file.
	 *
	 * @param filename the name of the file to be written.
	 */
	void writeFile(const std::string& filename);
};

} /* namespace outputWriter */

#endif /* SRC_OUTPUTWRITER_XMLWRITER_H_ */
