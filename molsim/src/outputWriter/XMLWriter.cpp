/*
 * XMLWriter.cpp
 *
 *  Created on: Dec 14, 2017
 *      Author: jan
 */

#include "XMLWriter.h"
#include "../ParticleType.h"
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

static autopas::log::Logger logg(log::Info, "../../molsim");

namespace outputWriter {

XMLWriter::XMLWriter() {
	inputFile = NULL;
}

XMLWriter::~XMLWriter() {
	// TODO Auto-generated destructor stub
}

void XMLWriter::initializeOutput() {

	inputFile = new input_t();

	// write all used types
	vector<int> types = ParticleType::getTypes();
	for (vector<int>::iterator it = types.begin(); it != types.end(); ++it) {
		int id = *it;
		ParticleType& pt = ParticleType::getType(id);
		particletype_t pt_xml(id, pt.getM(), pt.getSigma(), pt.getEpsilon());
		inputFile->types_input().push_back(pt_xml);
	}

}

void XMLWriter::writeFile(const std::string& filename) {

	logg.debug() <<"Started xml file writer"<<endl;
	std::ofstream file (filename.c_str());
	particle_input(file, *inputFile);
	delete inputFile;
}

void XMLWriter::plotParticle(ParticleMS& p) {

	std::array<double, 3> x = p.getR();
	vector_t pos(x[0], x[1], x[2]);

	std::array<double, 3> f = p.getF();
	vector_t force(f[0], f[1], f[2]);

	std::array<double, 3> v = p.getV();
	vector_t velocity(v[0], v[1], v[2]);

	single_t sinput = single_t(pos, velocity, p.getType());
	sinput.force(force);
	inputFile->single_input().push_back(sinput);
}

} /* namespace outputWriter */
