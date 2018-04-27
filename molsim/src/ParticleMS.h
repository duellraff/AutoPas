/*
 * ParticleMS.h
 *
 *  Created on: 23.02.2010
 *      Author: eckhardw
 */

#ifndef PARTICLEMS_H_
#define PARTICLEMS_H_

#include <list>
#include <array>
#include "utils/Vector.h"
#include "ParticleType.h"
#include "../mdutils.h"

using namespace autopas;

class ParticleMS : public Particle {

private:

	/** position r */
	std::array<double, 3> _r;

	/** velocity v */
	std::array<double, 3> _v;

	/** force f */
	std::array<double, 3> _f;

	/** molecule id for SoA	 */
	long unsigned _id;

	/** the force which was effective on this particle */
	std::array<double, 3> old_f;

	/** a constant force always effective on this particle */
	std::array<double, 3> const_f;

	/** the gravity acting on this particle */
	std::array<double, 3> ggrav;

	/** list of pointers to all neighbouring particles */
	std::array<ParticleMS*,4> direct_neighbours;
	std::array<ParticleMS*,4> diagonal_neighbours;

	/** type of the particle. Use it for whatever you want (e.g. to seperate
	 * molecules belonging to different bodies, matters, and so on)
	 */
	int type;
	int type_id;




public:
	ParticleMS(int type = 0, long unsigned ParticleID = 0);

	ParticleMS(const ParticleMS& other);

	ParticleMS(
			// for visualization, we need always 3 coordinates
			// -> in case of 2d, we use only the first and the second
			std::array<double, 3> x_arg,
	        std::array<double, 3> v_arg,
	        int type = 0, long unsigned ParticleID = 0);

	virtual ~ParticleMS();

	std::array<double, 3>& getR();

	std::array<double, 3>& getV();

	std::array<double, 3>& getF();

	std::array<double, 3> setR(std::array<double, 3>& r);

	std::array<double, 3> setV(std::array<double, 3>& v);

	std::array<double, 3> setF(std::array<double, 3>& f);

	std::array<double, 3>& getOldF();

	std::array<double, 3>& getConstF();

	std::array<double, 3>& getGrav();

	void addF(const std::array<double, 3> &f);

	void subF(const std::array<double, 3> &f);

	long unsigned getID();

	void actualiseNeighbours(ParticleMS* newptr);

	double getM();
	double getHpart();
	double getSigma();
	double getRtruncLJ();
	bool getFixed();

	int getType();
	int& getTypeID();

	std::array<ParticleMS*,4>& getDirectNeighbours();

	std::array<ParticleMS*,4>& getDiagonalNeighbours();

	bool operator==(ParticleMS& other);

	std::string toString();

    enum AttributeNames : int { id, posX, posY, posZ, velX, velY, velZ, forceX, forceY, forceZ,
    	oldforceX, oldforceY, oldforceZ, constforceX, constforceY, constforceZ, ggravX, ggravY, ggravZ,
		directNeigh1, directNeigh2, directNeigh3, directNeigh4,
		diagonalNeigh1, diagonalNeigh2, diagonalNeigh3, diagonalNeigh4, typeID};

};

std::ostream& operator<<(std::ostream& stream, ParticleMS& p);

inline std::array<double, 3>& ParticleMS::getR() {
	return _r;
}
inline std::array<double, 3>& ParticleMS::getV() {
	return _v;
}
inline std::array<double, 3>& ParticleMS::getF() {
	return _f;
}
inline std::array<double, 3> ParticleMS::setR(std::array<double, 3>& r) {
	 _r = r;
}
inline std::array<double, 3> ParticleMS::setV(std::array<double, 3>& v) {
	_v = v;
}
inline std::array<double, 3> ParticleMS::setF(std::array<double, 3>& f) {
	_f= f;
}

inline void ParticleMS::addF(const std::array<double, 3> &f) { _f = arrayMath::add(_f, f); }

inline void ParticleMS::subF(const std::array<double, 3> &f) { _f = arrayMath::sub(_f, f); }


inline std::array<double, 3>& ParticleMS::getOldF() {
	return old_f;
}

inline std::array<double, 3>& ParticleMS::getConstF() {
	return const_f;
}

inline std::array<double, 3>& ParticleMS::getGrav(){
	return ggrav;
}

inline double ParticleMS::getM() {
	return ParticleType::getM(type_id);
}
inline double ParticleMS::getHpart() {
	return ParticleType::getHpart(type_id);
}
inline double ParticleMS::getSigma() {
	return ParticleType::getSigma(type_id);
}
inline double ParticleMS::getRtruncLJ() {
	return ParticleType::getRtruncLJ(type_id);
}
inline bool ParticleMS::getFixed() {
	return ParticleType::getFixed(type_id);
}
inline long unsigned ParticleMS::getID() {
	return _id;
}

inline int ParticleMS::getType() {
	return type;
}

inline int& ParticleMS::getTypeID() {
	return type_id;
}

inline std::array<ParticleMS*,4>& ParticleMS::getDirectNeighbours() {
	return direct_neighbours;
}

inline std::array<ParticleMS*,4>& ParticleMS::getDiagonalNeighbours(){
	return diagonal_neighbours;
}


#endif /* ParticleMS_H_ */
