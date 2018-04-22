/*
 * ParticleType.h
 *
 *  Created on: Dec 12, 2017
 *      Author: jan
 */

#ifndef PARTICLETYPE_H_
#define PARTICLETYPE_H_

#include <unordered_map>
#include <vector>

/** @brief Class that defines properties of a particle type
 *
 */
class ParticleType {
private:
	double mass;
	double epsilon;
	double sigma;
	double Hpart;
	double RtruncLJ;
	bool fixed;

	static std::vector<ParticleType> types;

	// map each used id an number from 1 to type_count
	static std::unordered_map<int, int> type_map;
	// list of all added types
	static std::vector<int> type_list;
	static std::vector<float> epsilon_pairs;

public:
	/** Constructor
	 * @param mass mass
	 * @param epsilon epsilon for Lennard-Jones
	 * @param sigma sigma for Lennard-Jones
	 * @param RtruncLJ the distance at which the Lennard_JOnes calculation should be truncated
	 */
	ParticleType(double mass=1, double epsilon=5, double sigma=1, double RtruncLJ =3.0, bool fixed = false);
	virtual ~ParticleType();

	/*
	 * set the properties of a type
	 * properties can only set once
	 * @param type type to define
	 * @param mass mass of type
	 * @param epsilon epsilon of type
	 * @param sigma sigma of type
	 * @param RtruncLJ truncate distance for LJ
	 * @param fixed true if particle cannot be moved
	 */
	static void setType(int type, double mass, double epsilon, double sigma, double RtruncLJ, bool fixed);

	/*
	 * type of given id
	 * @param id id of type
	 * @return ParticleType of given id
	 */
	static ParticleType& getType(int id);

	/*
	 * get a list of all defined types
	 * @return int-vector of all defined ids
	 */
	static const std::vector<int>& getTypes();

	/*
	 * get id of type
	 * @param type type to find id of
	 * @return id of given type
	 */
	static int getID(int type);

	/*
	 * get mass of type
	 * @param id id of type
	 * @return mass of this type
	 */
	static double getM(int id);

	/*
	 * get mass of type
	 * @return mass of this type
	 */
	double getM();

	/*
	 * get epsilon of type
	 * @param id id of type
	 * @return epsilon of this type
	 */
	static double getEpsilon(int id);
	/*
	 * get epsilon of two particles of given type
	 * @param id1 id of first type
	 * @param id2 id of second type
	 * @return epsilon of this type
	 */
	static double getEpsilon(int id1, int id2);
	/*
	 * get epsilon of type
	 * @return epsilon of this type
	 */
	double getEpsilon();

	/*
	 * get sigma of type
	 * @param id id of type
	 * @return sigma of this type
	 */
	static double getSigma(int id);
	/*
	 * get sigma of type
	 * @return sigma of this type
	 */
	double getSigma();

	/*
	 * get the minimum that distance of two particles
	 * of a type so that they repel each other
	 * @param id id of type
	 * @return Hpart of this type
	 */
	static double getHpart(int id);
	/*
	 * get the minimum that distance of two particles
	 * of a type so that they repel each other
	 * @return Hpart of this type
	 */
	double getHpart();
	/*
	 * get RtruncLJ of type
	 * @param id id of type
	 * @return RtruncLJ of this type
	 */
	static double getRtruncLJ(int id);
	/*
	 * get RtruncLJ of type
	 * @return RtruncLJ of this type
	 */
	double getRtruncLJ();

	/*
	 * check if particle is movable
	 * @param id id of type
	 * @return true if particle is unmovable
	 */
	static bool getFixed(int id);
	/*
	 * check if particle is movable
	 * @return true if particle is unmovable
	 */
	bool getFixed();


};

inline ParticleType& ParticleType::getType(int id) {
	return types[id];
}
inline const std::vector<int>& ParticleType::getTypes() {
	return type_list;
}
inline int ParticleType::getID(int type) {
	return type_map[type];
}


inline double ParticleType::getM(int id) {
	return types[id].mass;
}
inline double ParticleType::getM() {
	return mass;
}

inline double ParticleType::getEpsilon(int id) {
	return types[id].epsilon;
}
inline double ParticleType::getEpsilon(int id1, int id2) {
	if (id1 == id2) {
		return types[id1].epsilon;
	}
	if (id1 > id2) {
		// skip first id1*(id1-1)/2 pairs
		// the following id2 pairs contain particles of type id2
		return epsilon_pairs[id1*(id1-1)/2 + id2];
	}

	// skip first id2*(id2-1)/2 pairs
	// the following id2 pairs contain particles of type id2
	return epsilon_pairs[id2*(id2-1)/2 + id1];
}
inline double ParticleType::getEpsilon() {
	return epsilon;
}

inline double ParticleType::getSigma(int id) {
	return types[id].sigma;
}
inline double ParticleType::getSigma() {
	return sigma;
}

inline double ParticleType::getHpart(int id) {
	return types[id].Hpart;
}
inline double ParticleType::getHpart() {
	return Hpart;
}

inline double ParticleType::getRtruncLJ(int id) {
	return types[id].RtruncLJ;
}
inline double ParticleType::getRtruncLJ() {
	return RtruncLJ;
}

inline bool ParticleType::getFixed(int id) {
	return types[id].fixed;
}
inline bool ParticleType::getFixed() {
	return fixed;
}

#endif /* PARTICLETYPE_H_ */
