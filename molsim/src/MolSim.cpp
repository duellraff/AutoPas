
#include "outputWriter/XYZWriter.h"
#include "outputWriter/VTKWriter.h"
#include "outputWriter/XMLWriter.h"
#include "outputWriter/CSVWriter.h"
#include "input/FileReader.h"
#include "GenerateContainer.h"
#include "ContainerProperties.h"
#include "ParticleMS.h"
#include "MolSimFunctors.h"

#include <list>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <typeinfo>
#include <sys/time.h>
#include <cppunit/ui/text/TestRunner.h>
#include "setting/setting.h"

#include "../mdutils.h"
#include "autopasIncludes.h"

#include <omp.h>

using namespace std;
using namespace autopas;

/**** forward declaration ****/
template<typename Container> void plotParticles(int iteration, Container* partcont);
template<typename Container> void plotParticlesXML(Container* partcont, string filename);
template<typename Container> void plotParticlesProfile(Container* partcont, int bucketNum, float sizeX, int iteration);

template<typename Container> void calculateR(Container* cont);
template<typename Container> void calculateV(Container* cont);
template<typename Container> void resetF(Container* cont);

template<class Container> void testEndaddForce(Container* cont, bool membrane_status, double time, double endtime);

double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;

int ContainerAlgorithm = 1;
int ForceComputationMethod = 1;

string outputname = "output/MD_vtk";
string profileFile = "";

static log::Logger logg(log::Info, "../molsim");

#define PROFILE_Y_INTERVAL 10000


/** \mainpage

\section req_sec Requirements

- xerces-c
- cppunit

Using linux is advised! 


\section args_sec Running Arguments

```
export OMP_NUM_THREADS=8
./MolSim filename
```
Run the simulation with configurations defined in 'filename'.
Set the environment variable **OMP_NUM_THREADS** to define the number of threads used.

```
./MolSim -test [TestSuite]
```
Run the testsuites. A specific suite can be passed. Otherwise all suites will be processed.


\subsection cfile_ssec Configuration File

Configuration files are xml file which specify needed setting.
These start with a **setting**-element which contains follwing property-elements:
- **outputname** - Base name of the output files
- **endfile** - Optional: save last particle states into a new xml input-file
- **inputfiles** - File containing information to create particles
- **frequency** - Write frequency of the output files
- **profileFile** - Optional: Write profile of y-velocity into this file
- **profileBucketsX** - Optional: Number of divisions in x-direction to profile y-velocity
- **delta_t** - The time step between iterations
- **t_end** - Time to simulate
- **sigma** - Sigma for Lennard-Jones calculation
- **epsilon** - Epsilon for Lennard-Jones calculation
- **b_factor** - Factor for Maxwell-Boltzmann Distribution
- **g_grav_x** - Acceleration due to gravity in x-direction
- **g_grav_y** - Acceleration due to gravity in y-direction
- **g_grav_z** - Acceleration due to gravity in z-direction
- **domainX** - x-Dimension of domain
- **domainY** - y-Dimension of domain
- **domainZ** - z-Dimension of domain
- **r_cutoff** - Maximum cutoff radius
- **CellType** - Defines the main iteration algorithm: 0 is DirectSum, 1 is Linked Cells (default if this parameter is not specified), 2 is Verlet lists
- **ForceComputationMethod** - Defines the used functor: 0 is no interactions, 1 is Lennard-Jones potential, 2 is Gravity-
- **bc_left** - Boundary condition at x=0: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **bc_lower** - Boundary condition at y=0: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **bc_upper** - Boundary condition at y=domainY: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **bc_right** - Boundary condition at x=domainX: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **bc_front** - Boundary condition at z=0: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **bc_back** - Boundary condition at z=domainZ: 1 is reflecting boundaries, 2 is periodic boundaries, else it is an outflow
- **thermostat** - Optional: Defines a thermostat that keeps the simulation at given temperature. This will override **b_factor**!
  - **initial** - inital temperature
  - **timestep** - number of timesteps after which the thermostat is applied
  - **ignoreY** - Optional: true if the thermostat does not affect the y-component of the velocity
  - **heating** - Optional: defines a temperature change over time
      - **target** - target temperature
      - **temperature_step** - step size in which the temperature should be changed
      - **timestep** - number of timesteps after which the temperature is changed

Example file: setting.xml

\subsection ifile_ssec Input File
Input files are xml files which specify the particles in given simulation.
These start with a particle_input-element which can contain following elements with subelements:
- **types_input** - Defines the properties of a type
  - **id** - The unique identifier of the type
  - **mass** - The mass of a particle of this type
  - **sigma** - The sigma used for the Lennard-Jones formula
  - **epsilon** - The epsilon used for the Lennard-Jones formula
  - **fixed** - Optional: true if particle cannot be moved
  - **RtruncLJ** - Optional: distance at which the Lennard-Jones force calculation should be truncated, if not specified, it will be set to Rcutoff
- **single_input** - Create a single particle
  - **coord** - Position of particle
  - **force** - Optional: Current force enacting on the particle
  - **velocity** - Initial velocity of the partilce
  - **type** - id of the type of the particle
- **cuboid_input** - Create a cuboid of particles
  - **coord** - Position of lower left front corner of the cuboid
  - **dimension** - Number of particles in each dimension
  - **mesh** - Spacing between each particle
  - **velocity** - Initial mean velocity of each particle
  - **type** - id of the type of the particle
- **sphere_input** - Create a sphere of particles
  - **coord** - Center of the sphere
  - **radius** - Number of particles as radius
  - **mesh** - Spacing between each particle
  - **velocity** - Initial mean velocity of each particle
  - **type** - id of the type of the particle
- **membrane_input** - Create a membrane of particles
  - **stiffness** - stiffness k of the membrane
  - **r_zero** - average bond length of a molecule pair
  - **force** - Optional: force which acts on certain molecules
  - **t_end_force** - Optional: time when the optional force stops to act
  - **coord_force** - mesh coordinates of the particles where the force acts
  - **coord** - Position of lower left front corner of the cuboid
  - **dimension** - Number of particles in each dimension
  - **mesh** - Spacing between each particle
  - **velocity** - Initial mean velocity of each particle
  - **type** - id of the type of the particle


Example files: eingabe-cuboid.xml

\section perf_sec Performance
The following graph compares the runtime of Linked cells algorithm to the use of a non-optimized algorithm.  
The tests were runned on 1000, 2000, 4000 and 8000 particles for 5s with a time step of 0.0005s.
 \htmlonly <style>div.image img[src="test.svg"]{width:40cm;}</style> \endhtmlonly 
 \image html performance-test/test.svg "Performance test"
 \image latex performance-test/test.eps "Performance test" width=40cm

\section thread_sec Parallelization
The following graph compares the Molecule-Updates per Second depending
on the number of threads used. The configuration file used was setting-liquid.xml
 \htmlonly <style>div.image img[src="mups.svg"]{width:40cm;}</style> \endhtmlonly 
 \image html performance-test/mups.svg "Performance test"
 \image latex performance-test/mups.eps "Performance test" width=40cm

*/



int main(int argc, char* argsv[]) {


	logg.info() <<"Hello from MolSim for PSE!"<<std::endl;
	if ((argc>1)&&(string(argsv[1]) == "-test")){           // arguments: -test [TestSuite]
		string testsuite = "";
		if (argc>2){
			testsuite = argsv[2];
		}
		logg.info() <<"Start Tests:"<<std::endl;
		CppUnit::TextUi::TestRunner runner;
		//runner.addTest(ParticleContainerTest::suite());
		//runner.addTest(LinkedContainerTest::suite());
		runner.run(testsuite);
		return 0;
	}
	else if (argc != 2) {  // if invalid arguments
		logg.info() <<"Errounous programme call! "<<std::endl;
		logg.info() << "./molsym xml-file"<<std::endl;
		logg.info() << "./molsym -test [TestSuite]"<<std::endl;
		return 1;
	}

	logg.info() <<"Read setting-file..."<<std::endl;
	unique_ptr<setting_t> s (setting (argsv[1], xml_schema::flags::dont_validate));
	delta_t = s->delta_t();
	end_time = s->t_end();

	float b_factor = s->b_factor();

	float domainX = s->domainX();
	float domainY = s->domainY();
	float domainZ = s->domainZ();
	float r_cutoff = s->r_cutoff();

	std::array<int, 6> BoundCond;

	if (s->ContainerType().present()) {
		ContainerAlgorithm = s->ContainerType().get();
	}

	if (s->ForceComputationMethod().present()) {
		ForceComputationMethod = s->ForceComputationMethod().get();
	}

	BoundCond[0] = s->bc_left();
	BoundCond[1] = s->bc_upper();
	BoundCond[2] = s->bc_right();
	BoundCond[3] = s->bc_lower();
	BoundCond[4] = s->bc_front();
	BoundCond[5] = s->bc_back();


	if ((BoundCond[0] == 2) != (BoundCond[2] == 2) ||
			(BoundCond[1] == 2) != (BoundCond[3] == 2) ||
			(BoundCond[4] == 2) != (BoundCond[5] == 2)) {
		logg.info() <<"Error: Opposite sides must also be set to periodic"<<std::endl;
		return 1;
	}


	int frequency = s->frequency();

	outputname = s->outputname();
	string endfile = "";
	if (s->endfile().present()) {
		endfile = s->endfile().get();
	}

	int profileBucketsX = 0;
	if (s->profileFile().present()) {
		profileFile = s->profileFile().get();
		logg.info() <<"Parse y-velocity profile into " << profileFile<<std::endl;

		if (s->profileBucketsX().present()) {
			profileBucketsX = s->profileBucketsX().get();
		}
	}

	float temperature = nan("");
	int T_timestep = 0;
	float TT_target = nan("");
	bool T_ignoreY = false;
	float TT_Tstep = nan("");
	int TT_timestep = 0;
	if (s->thermostat().present()) {
		thermo_t thermo = s->thermostat().get();
		temperature = thermo.initial();
		T_timestep = thermo.timestep();

		b_factor = sqrt(temperature);

		if (thermo.ignoreY().present()) {
			T_ignoreY = thermo.ignoreY().get();
		}

		if (thermo.heating().present()) {
			thermo_target_t thermo_target = thermo.heating().get();
			TT_target = thermo_target.target();
			TT_Tstep = thermo_target.temperature_step();
			TT_timestep = thermo_target.timestep();
		}
	}

	int threads = 1;
	#ifdef _OPENMP
	threads = omp_get_max_threads();
	#endif

	//Functor declarations
	GravityFunctor<ParticleMS, FullParticleCell<ParticleMS>> fgravity(r_cutoff);
	LennardJonesFunctor<ParticleMS, FullParticleCell<ParticleMS>> fLJ(r_cutoff);


	Functor<ParticleMS, FullParticleCell<ParticleMS>> *func;

	if (ForceComputationMethod == 2){
		func = &fgravity;
	}
	else{
		func = &fLJ;
	}
	MembraneFunctor<ParticleMS, FullParticleCell<ParticleMS>> fmembrane();


	std::array<double,3> BoxMin = {0.0, 0.0, 0.0};
	std::array<double,3> BoxMax = {domainX, domainY, domainZ};

	ParticleContainer<ParticleMS, FullParticleCell<ParticleMS>> *particles;
	DirectSum<ParticleMS, FullParticleCell<ParticleMS>> partDS(BoxMin, BoxMax, r_cutoff);
	LinkedCells<ParticleMS, FullParticleCell<ParticleMS>> partLC(BoxMin, BoxMax, r_cutoff);

	//chose the underlying data structure
	if (ContainerAlgorithm == 0){
		particles = &partDS;
	}
	else{
		particles = &partLC;
	}
	//generate the particle container and the corresponding properties
	ContProperties particlesprop = Generate(particles, delta_t, b_factor, &fLJ, r_cutoff, BoxMax, BoundCond, threads);


	logg.info() <<"Read input files:"<<std::endl;
	FileReader fileReader;
	for (setting_t::inputfiles_const_iterator i (s->inputfiles().begin());
	         i != s->inputfiles().end();
	         ++i)
	{
		const char* filename = (*i).c_str();
		logg.info() <<"Read " << filename << "..."<<std::endl;
		fileReader.readFile(&particlesprop,  filename, r_cutoff);
	}

	std::array<double, 3> g;
	g[0] = s->g_grav_x();
	g[1] = s->g_grav_y();
	g[2] = s->g_grav_z();
	particlesprop.setGravity(g);

	logg.info() <<"Total number: " << particlesprop.getNumParticles()<<std::endl;

	double current_time = start_time;

	int iteration = 0;
	if (!isnan(temperature)) {
			setTemperature(particles, temperature, T_ignoreY);
	}
	plotParticles(iteration, particles);

	// calculate the forces once at the beginning
	particles->iteratePairwiseAoS(func);

	//time measurement

	timeval start_timeofday;
	timeval end_timeofday;
	gettimeofday(&start_timeofday, 0);

	//MAIN LOOP
	// for this loop, we assume: current x, current f and current v are known
	while (current_time < end_time) {

		// calculate new x
		calculateR(particles);

		//update container (important for LinkedCells)
		particles->updateContainer();

		// calculate new f
		resetF(particles);
		particles->iteratePairwiseAoS(func);

		// calculate new v
		calculateV(particles);

		// test end of force up (only used in membrane simulation)

		testEndaddForce(particles, particlesprop.getMembraneStatus(), current_time, particlesprop.getTEndForce()) ;

		//adjust velocities according to thermostat
		if (!isnan(TT_target) && temperature != TT_target)
		{
			if (iteration % TT_timestep == 0) {
				if (TT_target - temperature < TT_Tstep) {
					temperature = TT_target;
				} else {
					temperature += TT_Tstep;
				}
			}
		}

		if (!isnan(temperature)) {
			if (iteration % T_timestep == 0) {
				setTemperature(particles, temperature, T_ignoreY);
			}
		}

		if (profileBucketsX > 0 && iteration % PROFILE_Y_INTERVAL == 0) {
			plotParticlesProfile(particles, profileBucketsX, domainX, iteration);
		}

		iteration++;
		logg.debug() <<"Iteration " << iteration << " at time "<<current_time<<" finished."<<std::endl;
		if (iteration % frequency == 0) {
			logg.info() <<"Plot iteration " << iteration<<std::endl;
			plotParticles(iteration, particles);
		}

		current_time += delta_t;
	}

	gettimeofday(&end_timeofday, 0);
	long long int computation_sec = end_timeofday.tv_sec - start_timeofday.tv_sec;
	long long int computation_usec = 1000000*computation_sec + end_timeofday.tv_usec - start_timeofday.tv_usec;

	//get Molecule update rate per second
	double MUPS = 1000000.0*(iteration * particlesprop.getNumParticles())/computation_usec;
	logg.info() <<"Molecule update rate per second: " << MUPS<<std::endl;
	logg.info() <<"Runtime: " << (computation_usec / 1000.0) << "s"<<std::endl;


	if (endfile.compare("")) {
		plotParticlesXML(particles, endfile);
	}

	logg.info() <<"output written. Terminating..."<<std::endl;
	return 0;
}

/**
 * plot the y-profile to a csv file
 */
template<typename Container> void plotParticlesProfile(Container* particles, int bucketNum, float sizeX, int iteration) {
	outputWriter::CSVWriter writer(bucketNum, sizeX);

	writer.plotParticles(particles, profileFile, iteration);
}

/**
 * plot the particles to a vtk-file
 */
template<typename Container> void plotParticles(int iteration, Container* particles) {

	outputWriter::VTKWriter writer;
	long unsigned i = 0;
	for (auto pp = particles->begin(); pp.isValid(); ++pp) {
		++i;
	}
	writer.initializeOutput(i);

	auto iterator = particles->begin();
	while (iterator.isValid()) {
		ParticleMS& p = *iterator;

		writer.plotParticle(p);

		++iterator;
	}

	writer.writeFile(outputname, iteration);
}

/**
 * plot the particles to a xml-file
 */
template<typename Container> void plotParticlesXML(Container* particles, string filename) {

	outputWriter::XMLWriter writer;
	writer.initializeOutput();

	auto iterator = particles->begin();
	while (iterator.isValid()) {
		ParticleMS& p = *iterator;

		writer.plotParticle(p);

		++iterator;
	}

	writer.writeFile(filename);
}




template<typename Container> void calculateR(Container* cont ){

	for (auto p = cont->begin(); p.isValid(); ++p){
			if (p->getFixed()) {
				continue;
			}

			std::array<double, 3> x = p->getR();

			// calculate the position in the next step
			std::array<double, 3> new_x = arrayMath::add(x,
					arrayMath::add(arrayMath::mulScalar(p->getV(), delta_t),
							arrayMath::mulScalar(p->getF(), delta_t * delta_t / (2 * p->getM()))));

			assert(!std::isnan(new_x[0]) && !std::isnan(new_x[1]) && !std::isnan(new_x[2]));
			p->setR(new_x);
	}

}


template<typename Container> void calculateV(Container* cont ){
	// for all particles
	for (auto p = cont->begin(); p.isValid(); ++p){
		if (p->getFixed()) {
			continue;
		}

		 std::array<double, 3> v = p->getV();

		// calculate the velocity in the next step
		v = arrayMath::add(v, arrayMath::mulScalar(arrayMath::add(p->getOldF(), p->getF()), delta_t / (2 * p->getM())));


		assert(!std::isnan(v[0]) && !std::isnan(v[1]) && !std::isnan(v[2]));
		p->setV(v);
	}

}

template<typename Container> void resetF(Container* cont ){
	//for all Particles
	for (auto p = cont->begin(); p.isValid(); ++p){
		//save current force to old force
		p->getOldF() = p->getF();

		//reset forces to the given constant force
		p->setF(p->getConstF());
	}

}


template<class Container> void testEndaddForce(Container* cont, bool membrane_status, double time, double endtime){
	if (membrane_status && (endtime>time)){

		for(auto p = cont->begin(); p.isValid(); ++p){
			//reset constant force
			p->getConstF() = {0.0, 0.0, 0.0};
			p->getConstF() = arrayMath::mulScalar(p->getGrav(), p->getM());
		}
	}
}









