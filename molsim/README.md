# PSE Molekulardynamik Gruppe 13
Bachelor-Praktikum - Scientific Computing (PSE) Molekulardynamik.  
This is the repository of group 13.

## Requirements
- xerces-c
- cppunit
- log4cxx

Using linux is advised!

## Running Arguments
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

### Configuration file
Configuration files are xml file which specify needed setting.
These start with a setting-element which contains follwing property-elements:
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

Example file: [setting.xml](setting.xml)

### Input file
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

Example files: 
[eingabe-cuboid.xml](input/eingabe-cuboid.xml)

## Member
Raffael (@ga78jey)  
Jan (@ga58muk)

