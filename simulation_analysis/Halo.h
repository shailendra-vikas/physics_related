#ifndef _HALO_H
#define _HALO_H
#include "Common.h"
#include <set>

typedef std::set<unsigned int> UIntSet;
// Enumeration for the halo properties.
enum haloProperties {
	MASS_TOTAL,
	X_GRP,
	Y_GRP,
	Z_GRP,
	X_DEN,
	Y_DEN,
	Z_DEN,
	MASS_B,
	MASS_GAS,
	MASS_COLD_GAS,
	V_STAR,
	Z_STAR,
	Z_GAS,
	SFR
};

/**
	This class is to contain the Halo related info
*/
class Halo {
	public:
		// Contructor
		Halo();
		// Destructor
		~Halo();

		// Access methods for the Halo Information
		inline const double& getTotalMass() const{
				return _properties[MASS_TOTAL];
		}

		inline const double * getGrpCord() const{
				return &_properties[X_GRP];
		}

		inline const double * getDenCord() const{
				return &_properties[X_DEN];
		}
		
		inline const double& getMassB() const{
				return _properties[MASS_B];
		}

		inline const double& getMassGas() const{
				return _properties[MASS_GAS];
		}

		inline const double& getMassColdGas() const{
				return _properties[MASS_COLD_GAS];
		}

		inline const double& getMassVStar() const{
				return _properties[V_STAR];
		}

		inline const double& getZStar() const{
				return _properties[Z_STAR];
		}

		inline const double& getZGas() const{
				return _properties[Z_GAS];
		}

		inline const double& getSfr() const{
				return _properties[SFR];
		}
		
		inline void insertParticle(const unsigned int& id,const int & type){
			_particles[type].insert(id);
		}

		// Get the virial radius
		const double getVirialRadius() const;
		
		// Print the Halo
		friend std::ostream& operator<<(std::ostream& os,const Halo &halo);
		
		friend class HaloGroup;
	private:
		double _properties[14];
		UIntSet _particles[TYPES_OF_PARTICLE];
};
#endif
