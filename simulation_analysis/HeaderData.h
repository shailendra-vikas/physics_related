#ifndef _HEADERDATA_H
#define _HEADERDATA_H
#include "Common.h"
#include <limits>
#include <string>

/**
	This file keeps the header information of the the files
*/
class HeaderData{
	public:
		//Constructor
		HeaderData(const std::string& name);

		//Destructor
		~HeaderData();
		
		// Access methods for various attribute
		inline std::string getName() const{
			return _filename;
		}

		inline bool isGasPresent() const{
			return _npart[0];
		}

		inline bool isHaloPresent() const{
			return _npart[1];
		}

		inline bool isDiskPresent() const{
			return _npart[2];
		}

		inline bool isBulgePresent() const{
			return _npart[3];
		}

		inline bool isStarPresent() const{
			return _npart[4];
		}

		inline bool isBndryPresent() const{
			return _npart[5];
		}

		inline const unsigned int getNumberForType(const int &type) const {
			return _npart[type];
		}

		inline const double getMassForType(const int &type) const{
			return _mass[type];
		}

		//If there is any mass section present or all the masses are in header
		bool isMassSectionPresent() const;

		//Setter methods for various attribute (no of particles for types ,mass)
		void setNumberParticles(unsigned int*& npart);

		void setMass(double*& mass);

		// Method for printing 
		friend std::ostream& operator<<(std::ostream& os,const HeaderData &headerData);
		
	private:
		const std::string _filename;
		unsigned int *_npart;
		double *_mass;
};



#endif
