#ifndef _COMMONDETAILS_H
#define _COMMONDETAILS_H
#include "Common.h"

/**
	This class store the detail of the particles ( mass, position, velocity)  and provides methods to access the details
*/

class CommonDetail{
	public:
		// Contructor for the most basic details
		CommonDetail(const double  &mass,const double  &xpos,const double &ypos,const double &zpos);
		//Contructor for the complete details
		CommonDetail(const double  &mass,const double  &xpos,const double &ypos,const double &zpos ,const double  &xvel,const double &yvel,const double &zvel);
		// Constructor  .... Hopefully should never be used
		CommonDetail(const CommonDetail& cd);
		
		// Destructor
		~CommonDetail();
		
		// Access methods for details
		inline const double mass() const{
			return this->_mass;
		}
		inline const double X() const{
			return this->_pos[0];
		}
		inline const double Y() const{
			return this->_pos[1];
		}
		inline const double Z() const{
			return this->_pos[2];
		}
		inline const double Vx() const{
			return this->_vel[0];
		}
		inline const double Vy() const{
			return this->_vel[1];
		}
		inline const double Vz() const{
			return this->_vel[2];
		}

		// Method to print the Common details
		friend std::ostream& operator<<(std::ostream& os,const CommonDetail &detail);
		
		// Equal operator for the class
		CommonDetail& operator=(const CommonDetail& cd);

	private:
		double _mass;
		double* _pos;
		double* _vel;
};

#endif

