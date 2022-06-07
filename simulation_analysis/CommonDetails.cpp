#include "CommonDetails.h"

// Contructor for the most basic details
CommonDetail::CommonDetail(const double  &mass,const double  &xpos,const double &ypos,const double &zpos):_mass(mass),_vel(0){
	_pos=new double[3];
	_pos[0]=xpos;
	_pos[1]=ypos;
	_pos[2]=zpos;
}

//Contructor for the complete details
CommonDetail::CommonDetail(const double  &mass,const double  &xpos,const double &ypos,const double &zpos ,const double  &xvel,const double &yvel,const double &zvel):_mass(mass){
	_pos=new double[3];
	_pos[0]=xpos;
	_pos[1]=ypos;
	_pos[2]=zpos;
	_vel=new double[3];
	_vel[0]=xvel;
	_vel[1]=yvel;
	_vel[2]=zvel;
}

// Constructor  .... Hopefully should never be used
CommonDetail::CommonDetail(const CommonDetail &cd){
	_mass=cd._mass;
	//If positions are mentioned than copy it
	if(cd._pos){
		_pos = new double[3];
		for(int i=0;i<3;i++)
			_pos[i]=cd._pos[i];
	}
	else
		_pos=0;
		
	//If velocity are mentioned than copy it		
	if(cd._vel){
		_vel = new double[3];
		for(int i=0;i<3;i++)
			_vel[i]=cd._vel[i];
	}
	else
		_vel=0;
}

// Method to print the Common details
std::ostream& operator<<(std::ostream& os,const CommonDetail &detail){
	std::cout << "Mass:" << detail._mass << std::endl;
	if(detail._pos)
		std::cout << "Position X:" << detail._pos[0] << " Y:" << detail._pos[1] << " Z:" << detail._pos[2] << std::endl;
	if(detail._vel)
		std::cout << "Velocity Vx:" << detail._vel[0] << " Vy:" << detail._vel[1] << " Vz:" << detail._vel[2] << std::endl;
	return os;
}

// Equal operator for the class
CommonDetail& CommonDetail::operator=(const CommonDetail& cd){
	if(this==&cd)
		return (*this);
	// delete the existing positions ans velocity	
	if(_pos)	
		delete[] _pos;
	if(_vel)
		delete[] _vel;

	_mass=cd._mass;
	//If positions are mentioned than copy it	
	if(cd._pos){
		_pos = new double[3];
		for(int i=0;i<3;i++)
			_pos[i]=cd._pos[i];
	}
	else
		_pos=0;
		
	//If velocity are mentioned than copy it		
	if(cd._vel){
		_vel = new double[3];
		for(int i=0;i<3;i++)
			_vel[i]=cd._vel[i];
	}
	else
		_vel=0;

	return (*this);
}

// Destructor
CommonDetail::~CommonDetail(){
	if(_pos)
		delete[] _pos;
	if(_vel)
		delete[] _vel;
}
