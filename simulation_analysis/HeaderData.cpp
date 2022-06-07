#include "HeaderData.h"

//Constructor
HeaderData::HeaderData(const std::string& name):_filename(name){
}
//Destructor
HeaderData::~HeaderData(){
	if(_mass)
		delete[] _mass;
	if(_npart)
		delete[] _npart;
}

//If there is any mass section present or all the masses are in header
bool HeaderData::isMassSectionPresent() const{
	std::numeric_limits<double> double_info;
	for(int i=0;i<TYPES_OF_PARTICLE;i++)
		if(_npart[i]!= 0 && _mass[i] > double_info.epsilon()) return 1;
	return 0;
}

//Setter methods for various attribute (no of particles for types ,mass)
void HeaderData::setNumberParticles(unsigned int*& npart){
		_npart=npart;
}
void HeaderData::setMass(double*& mass){
		_mass=mass;
}

// Method for printing 
std::ostream& operator<<(std::ostream& os,const HeaderData &headerData){
	os<< "Header info for " << headerData._filename << std::endl;
	for(int i=0;i<TYPES_OF_PARTICLE;i++) {
		os << "For Type:" << i << " Mass:" << headerData._mass[i] << " Numbers:" << headerData._npart[i] << std::endl;
	}
	return os;
}

		

