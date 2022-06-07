#include "Halo.h"

Halo::Halo(){
}

Halo::~Halo(){
}

const double Halo::getVirialRadius() const{
	//do the calculation
	return 1000;
}

// Print the Halo
std::ostream& operator<<(std::ostream& os,const Halo &halo){
	os << "MASS_TOTAL:" << halo._properties[MASS_TOTAL] 
		<< " X_GRP:" << halo._properties[X_GRP] 
		<< " Y_GRP:" << halo._properties[Y_GRP] 
		<< " Z_GRP:" << halo._properties[Z_GRP] 
		<< " X_DEN:" << halo._properties[X_DEN] 
		<< " Y_DEN:" << halo._properties[Y_DEN] 
		<< " Z_DEN:" << halo._properties[Z_DEN] 
		<< " MASS_B:" << halo._properties[MASS_B] 
		<< " MASS_GAS:" << halo._properties[MASS_GAS] 
		<< " MASS_COLD_GAS:" << halo._properties[MASS_COLD_GAS] 
		<< " V_STAR:" << halo._properties[V_STAR] 
		<< " Z_STAR:" << halo._properties[Z_STAR] 
		<< " Z_GAS:" << halo._properties[Z_GAS] 
		<< " SFR:" << halo._properties[SFR] 
		<< std::endl;
	return os;
}
