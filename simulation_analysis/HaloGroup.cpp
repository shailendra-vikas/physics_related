#include "HaloGroup.h"

HaloGroup::HaloGroup(const std::string& name):_filename(std::string(name).c_str()){
	std::ostringstream os;
        os << name;

	//open stream 
	std::ifstream inputFile(os.str().c_str(),std::ifstream::in);

        if(!inputFile.good()){
        	std::cerr << "ERROR: Could nt open file" << std::endl;
	}

	inputFile >> _noOfHalos;
	inputFile.close();
}

HaloGroup::~HaloGroup(){
	HaloMap::iterator halo;
	for(halo=_haloMap.begin();halo!=_haloMap.end();halo++){
		delete halo->second;
	}
}

void HaloGroup::loadHalo(){
	std::ifstream inputFile(_filename.c_str(),std::ifstream::in);
        if(!inputFile.good()){
        	std::cerr << "ERROR: Could nt open file" << std::endl;
	}
	int nos;
	inputFile >> nos;
	unsigned int ID;
	for(int j=0;j<nos;j++){
		Halo &halo = *(new Halo());
		for(int i=0;i<7;i++)
			inputFile >> halo._properties[i];
		inputFile >> ID;
		for(int i=7;i<14;i++)
			inputFile >> halo._properties[i];
			
		HaloMap::iterator it = _haloMap.find(ID);

		_haloMap[ID]=&halo;
	}
	inputFile.close();
	return;
}

const HaloMap& HaloGroup::getHalos() const{
	return _haloMap;
}

void HaloGroup::makeCumulMassPlot() const{
	std::ofstream outputFile("N_vs_m",std::ofstream::out);
	double Mmin=.1,Mmax=1000;
	int bins=40;
	unsigned int counts[bins];
	double massRange[bins+1],m;
	for(int i=0;i<bins;i++){
		massRange[i]=Mmin*pow(10.0,i*log10(Mmax/Mmin)/bins);
		counts[i]=0;
	}	
	massRange[bins]=Mmax;
	
	HaloMap::const_iterator halo;
	for(halo=_haloMap.begin();halo != _haloMap.end();halo++){
		m=(*(halo->second)).getTotalMass();
		if(m < Mmin || m > Mmax)
			std::cerr << " Discarding some values for mass bin m= " << m << std::endl;
		else{
			for(int i=0;i<bins;i++)
				if(m<massRange[i+1]){
					counts[i]++;
					break;
				}
		}
	}
	
	for(int i=0;i<bins;i++)
		outputFile << massRange[i] << " " << massRange[i+1] << " " << counts[i] << std::endl;
		
	outputFile.close();	
}

void HaloGroup::makeCorrelation() const {
	ChainMesh cm(32,37250,1);
	HaloMap::const_iterator halo;
	for(halo=_haloMap.begin();halo != _haloMap.end();halo++){
		const double * cord= (*(halo->second)).getDenCord();
		CommonDetail *detail= new CommonDetail(1,cord[0],cord[1],cord[2]);
		// making fake data so that correlation function can be calulated
		cm.placeParticle(halo->first,(*detail),0);
	}
	cm.createCorrelationFunction(0,10,10000,30);
}

std::ostream& operator<<(std::ostream& os,const HaloGroup& haloGrp){
	os << " Name:" << haloGrp._filename << std::endl;
	os << " Number of halos in file:" << haloGrp._noOfHalos << std::endl;
	return os;
}

void HaloGroup::createOutput() const{
	std::ofstream outputFile("output/Halo_particle.txt",std::ofstream::out|std::ofstream::binary);
	unsigned int halomap_size,particles_size;
	halomap_size = _haloMap.size();
	outputFile.write(reinterpret_cast<char*>(&halomap_size),sizeof(unsigned int));
//	outputFile << _haloMap.size() << std::endl;
	HaloMap::const_iterator halo;
	for(halo=_haloMap.begin();halo != _haloMap.end();halo++){
//		outputFile << " " << halo->first ;
		outputFile.write(reinterpret_cast<const char*>(&halo->first),sizeof(unsigned int));
		for(int i=0;i<TYPES_OF_PARTICLE;i++){
			UIntSet &particleIds =halo->second->_particles[i];
//			outputFile <<" "<< particleIds.size();
			particles_size=particleIds.size();
			outputFile.write(reinterpret_cast<char*>(&particles_size),sizeof(unsigned int));
			UIntSet::iterator particleIt;
			for(particleIt = particleIds.begin(); particleIt !=particleIds.end();particleIt++){
//				outputFile << " " << *particleIt;
				outputFile.write(reinterpret_cast<const char*>(&(*particleIt)),sizeof(unsigned int));
			}	
//			outputFile << " " << particleIds.size();
			outputFile.write(reinterpret_cast<char*>(&particles_size),sizeof(unsigned int));
		}
//		outputFile << std::endl;
	}
}
