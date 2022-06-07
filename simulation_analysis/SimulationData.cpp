#include "SimulationData.h"

SimulationData::SimulationData(const std::string& name,unsigned int const & filecount):_filename(std::string(name).c_str()),_filecount(filecount),_haloGroup(0){
	_headerList=new HeaderData*[filecount];
	for(int i=0;i<(int)_filecount;i++){
                std::ostringstream os;
                os << name;
                if(_filecount !=1)
                        os << "." <<i;

                //open stream 
                std::ifstream inputFile(os.str().c_str(),std::ifstream::in|std::ifstream::binary);

                if(!inputFile.good()){
                        std::cerr << "ERROR: Could nt open file" << std::endl;
                }

		int sizeStart,sizeEnd;
		headerBlock header;

	        inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
        	inputFile.read(reinterpret_cast<char*>(&header) ,sizeof(headerBlock));
	        inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
		if(sizeStart != sizeEnd)
                	std::cerr << "ERROR:Size of the Block at the end does not match with size at starting" << std::endl;

		_headerList[i]= new HeaderData(os.str().c_str());

		//Load Simulation Specific Fields
		if(i==0){
			_time=header.time;
			_redshift=header.redshift;
			_boxSize=header.BoxSize;
			_omega0=header.Omega0;
			_omegaLambda=header.OmegaLambda;
			_hubbleParam=header.HubbleParam;
			_flag_sfr=header.flag_sfr;
			_flag_feedback=header.flag_feedback;
			_flag_cooling=header.flag_cooling;
			for(int j=0;j<TYPES_OF_PARTICLE;j++)
				_nTotalPart[j]=header.npartTotal[j];
		}

		//Create Headers
		unsigned int *npart=new unsigned int[TYPES_OF_PARTICLE];
		double * mass = new double[TYPES_OF_PARTICLE];
		for(int j=0;j<TYPES_OF_PARTICLE;j++){
			npart[j]=header.npart[j];
			mass[j]=header.mass[j];
		}

		(*_headerList[i]).setMass(mass);
		(*_headerList[i]).setNumberParticles(npart);

		inputFile.close();
	}
	_ptr_mesh= new ChainMesh(128,_boxSize,TYPES_OF_PARTICLE);
}

SimulationData::~SimulationData(){
	if(_headerList){
		for(int j=0;j<(int)_filecount;j++)
			delete _headerList[j];
		delete[] _headerList;
	}
	if(_ptr_mesh)
		delete _ptr_mesh;
		
	if(_haloGroup)
		delete _haloGroup; 
}

void SimulationData::loadGasWoVel(){
	_gasVel=0;
	loadParticles(0,_gasVel);
}
void SimulationData::loadGasWVel(){
	_gasVel=1;
	loadParticles(0,_gasVel);
}

void SimulationData::loadHaloWoVel(){
	_haloVel=0;
	loadParticles(1,_haloVel);
}
void SimulationData::loadHaloWVel(){
	_haloVel=1;
	loadParticles(1,_haloVel);
}


void SimulationData::loadDiskWoVel(){
	_diskVel=0;
	loadParticles(2,_diskVel);
}
void SimulationData::loadDiskWVel(){
	_diskVel=1;
	loadParticles(2,_diskVel);
}

void SimulationData::loadBulgeWoVel(){
	_bulgeVel=0;
	loadParticles(3,_bulgeVel);
}
void SimulationData::loadBulgeWVel(){
	_bulgeVel=1;
	loadParticles(3,_bulgeVel);
}


void SimulationData::loadStarsWoVel(){
	_starsVel=0;
	loadParticles(4,_starsVel);
}
void SimulationData::loadStarsWVel(){
	_starsVel=1;
	loadParticles(4,_starsVel);
}

void SimulationData::loadBndryWoVel(){
	_bndryVel=0;
	loadParticles(5,_bndryVel);
}
void SimulationData::loadBndryWVel(){
	_bndryVel=1;
	loadParticles(5,_bndryVel);
}

std::ostream& operator<<(std::ostream& os,const SimulationData& simData){
	os << " Info for simulation " << simData._filename << std::endl;
	os << " Number of data files:" << simData._filecount <<  std::endl;
	os << " Time:"<< simData._time << " RedShift:" << simData._redshift << std::endl;
	os << " BoxSize:" << simData._boxSize << std::endl;
	os << " Paramters are Omeaga0:" << simData._omega0 << " OmeagaLambda:" << simData._omegaLambda << " Hubble Constant:" << simData._hubbleParam << std::endl;
	for(int i=0;i<TYPES_OF_PARTICLE;i++){
		os << " For Type:" << i << " Number of Partical:" << simData._nTotalPart[i] << std::endl;
	}
	os << "Flags are SFR:" << simData._flag_sfr << " FEEDBACK:" << simData._flag_feedback << " COOLING:" << simData._flag_cooling << std::endl;
	return os;
};

void SimulationData::loadParticles(const int &type,const bool &vel){
	for(int i=0;i<(int)_filecount;i++){
		std::ifstream  inputFile((*_headerList[i]).getName().c_str(),std::ifstream::in|std::ifstream::binary);		
	   	int sizeStart,sizeEnd,sizeoffloat=sizeof(float),sizeofint=sizeof(int),intdata;
		float floatdata;
		unsigned int numberOfPart =(*_headerList[i]).getNumberForType(type);
		particle_data *particles =  new particle_data[numberOfPart];

		//skip header
		{
	                headerBlock header;
        		inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
	                inputFile.read(reinterpret_cast<char*>(&header) ,sizeof(headerBlock));
        		inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
		}
		//Read Position	
       		inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
		for(int j=0; j<type;j++)
			inputFile.ignore(3*sizeoffloat*(*_headerList[i]).getNumberForType(j));
		for(unsigned int j=0;j< numberOfPart;j++){
	                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
        	        particles[j].Pos[0]=floatdata;
	                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
        	        particles[j].Pos[1]=floatdata;
	                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
        	        particles[j].Pos[2]=floatdata;
		}
		for(int j=type+1; j<TYPES_OF_PARTICLE;j++)
			inputFile.ignore(3*sizeoffloat*(*_headerList[i]).getNumberForType(j));
        	inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
		if(sizeStart != sizeEnd)
                	std::cerr << "ERROR:Size of the Block at the end does not match with size at starting for Reading Position" << std::endl;
		//Read Velocity	
       		inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
		if(vel){
			for(int j=0; j<type;j++)
				inputFile.ignore(3*sizeoffloat*(*_headerList[i]).getNumberForType(j));
			for(unsigned int j=0;j< numberOfPart;j++){
		                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
	        	        particles[j].Vel[0]=floatdata;
		                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
	        	        particles[j].Vel[1]=floatdata;
		                inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
	        	        particles[j].Vel[2]=floatdata;
			}
			for(int j=type+1; j<TYPES_OF_PARTICLE;j++)
				inputFile.ignore(3*sizeoffloat*(*_headerList[i]).getNumberForType(j));

		}
		else{
			inputFile.ignore(sizeStart);
		}
        	inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
		if(sizeStart != sizeEnd)
                	std::cerr << "ERROR:Size of the Block at the end does not match with size at starting for Reading Velocity" << std::endl;
		//Read id	
       		inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
		for(int j=0; j<type;j++)
			inputFile.ignore(sizeofint*(*_headerList[i]).getNumberForType(j));
		for(unsigned int j=0;j< numberOfPart;j++){
	                inputFile.read(reinterpret_cast<char*>(&intdata),sizeofint);
        	        particles[j].id=intdata;
		}
		for(int j=type+1; j<TYPES_OF_PARTICLE;j++)
			inputFile.ignore(sizeofint*(*_headerList[i]).getNumberForType(j));
        	inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
		if(sizeStart != sizeEnd)
                	std::cerr << "ERROR:Size of the Block at the end does not match with size at starting for Reading Id " << std::endl;
		//Read Position	
		//Read mass	
		if((*_headerList[i]).isMassSectionPresent()){
       			inputFile.read(reinterpret_cast<char*>(&sizeStart) ,sizeof(int));
			for(int k=0;k<type;k++)
				if((*_headerList[i]).getMassForType(k) < std::numeric_limits<double>::min())
					inputFile.ignore((*_headerList[i]).getNumberForType(k)*sizeoffloat);

			if((*_headerList[i]).getMassForType(type) < std::numeric_limits<double>::min())
				for(unsigned int j=0;j< numberOfPart;j++){
	                		inputFile.read(reinterpret_cast<char*>(&floatdata),sizeoffloat);
		        	        particles[j].Mass=floatdata;
				}
			else
				for(unsigned int j=0;j< numberOfPart;j++)
		        	        particles[j].Mass=(*_headerList[i]).getMassForType(type);

			for(int k=type+1;k<TYPES_OF_PARTICLE;k++)
				if((*_headerList[i]).getMassForType(k) < std::numeric_limits<double>::min())
					inputFile.ignore((*_headerList[i]).getNumberForType(k)*sizeoffloat);

        		inputFile.read(reinterpret_cast<char*>(&sizeEnd) ,sizeof(int));
			if(sizeStart != sizeEnd)
                		std::cerr << "ERROR:Size of the Block at the end does not match with size at starting for Reading Mass" << std::endl;
		}
		inputFile.close();
		std::srand((unsigned int)time(0));
		for(unsigned int j=0;j< numberOfPart;j++){
		#ifdef DOWNSAMPLE
			if(((double)rand())/(RAND_MAX+1.0) < 0.005){
		#else
			if(1){
		#endif
				CommonDetail *detail;
				if(vel)
					detail=new CommonDetail(particles[j].Mass,particles[j].Pos[0],particles[j].Pos[1],particles[j].Pos[2],particles[j].Vel[0],particles[j].Vel[1],particles[j].Vel[2]);
				else
					detail=new CommonDetail(particles[j].Mass,particles[j].Pos[0],particles[j].Pos[1],particles[j].Pos[2]);
				(*_ptr_mesh).placeParticle(particles[j].id,(*detail),type);
			}
		}
		delete[] particles;
	}
}

void SimulationData::getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int &type,const bool &cyclic,UIntSet& neighbourPart) const{
	return (*_ptr_mesh).getNeighbourParticles(detail,minDist,maxDist,type,cyclic,neighbourPart);

}

void SimulationData::getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int &type,const bool &cyclic,PartDetPtrList& neighbourPart) const{
	return (*_ptr_mesh).getNeighbourParticles(detail,minDist,maxDist,type,cyclic,neighbourPart);

}

void SimulationData::createCorrelationFunction(const int &type,const double &r_min,const double &r_max,const int &nBin) const{
	(*_ptr_mesh).createCorrelationFunction(type,r_min,r_max,nBin);
}



void SimulationData::addHaloGroup(const std::string& name){
	if(_haloGroup)
		delete _haloGroup;
	_haloGroup = new HaloGroup(name);
}

void SimulationData::loadHaloGroups(){
	if(!_haloGroup){
		std::cerr << " Error: add Halo group first " << std::endl;
		return;
	}
	(*_haloGroup).loadHalo();
}

void SimulationData::linkParticlesToHalo(){
	if(!_haloGroup){
		std::cerr << " Error: add Halo group first " << std::endl;
		return;
	}

	HaloMap::const_iterator halo;
	const HaloMap &haloMap =(*_haloGroup).getHalos();
	std::cout << " The number of halos:" << haloMap.size() << std::endl;
	for(halo=haloMap.begin();halo != haloMap.end() ; halo++){
		const double *cord = (*(halo->second)).getDenCord();
		CommonDetail& fakeParticle = *(new CommonDetail(0,cord[0],cord[1],cord[2]));
		for(int j=0;j<TYPES_OF_PARTICLE;j++){
			UIntSet neighbours;
			getNeighbourParticles(fakeParticle,0,(*(halo->second)).getVirialRadius(),j,(bool)1,neighbours);
			UIntSet::const_iterator id;
			for(id = neighbours.begin();id != neighbours.end();id++)
				(halo->second)->insertParticle(*id,j);	
		}
		delete &fakeParticle;
	}
}


void SimulationData::linkParticlesToHaloInFile(){
	if(!_haloGroup){
		std::cerr << " Error: add Halo group first " << std::endl;
		return;
	}
	std::ofstream outputFile("output/Halo_particle.txt",std::ofstream::out|std::ofstream::binary);
	unsigned int halomap_size,particles_size;
	HaloMap::const_iterator halo;

	const HaloMap &haloMap =(*_haloGroup).getHalos();
	halomap_size = haloMap.size();
	outputFile.write(reinterpret_cast<char*>(&halomap_size),sizeof(unsigned int));
	for(halo=haloMap.begin();halo != haloMap.end() ; halo++){
		const double *cord = (*(halo->second)).getDenCord();
		CommonDetail& fakeParticle = *(new CommonDetail(0,cord[0],cord[1],cord[2]));
		
		outputFile.write(reinterpret_cast<const char*>(&halo->first),sizeof(unsigned int));
		for(int j=0;j<TYPES_OF_PARTICLE;j++){
			UIntSet potential_neighbours,neighbours;
			getNeighbourParticles(fakeParticle,0,(*(halo->second)).getVirialRadius(),j,(bool)1,potential_neighbours);
			UIntSet::const_iterator id;
			for(id = potential_neighbours.begin();id != potential_neighbours.end();id++)
				neighbours.insert(*id);
			particles_size=neighbours.size();
			outputFile.write(reinterpret_cast<char*>(&particles_size),sizeof(unsigned int));
			for(id = neighbours.begin();id != neighbours.end();id++)
				outputFile.write(reinterpret_cast<const char*>(&(*id)),sizeof(unsigned int));			
			outputFile.write(reinterpret_cast<char*>(&particles_size),sizeof(unsigned int));				
		}
		delete &fakeParticle;
	}
	outputFile.close();
}



void SimulationData::generateHaloProfile() const{
	if(!_haloGroup){
		std::cerr << " Error: add Halo group first " << std::endl;
		return;
	}
	
	double Rmin=1,Rmax=1000;
	int bins=30;
	double density[bins];
	double range[bins+1];
	for(int i=0;i<bins;i++){
		range[i]=Rmin*pow(10.0,i*log10(Rmax/Rmin)/bins);
		density[i]=0;
	}	
	range[bins]=Rmax;
	std::ofstream outputFile("HaloProfile.txt",std::ofstream::out);
	HaloMap::const_iterator halo;
	const HaloMap &haloMap =(*_haloGroup).getHalos();
	int countHalo;
	for(int k=0;k<bins;k++){
		countHalo=0;
		for(halo=haloMap.begin();halo != haloMap.end() ; halo++){
			if((*(halo->second)).getTotalMass() < 100)
				continue;
			countHalo++;
			const double *cord = (*(halo->second)).getDenCord();
			CommonDetail& fakeParticle = *(new CommonDetail(0,cord[0],cord[1],cord[2]));

			double shellMassforHalo=0;
			for(int j=0;j<TYPES_OF_PARTICLE;j++){
				UIntSet neighbours;
				getNeighbourParticles(fakeParticle,range[k],range[k+1],j,(bool)1,neighbours);
				UIntSet::const_iterator id;
				for(id = neighbours.begin();id != neighbours.end();id++){
					const CommonDetail& detail=(*_ptr_mesh).getDetail(*id,j);
					shellMassforHalo += detail.mass();
				}
				density[k]+=shellMassforHalo/(4.18879*(pow(range[k+1],3)-pow(range[k],3)));
			}	
		}
		density[k]/= countHalo;
		outputFile << range[k] << " "<< range[k+1] << " " << density[k] << std::endl;
		std::cout << range[k] << " "<< range[k+1] << " " << density[k] << " " << countHalo << std::endl;
	}
	
	outputFile.close();	
}

void SimulationData::generateHaloProfileForEqNo() const{
	int N=1000;
	int bins=40;
	double density[bins];
	std::ofstream outputFile("HaloProfileForNo.txt",std::ofstream::out);	
	int countHalo =0;
	HaloMap::const_iterator halo;
	const HaloMap &haloMap =(*_haloGroup).getHalos();
	for(halo=haloMap.begin();halo != haloMap.end() ; halo++){
		if((*(halo->second)).getTotalMass() < 100)
			continue;
		countHalo++;	
		const double *cord = (*(halo->second)).getDenCord();
		CommonDetail& fakeParticle = *(new CommonDetail(0,cord[0],cord[1],cord[2]));
		PartDetPtrList neighbours;
		for(int j=0;j<TYPES_OF_PARTICLE;j++)
			getNeighbourParticles(fakeParticle,0,(*(halo->second)).getVirialRadius(),j,(bool)1,neighbours);
		neighbours.sort(ParticleDetailsForHalo::compare());
		
		PartDetPtrList::iterator current= neighbours.begin();
		double mindist=0,maxdist;
		for(int i=0;i<bins;i++){
			int particlecount=0;
			double massaccum =0;
			while(current != neighbours.end() && particlecount <= N){
				massaccum += (*(*current)->getDetail()).mass();
				maxdist = (*current)->getDistance();
				particlecount++;
				current++;
			}
			density[i]+=massaccum/(4.18879*(pow(maxdist,3)-pow(mindist,3)));
			mindist=maxdist;
		}
		
		for(current=neighbours.begin(); current != neighbours.end(); current++)
			delete *current;
	}
	for(int i=0;i<bins;i++){
		density[i] /= countHalo;
		outputFile << density[i] << " " << i*N << std::endl; 
	}	
	outputFile.close();	
}

