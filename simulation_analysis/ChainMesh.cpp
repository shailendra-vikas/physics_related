#include "ChainMesh.h"
 
const unsigned int ChainMesh::_instatiated=1;


ChainMesh::ChainMesh(const unsigned int& dev,const double& sizeBox,const int& types):_number_of_dev(dev),_sizeOfCell(sizeBox/dev),_types(types){
	_3DarrayOfMap= new IdDetailMap ***[_number_of_dev];

	for(int i=0;i<_number_of_dev;i++)
		_3DarrayOfMap[i]= new IdDetailMap **[_number_of_dev];

	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			_3DarrayOfMap[i][j] = new IdDetailMap*[_number_of_dev];

	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			for(int k=0;k<_number_of_dev;k++)
				_3DarrayOfMap[i][j][k] = new IdDetailMap[_types];
}

ChainMesh::~ChainMesh(){
	//Delete IdCellMap
	IdCellMap::iterator idCellEle;
	for(idCellEle=_idCellMap.begin(); idCellEle != _idCellMap.end(); idCellEle++)
		delete idCellEle->second;
		
	
	// Delete all the CommonDetails before deleting Chain Mesh
	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			for(int k=0;k<_number_of_dev;k++)
				for(int l=0;l<_types;l++){
					IdDetailMap::iterator particle;
					for(particle=_3DarrayOfMap[i][j][k][l].begin();particle!=_3DarrayOfMap[i][j][k][l].end();particle++)
						delete (*particle).second;
				}

	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			for(int k=0;k<_number_of_dev;k++)
				delete [] _3DarrayOfMap[i][j][k];

		
	for(int i=0;i<(int)_number_of_dev;i++)
		for(int j=0;j<(int)_number_of_dev;j++)
			delete [] _3DarrayOfMap[i][j];
	
	for(int i=0;i<(int)_number_of_dev;i++)
		delete [] _3DarrayOfMap[i];

	delete [] _3DarrayOfMap;	
}

void ChainMesh::placeParticle(const IdType &id,CommonDetail &detail,const int &type){
	int nthx,nthy,nthz;
	getLocationOfCell(detail,nthx,nthy,nthz);
	_3DarrayOfMap[nthx][nthy][nthz][type][id]=&detail;
	_idCellMap.insert(std::pair<IdType,CellDetailForId*>(id,new CellDetailForId(type,nthx,nthy,nthz)));
}

void ChainMesh::getCellCord(const IdType &id,const int& type,int& nthx,int& nthy,int& nthz) const{
	std::pair<IdCellConstIt,IdCellConstIt> idCellElements;
	idCellElements =_idCellMap.equal_range(id);
	for(IdCellConstIt idCellElement = idCellElements.first;idCellElement != idCellElements.second;++idCellElement)
		if((*(idCellElement->second)).getType() ==type){
			const int (&cord)[3] =(*(idCellElement->second)).getCord();
			nthx=cord[0];
			nthy=cord[1];
			nthz=cord[2];
			return;
		}	
	std::cerr << " Couldn't find the particle for id:" <<id << " Type:" << type << std::endl;
}

void ChainMesh::distanceToCellWall(const CommonDetail &detail,double &min,double &max) const{
	int nthx,nthy,nthz;
	getLocationOfCell(detail,nthx,nthy,nthz);
	double distToX= detail.X()-nthx*_sizeOfCell;
	distToX=std::min(distToX,_sizeOfCell-distToX);
	double distToY= detail.Y()-nthy*_sizeOfCell;
	distToY=std::min(distToY,_sizeOfCell-distToY);
	double distToZ= detail.Z()-nthz*_sizeOfCell;
	distToZ=std::min(distToZ,_sizeOfCell-distToZ);

	min= std::min(distToX,std::min(distToY,distToZ));
	max= sqrt((_sizeOfCell-distToX)*(_sizeOfCell-distToX)+(_sizeOfCell-distToY)*(_sizeOfCell-distToY)+(_sizeOfCell-distToZ)*(_sizeOfCell-distToZ));
	return;
}

void ChainMesh::getLocationOfCell(const CommonDetail &detail,int &nthx,int &nthy,int &nthz) const{
	nthx=(int)std::floor((detail.X()-A_VERY_SMALL_NUMBER)/_sizeOfCell);
	nthy=(int)std::floor((detail.Y()-A_VERY_SMALL_NUMBER)/_sizeOfCell);
	nthz=(int)std::floor((detail.Z()-A_VERY_SMALL_NUMBER)/_sizeOfCell);
	return;
}

const IdDetailMap& ChainMesh::getCellParticles(const int nthx,const int nthy,const int nthz,const int &type) const{
	return _3DarrayOfMap[nthx][nthy][nthz][type];
}

void ChainMesh::getNeighbourCells(const CommonDetail &detail,const double &mindist,const double &maxdist,const bool& cyclic,CordinateSet &neighbours) const{
	int nthx,nthy,nthz;
	getLocationOfCell(detail,nthx,nthy,nthz);

	double closestDistance,furthestDistance;
	distanceToCellWall(detail,closestDistance,furthestDistance);


	int maxCellDist= (int)std::ceil((maxdist-closestDistance)/_sizeOfCell);
	int minCellDist= (int)std::floor((mindist-furthestDistance)/(1.733*_sizeOfCell));

	for(int i=nthx-maxCellDist;i<= nthx+maxCellDist;i++)
		for(int j=nthy-maxCellDist;j<= nthy+maxCellDist;j++)
			for(int k=nthz-maxCellDist;k<= nthz+maxCellDist;k++){
				if(	minCellDist >= 0 && 
					i >= nthx-minCellDist+1 && 
					i <= nthx+minCellDist-1 &&
					j >= nthy-minCellDist+1 && 
					j <= nthy+minCellDist-1 &&
					k >= nthz-minCellDist+1 && 
					k <= nthz+minCellDist-1 
				)
					continue;
				if(cyclic){
					int * cord= new int[3];
					cord[0]=makeCordCyclic(i);
					cord[1]=makeCordCyclic(j);
					cord[2]=makeCordCyclic(k);
					neighbours.insert(cord);
				}	
				else
					if(i>=0 && i<(int)_number_of_dev && j>=0 && j<(int)_number_of_dev && k>=0 && k<(int)_number_of_dev){
						int * cord= new int[3];
						cord[0]=i;
						cord[1]=j;
						cord[2]=k;
						neighbours.insert(cord);
					}	
			}
};


void ChainMesh::getAllCellNos(CordinateSet& cells) const{
	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			for(int k=0;k<_number_of_dev;k++){
				int * cord= new int[3];
				cord[0]=i;
				cord[1]=j;
				cord[2]=k;
				cells.insert(cord);
			}	
}

void ChainMesh::getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int &type,const bool &cyclic,UIntSet &neighbourPart)const {
	CordinateSet neighbourCells;
	getNeighbourCells(detail,minDist,maxDist,cyclic,neighbourCells);
	CordinateSet::const_iterator cell;
	double dist;
	for(cell= neighbourCells.begin();cell != neighbourCells.end();cell++){
		const IdDetailMap &cellParticles = getCellParticles((*cell)[0],(*cell)[1],(*cell)[2],type);		
		IdDetailMap::const_iterator cellParticle;

		for(cellParticle=cellParticles.begin();cellParticle!=cellParticles.end();cellParticle++){
			if(isParticleInNeighbour(detail,*(cellParticle->second),minDist,maxDist,cyclic,dist))
				neighbourPart.insert(cellParticle->first);	
		}
	}
	disposeCordinateSet(neighbourCells);
}

void ChainMesh::getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int &type,const bool &cyclic,PartDetPtrList &neighbourPart)const {
	CordinateSet neighbourCells;
	getNeighbourCells(detail,minDist,maxDist,cyclic,neighbourCells);
	CordinateSet::const_iterator cell;
	double dist;
	for(cell= neighbourCells.begin();cell != neighbourCells.end();cell++){
		const IdDetailMap &cellParticles = getCellParticles((*cell)[0],(*cell)[1],(*cell)[2],type);		
		IdDetailMap::const_iterator cellParticle;

		for(cellParticle=cellParticles.begin();cellParticle!=cellParticles.end();cellParticle++){
			if(isParticleInNeighbour(detail,*(cellParticle->second),minDist,maxDist,cyclic,dist))
				neighbourPart.push_back(new ParticleDetailsForHalo(cellParticle->first,cellParticle->second,dist));	
		}
	}
	disposeCordinateSet(neighbourCells);
}

double ChainMesh::getDistanceForNoParticle(const CommonDetail& detail,const int& N,const int& error,const double& rmin) const {
	double max;
	double min;
	int count=0;

	max= (rmin < TINY)?0.1:rmin;
	// Get the range in which rmax would lie
	do{
		count=0;
		UIntSet neighbourPart;
		min=max;
		max=max*2;
		for(int i=0;i< _types;i++)
			getNeighbourParticles(detail,rmin,max,i,1,neighbourPart);
		
		count= neighbourPart.size();
		if(min < TINY && count >N)
			std::cerr << " The rmax value is wrong" << std::endl;
	}while(count < N);

	// Now detail research for the rmax to have number of particle N
	do{
		UIntSet neighbourPart;
		count=0;

		for(int i=0;i< _types;i++)
			getNeighbourParticles(detail,rmin,(max+min)/2,i,1,neighbourPart);
		
		count= neighbourPart.size();
		
		if(count < N-error )
			min =(max+min)/2;
		if(count > N+error )			
			max= (max+min)/2;
		
	}while(count < N-error && count > N+error);
	return (max+min)/2;
}

const bool ChainMesh::isParticleInNeighbour(const CommonDetail& detail1,const CommonDetail& detail2,const double &minDist,const double &maxDist, const bool &cyclic,double &dist) const{
	dist = getDistance(detail1,detail2,cyclic);
	return (dist<maxDist && dist >minDist); 
}

const double ChainMesh::getDistance(const CommonDetail& detail1,const CommonDetail& detail2,const bool &cyclic) const{

	double dist=0;
	if(!cyclic){
		dist = std::sqrt(std::pow(detail1.X()-detail2.X(),2)+std::pow(detail1.Y()-detail2.Y(),2)+std::pow(detail1.Z()-detail2.Z(),2));
		return dist;
	}
	double temp;
	double boxSize=_sizeOfCell*_number_of_dev;
	temp=fabs(detail1.X()-detail2.X());
	dist+=(temp>(boxSize-temp))?(boxSize-temp)*(boxSize-temp):temp*temp;
	temp=fabs(detail1.Y()-detail2.Y());
	dist+=(temp>(boxSize-temp))?(boxSize-temp)*(boxSize-temp):temp*temp;
	temp=fabs(detail1.Z()-detail2.Z());
	dist+=(temp>(boxSize-temp))?(boxSize-temp)*(boxSize-temp):temp*temp;
	dist = std::sqrt(dist);
	return dist;
}

unsigned int ChainMesh::getTotalNumber(const int type) const{
	unsigned int nTotal=0;
	for(int i=0;i<_number_of_dev;i++)
		for(int j=0;j<_number_of_dev;j++)
			for(int k=0;k<_number_of_dev;k++){
				const IdDetailMap &particles = getCellParticles(i,j,k,type);
				nTotal += particles.size();
			}	
	return nTotal;
}

void ChainMesh::createCorrelationFunction(const int &type,const double &r_min,const double &r_max,const int &nBin)const {
	std::ofstream outputFile("Correlation2.txt",std::ofstream::out);
	unsigned int nTotal = getTotalNumber(type);
	double boxSize=_sizeOfCell*_number_of_dev;
	for(double r=r_min;r<r_max;r*=std::pow(10.0,std::log10(r_max/r_min)/(double)nBin)){
		unsigned int pair=0;
		CordinateSet allCells;
		getAllCellNos(allCells);
		CordinateSet::const_iterator cell;
		double r2=r*std::pow(10.0,std::log10(r_max/r_min)/(double)nBin);

		for(cell=allCells.begin();cell!=allCells.end();cell++){
			const IdDetailMap& particles = getCellParticles((*cell)[0],(*cell)[1],(*cell)[2],type);
			IdDetailMap::const_iterator particle;
			for(particle=particles.begin();particle != particles.end();particle++){
				UIntSet neighbours;
				getNeighbourParticles(*(particle->second),r,r2,type,(bool)1,neighbours);
				pair += neighbours.size();
			}
		}
		
		disposeCordinateSet(allCells);
		double dv2=4*3.14*(r2*r2*r2-r*r*r)/3;
		double onePlusSai= pair*std::pow((boxSize/nTotal),2)*boxSize/dv2;
		outputFile <<  r << "     " << r2 << "     " << onePlusSai << std::endl; 
		std::cout  <<  r << "     " << r2 << "     " << onePlusSai << std::endl; 
		std::cout << "Total:" << nTotal << std::endl;
	}
	outputFile.close();
}

const CommonDetail&  ChainMesh::getDetail(const IdType &id,const int &type) const{
	int nthx,nthy,nthz;
	getCellCord(id,type,nthx,nthy,nthz);
	IdDetailMap::const_iterator element= _3DarrayOfMap[nthx][nthy][nthz][type].find(id);
	if(element == _3DarrayOfMap[nthx][nthy][nthz][type].end())
		std::cerr << "Error: Could not find Particle:" << id << " in Mesh" << std::endl;
	return *(element->second);
}

int ChainMesh::makeCordCyclic(int i) const{
	if(i < 0)
		return i+_number_of_dev;
	else 
		if(i > _number_of_dev-1)
			return i-_number_of_dev ;
		else
			return i;	
}

void ChainMesh::disposeCordinateSet(CordinateSet& cordSet) const{
	CordinateSet::iterator cord;
	for(cord=cordSet.begin();cord != cordSet.end();cord++)
		delete[] (*cord);
}
