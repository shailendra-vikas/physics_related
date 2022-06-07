#ifndef _CHAIN_MESH_H
#define _CHAIN_MESH_H
#define A_VERY_SMALL_NUMBER .0001			//deduce this amount from the coordiantes to make the indices in bound
#include <map>
#include <set>
#include <list>
#include "CommonDetails.h"
#include "Halo.h"
#include <fstream>
/**
	This is header file for implementation of Chaining Mesh.
	In Chaining Mesh the box is devided in my parts. You have to specify how many division per side you want.
	In every such small division we devide it in number of types.
	Every cell contains a IdDetailMap
	IdDetailMap maps id and CommonDetail
	
	It also contains a multimap idCellMap to store id and CellDetailForId.
	CellDetailForId class contains position of cell for that particle and type.
	Hence for a given Id you can get the cell info fast and fetch all other particle in that cell.
	
	The data structure is optimized to get the neighbour particle  in a given distance fast. The most 
	useful method to access data is getNeighbourParticles.
	
	There are two getNeighbourParticles methods. 
	One return the set of ids of the particle.
	Second one return of list of ParticleDetailsForHalo.
	ParticleDetailsForHalo contains some calculated info about the particle like distance which is helpful in 
	making plot as a function of distance.
*/

/**
	This class contains the details of the cell of a particle
*/
class CellDetailForId{
	public:
		//Constructor and destructor
		CellDetailForId(const int& type,const int& nthx,const int& nthy,const int& nthz):_type(type){
		_nth[0]=nthx;
		_nth[1]=nthy;
		_nth[2]=nthz;
		}
		~CellDetailForId(){}

		inline const int& getType(){
			return _type;
		}
		inline const int (&getCord())[3]{
			return _nth;
		}

	private:
		int _type;
		int _nth[3];
};
/**
	This class contains calculated and other details of the particles
*/
class ParticleDetailsForHalo{
	public:
		ParticleDetailsForHalo(const IdType& id,CommonDetail* const &detail,const double& distance):_id(id),_detail(detail),_distance(distance){
		}
		~ParticleDetailsForHalo(){}
		
		
		inline const double& getDistance(){
			return _distance;
		} 
		
		inline CommonDetail* const getDetail(){
			return _detail;
		}

		class compare{
			public:
				bool operator()( ParticleDetailsForHalo* const &detail1,ParticleDetailsForHalo* const &detail2){
					return (*detail1).getDistance() < (*detail2).getDistance();
				}
};
				
	private:
		IdType _id;
		CommonDetail* _detail;
		double _distance;
};

typedef std::map<IdType,CommonDetail*> IdDetailMap;
typedef std::multimap<IdType,CellDetailForId*> IdCellMap;
typedef IdCellMap::const_iterator IdCellConstIt;
typedef std::set<int*> CordinateSet;
typedef std::list<ParticleDetailsForHalo*> PartDetPtrList;

/**
	This is the main class that implements ChainMesh
*/
class ChainMesh{
	public:
		//Public methods
		//Constructors and destructors
		ChainMesh(const unsigned int& dev,const double& sizeBox,const int& types);
		~ChainMesh();

		// Methods to access the Id Detail maps
		void placeParticle(const IdType& id,CommonDetail &detail,const int &type);
		void getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int&type,const bool &cyclic,UIntSet &neighbourPart) const;
		void getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int&type,const bool &cyclic,PartDetPtrList &neighbourPart) const;
		const CommonDetail& getDetail(const IdType &id,const int &type) const;
		
		//Generic methods for particle details
		const bool isParticleInNeighbour(const CommonDetail& p1,const CommonDetail& p2,const double &minDist,const double &maxDist,const bool &cyclic,double &dist) const;
		const double getDistance(const CommonDetail& detail1,const CommonDetail& detail2,const bool &cyclic) const;

		// Generic methods for whole datasets
		unsigned int getTotalNumber(const int type) const;
		// It calculate correlation function .... Ideally this function should not exist in this class
		void createCorrelationFunction(const int &type,const double &r_min,const double &r_max,const int &nBin)const;
		
		//Mathod to get distance for a given number of particle
		double getDistanceForNoParticle(const CommonDetail& detail,const int& N,const int& error,const double& rmin) const;
	private:
		//private methods

		void getCellCord(const IdType &id,const int& type,int& nthx,int& nthy,int& nthz) const;

		void distanceToCellWall(const CommonDetail &detail,double &min,double &max) const;

		void getLocationOfCell(const CommonDetail &detail,int &nthx,int &nthy,int &nthz) const;

		void getNeighbourCells(const CommonDetail &detail,const double &mindist,const double &maxdist,const bool &cyclic,CordinateSet &neighbours) const;

		const IdDetailMap& getCellParticles(const int nthx,const int nthy,const int nthz,const int &type) const;

		void getAllCellNos(CordinateSet& cells) const;
		
		int makeCordCyclic(int i) const;

		void disposeCordinateSet(CordinateSet& cordSet) const;
		//private varibles
		const int _number_of_dev;
		const double _sizeOfCell;
		const int _types;
		const static unsigned int _instatiated;
		IdCellMap _idCellMap;
		IdDetailMap  ****_3DarrayOfMap;
};

#endif
