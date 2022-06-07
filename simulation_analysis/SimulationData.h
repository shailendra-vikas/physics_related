#ifndef _SIMULATIONDATA_H
#define _SIMULATIONDATA_H
//#define DOWNSAMPLE
#include "HeaderData.h"
#include "ChainMesh.h"
#include "HaloGroup.h"
#include <fstream>
#include <sstream>
#include <ctime>
#include <cstdlib>

struct headerBlock {
  unsigned int      npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[6];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
};

typedef struct particle_data {
       int id;
  double  Pos[3];
  double  Vel[3];
  double  Mass, Rho, Temp, Ne, Hsml;
} particle_data;

/*
    Represent the Simulation Data. Various method defined are the typical queries that the analysis would be doing. The name of the methods are self explainorties.
    The simulation assumes the periodic boundary condition. Estimate distance assuming the periodic boundary distance.

    The programe try to devide the simulation in small cubes so that calculating neighbour particle only consider particle from only the palussible cubes. 
    The periodic boundary condition make the estimation of possible cubes non trivial task. The datastructure is transparent to the user.
*/

class SimulationData{
	public:
		SimulationData(const std::string& name,const unsigned int & filecount);
		~SimulationData();

		inline bool isGasPresent() const{
			return _nTotalPart[0];			
		}
		void loadGasWoVel();
		void loadGasWVel();

		inline bool isHaloPresent() const{
			return _nTotalPart[1];			
		}
		void loadHaloWoVel();
		void loadHaloWVel();

		inline bool isDiskPresent() const{
			return _nTotalPart[2];			
		}
		void loadDiskWoVel();
		void loadDiskWVel();

		inline bool isBulgePresent() const{
			return _nTotalPart[3];			
		}
		void loadBulgeWoVel();
		void loadBulgeWVel();

		inline bool isStarPresent() const{
			return _nTotalPart[4];			
		}
		void loadStarsWoVel();
		void loadStarsWVel();

		inline bool isBndryPresent() const{
			return _nTotalPart[5];			
		}
		void loadBndryWoVel();
		void loadBndryWVel();

		friend std::ostream& operator<<(std::ostream& os,const SimulationData& simData);

		void getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int&type,const bool &cyclic,UIntSet& cells) const;

		void getNeighbourParticles(const CommonDetail& detail,const double & minDist,const double & maxDist,const int &type,const bool &cyclic,PartDetPtrList& neighbourPart) const;
		
		void createCorrelationFunction(const int &type,const double &r_min,const double &r_max,const int &nBin)const;
		
		void addHaloGroup(const std::string& name);
			
		void loadHaloGroups();
		
		void linkParticlesToHalo();

		void linkParticlesToHaloInFile();
		
		void generateHaloProfile() const;
		
		void generateHaloProfileForEqNo() const;
		
		inline HaloGroup* getHaloGroup() const{
			return _haloGroup;
		}

	private:
		HeaderData **_headerList;
		const std::string _filename;
		unsigned int _filecount;
		bool _gasVel,_haloVel,_diskVel,_bulgeVel,_starsVel,_bndryVel;
		double _time,_redshift;
		double _boxSize;
		double _omega0,_omegaLambda,_hubbleParam;
		unsigned int _nTotalPart[TYPES_OF_PARTICLE];
		int _flag_sfr,_flag_feedback,_flag_cooling;
		ChainMesh *_ptr_mesh;
		HaloGroup *_haloGroup;
		void loadParticles(const int &type,const bool &vel);
};


#endif
