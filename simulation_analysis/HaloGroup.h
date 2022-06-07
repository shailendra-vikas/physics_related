#ifndef _HALOGROUP_H
#define _HALOGROUP_H

#include <fstream>
#include <sstream>
#include <map>
#include "Halo.h"
#include "ChainMesh.h"


typedef std::map<unsigned int,Halo*> HaloMap;
/**
	This class is for various method which can act on the Halo and header info for the HaloGrp files	
*/
class HaloGroup{
	public:
		HaloGroup(const std::string& name);
		~HaloGroup();

		void loadHalo();

		const HaloMap &getHalos() const;
		
		void makeCumulMassPlot() const;

		friend std::ostream& operator<<(std::ostream& os,const HaloGroup& haloGrp);
		
		void makeCorrelation() const;
		
		void createOutput() const;

	private:
		const std::string _filename;
		int _noOfHalos;
		HaloMap _haloMap;
};
#endif
