#include "SimulationData.h"

using namespace std;

/*
    Load the all particle type stars, dark matter and gas.
*/
void loadGasHaloStar(SimulationData& simulation){
	simulation.loadGasWoVel();
	simulation.loadHaloWoVel();
	simulation.loadStarsWoVel();

}

/*
    Load the Halo groups information from the simulation
*/
void loadHalo(SimulationData& simulation){
	simulation.addHaloGroup("../D4/GRP/grpcens_D4_033");
	simulation.loadHaloGroups();
}

/*
    Make the connection between the simulation particle to the halo catalog.
*/
void makeFileOfHalosWithParticle(SimulationData& simulation){
	simulation.linkParticlesToHaloInFile();
}

int main(int argc,char** argv){
	SimulationData simulation("../D4/DAT/snap_D4_033",4);
	loadGasHaloStar(simulation);
	loadHalo(simulation);
	makeFileOfHalosWithParticle(simulation);

	simulation.createCorrelationFunction(1,10,10000,30);
	simulation.generateHaloProfile();
	simulation.generateHaloProfileForEqNo();


/* This is for the halo correlation and halo cumulative plot

	HaloGroup grp("../D4/GRP/grpcens_D4_044");
	grp.loadHalo();
	grp.makeCumulMassPlot();
	grp.makeCorrelation();

*/

return 0;
}

