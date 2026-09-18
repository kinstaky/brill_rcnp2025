#pragma once

#include <TTree.h>

namespace brill {

struct GaggEvent {
	int num;
	int index[64];
	double energy[64], integral[64], time[64];
	bool saturation[64];
	double camplitude[64];
	// free-tau fitting
	double tau[64], base[64], amplitude[64], t0[64];
	// fix-tau fitting
	double ftau[64], fbase[64], famplitude[64], ft0[64];
	// meta
	int run, entry;
};

void SetupInput(TTree *tree, GaggEvent &event, const std::string &prefix = "");
void SetupOutput(TTree *tree, GaggEvent &event);

}