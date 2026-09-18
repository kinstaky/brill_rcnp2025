#include "include/event/gagg_event.h"

namespace brill {

void SetupInput(
	TTree *tree,
	GaggEvent &event,
	const std::string &prefix
) {
	tree->SetBranchAddress((prefix+"n").c_str(), &event.num);
	tree->SetBranchAddress((prefix+"index").c_str(), event.index);
	tree->SetBranchAddress((prefix+"energy").c_str(), event.energy);
	tree->SetBranchAddress((prefix+"integral").c_str(), event.integral);
	tree->SetBranchAddress((prefix+"time").c_str(), event.time);
	tree->SetBranchAddress((prefix+"saturation").c_str(), event.saturation);
	tree->SetBranchAddress((prefix+"camp").c_str(), event.camplitude);
	// free-tau fitting
	tree->SetBranchAddress((prefix+"tau").c_str(), event.tau);
	tree->SetBranchAddress((prefix+"base").c_str(), event.base);
	tree->SetBranchAddress((prefix+"amp").c_str(), event.amplitude);
	tree->SetBranchAddress((prefix+"t0").c_str(), event.t0);
	// fix-tau fitting
	tree->SetBranchAddress((prefix+"ftau").c_str(), event.ftau);
	tree->SetBranchAddress((prefix+"fbase").c_str(), event.fbase);
	tree->SetBranchAddress((prefix+"famp").c_str(), event.famplitude);
	tree->SetBranchAddress((prefix+"ft0").c_str(), event.ft0);
	// meta
	tree->SetBranchAddress((prefix+"run").c_str(), &event.run);
	tree->SetBranchAddress((prefix+"entry").c_str(), &event.entry);
}


void SetupOutput(
	TTree *tree,
	GaggEvent &event
) {
	tree->Branch("n", &event.num, "n/I");
	tree->Branch("index", event.index, "index[n]/I");
	tree->Branch("energy", event.energy, "energy[n]/D");
	tree->Branch("integral", event.integral, "integral[n]/D");
	tree->Branch("time", event.time, "time[n]/D");
	tree->Branch("saturation", event.saturation, "saturation[n]/O");
	tree->Branch("camplitude", event.camplitude, "camp[n]/D");
	// free-tau fitting
	tree->Branch("tau", event.tau, "tau[n]/D");
	tree->Branch("base", event.base, "base[n]/D");
	tree->Branch("amplitude", event.amplitude, "amp[n]/D");
	tree->Branch("t0", event.t0, "t0[n]/D");
	// fix-tau fitting
	tree->Branch("ftau", event.ftau, "ftau[n]/D");
	tree->Branch("fbase", event.fbase, "fbase[n]/D");
	tree->Branch("famplitude", event.famplitude, "famp[n]/D");
	tree->Branch("ft0", event.ft0, "ft0[n]/D");
	// meta
	tree->Branch("run", &event.run, "run/I");
	tree->Branch("entry", &event.entry, "entry/I");
}


} // namespace name
