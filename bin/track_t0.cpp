#include <cstdio>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <TChain.h>
#include <TCutG.h>
#include <TFile.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"
#include "include/t0_utils.h"

constexpr int kMaxHit = 8;

struct ParticleCutInfo {
	std::string particle;
	int charge = 0;
	int mass = 0;
	std::unique_ptr<TCutG> cut;
};

const char *const kD1D2Particles[] = {"1H", "2H", "3H", "4He", "6He", "7Li", "10Be", "12Be"};

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int ElementCharge(const std::string &element) {
	if (element == "H") return 1;
	if (element == "He") return 2;
	if (element == "Li") return 3;
	if (element == "Be") return 4;
	if (element == "B") return 5;
	if (element == "C") return 6;
	if (element == "N") return 7;
	if (element == "O") return 8;
	return 0;
}

std::string ParticleElement(const std::string &particle) {
	size_t index = 0;
	while (index < particle.size() && std::isdigit(static_cast<unsigned char>(particle[index]))) {
		++index;
	}
	return particle.substr(index);
}

bool ParseParticleName(const std::string &particle, int &charge, int &mass) {
	size_t index = 0;
	while (index < particle.size() && std::isdigit(static_cast<unsigned char>(particle[index]))) {
		++index;
	}
	if (index == 0 || index >= particle.size()) return false;
	mass = std::stoi(particle.substr(0, index));
	charge = ElementCharge(particle.substr(index));
	return charge > 0;
}

bool CutFileExists(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	bool tail
) {
	std::ifstream fin(brill::CutFilePath(workspace, slice, particle, tail));
	return fin.good();
}

int LoadTailCut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	std::set<std::string> &loaded,
	std::unique_ptr<TCutG> &cut
) {
	if (CutFileExists(workspace, slice, particle, true)) {
		return brill::ParseCutFile(workspace, slice, particle, true, cut);
	}
	std::string element = ParticleElement(particle);
	if (loaded.count(element) == 1) return 0;
	if (
		!element.empty()
		&& element != particle
		&& CutFileExists(workspace, slice, element, true)
	) {
		loaded.insert(element);
		return brill::ParseCutFile(workspace, slice, element, true, cut);
	}
	std::cerr
		<< "Error: Missing tail cut file for " << slice
		<< " particle " << particle << ".\n";
	return -1;
}

int LoadStoppedCut(
	const std::string &workspace,
	const std::string &slice,
	const std::string &particle,
	std::unique_ptr<TCutG> &cut
) {
	if (!CutFileExists(workspace, slice, particle, false)) {
		std::cerr
			<< "Error: Missing stopped cut file for " << slice
			<< " particle " << particle << ".\n";
		return -1;
	}
	return brill::ParseCutFile(workspace, slice, particle, false, cut);
}

void BuildTailCuts(
	const std::string &workspace,
	const std::string &slice,
	const char *const *particles,
	size_t particle_count,
	std::vector<ParticleCutInfo> &cuts
) {
	cuts.clear();
	std::set<std::string> loaded;
	std::set<std::string> loaded_tail;
	for (size_t i = 0; i < particle_count; ++i) {
		std::string particle = particles[i];
		if (loaded.count(particle) != 0) continue;
		loaded.insert(particle);
		ParticleCutInfo cut_info;
		cut_info.particle = particle;
		if (!ParseParticleName(cut_info.particle, cut_info.charge, cut_info.mass)) continue;
		if (LoadTailCut(workspace, slice, cut_info.particle, loaded_tail, cut_info.cut)) {
			throw std::runtime_error("load tail cut failed");
		}
		cuts.push_back(std::move(cut_info));
	}
}

void BuildStoppedCuts(
	const std::string &workspace,
	const std::string &slice,
	const char *const *particles,
	size_t particle_count,
	std::vector<ParticleCutInfo> &cuts
) {
	cuts.clear();
	for (size_t i = 0; i < particle_count; ++i) {
		ParticleCutInfo cut_info;
		cut_info.particle = particles[i];
		if (!ParseParticleName(cut_info.particle, cut_info.charge, cut_info.mass)) continue;
		if (LoadStoppedCut(workspace, slice, cut_info.particle, cut_info.cut)) {
			throw std::runtime_error("load stopped cut failed");
		}
		cuts.push_back(std::move(cut_info));
	}
}

int main(int argc, char **argv) {
	cxxopts::Options options("track_t0", "Track and identify T0 particles.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		("t,trigger", "Trigger type.", cxxopts::value<std::string>(), "trigger")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		);

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		PrintUsage(options);
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		PrintUsage(options);
		return 1;
	}

	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		return 1;
	}
	if (result.count("trigger")) {
		config.root.trigger = result["trigger"].as<std::string>();
	}
	const int run = result["run"].as<int>();
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	std::vector<ParticleCutInfo> d1d2_stop_cuts;
	// std::vector<ParticleCutInfo> d1d2_tail_cuts;
	try {
		BuildStoppedCuts(
			config.root.workspace, "d1d2", kD1D2Particles,
			sizeof(kD1D2Particles) / sizeof(kD1D2Particles[0]), d1d2_stop_cuts
		);
		// BuildTailCuts(
		// 	config.root.workspace, "t0d1d2", kD1D2Particles,
		// 	sizeof(kD1D2Particles) / sizeof(kD1D2Particles[0]), d1d2_tail_cuts
		// );
	} catch (const std::runtime_error &) {
		return -1;
	}

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);

	TChain chain1("tree");
	chain1.Add(TString::Format(
		"%s/t0d1_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	TChain chain2("tree");
	chain2.Add(TString::Format(
		"%s/t0d2_%s%04d.root",
		match_dir.c_str(),
		trigger_infix.c_str(),
		run
	));
	chain1.AddFriend(&chain2, "d2");

	brill::DssdMatchEvent d1_event;
	brill::DssdMatchEvent d2_event;
	brill::SetupInput(&chain1, d1_event);
	brill::SetupInput(&chain1, d2_event, "d2.");

	TString output_path = TString::Format(
		"%s/t0_%s%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.telescope).c_str(),
		trigger_infix.c_str(),
		run
	);
	TFile opf(output_path, "recreate");
	TTree opt("tree", "tracked t0");
	brill::T0Event t0_event;
	brill::Reset(t0_event);
	brill::SetupOutput(&opt, t0_event);

	const brill::SiliconDetectorConfig *t0d1_config = config.FindDetector("t0d1");
	const brill::SiliconDetectorConfig *t0d2_config = config.FindDetector("t0d2");

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Tracking T0   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		brill::Reset(t0_event);
		t0_event.run = run;
		t0_event.entry = entry;
		int &num = t0_event.num;

		for (int i = 0; i < d1_event.num; ++i) {
			for (int j = 0; j < d2_event.num; ++j) {
				// check tracking window
				if (!brill::t0::IsInTrackWindow(
					d1_event.front_strip[i],
					d1_event.back_strip[i],
					d2_event.front_strip[j],
					d2_event.back_strip[j]
				)) continue;
				// check PID cuts
				for (const auto &cut : d1d2_stop_cuts) {
					if (!cut.cut->IsInside(d2_event.energy[j], d1_event.energy[i])) continue;
					t0_event.layer[num] = 1;
					t0_event.flag[num] = 0x3;
					t0_event.charge[num] = cut.charge;
					t0_event.mass[num] = cut.mass;
					t0_event.energy[num][0] = d1_event.energy[i];
					t0_event.energy[num][1] = d2_event.energy[j];
					t0_event.time[num][0] = d1_event.time[i];
					t0_event.time[num][1] = d2_event.time[j];
					brill::t0::GetPixelPosition(
						t0d1_config,
						d1_event.front_strip[i],
						d1_event.back_strip[j],
						t0_event.x[num][0],
						t0_event.y[num][0],
						t0_event.z[num][0]
					);
					brill::t0::GetPixelPosition(
						t0d2_config,
						d2_event.front_strip[i],
						d2_event.back_strip[j],
						t0_event.x[num][1],
						t0_event.y[num][1],
						t0_event.z[num][1]
					);
					t0_event.last[t0_event.num][0] = i;
					t0_event.last[t0_event.num][1] = j;
					t0_event.num++;
					break;
				}
			}
		}

		opt.Fill();
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	opt.Write();
	opf.Close();
	return 0;
}