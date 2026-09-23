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
#include "include/event/gagg_event.h"
#include "include/utils.h"
#include "include/t0_utils.h"
#include "include/t0/dssd.h"

constexpr int kMaxHit = 8;

struct ParticleCutInfo {
	std::string particle;
	int charge = 0;
	int mass = 0;
	std::unique_ptr<TCutG> cut = nullptr;
};

struct GAGGStraightCut {
	std::string particle;
	brill::t0::GAGGStraightParameters param;
	double min_ef;
	double max_ef;
	double max_gagg_energy;
};

const char *const kD1D2StopParticles[] = {"1H", "2H", "3H", "4He", "6He", "7Li", "10Be", "12Be"};
const char *const kD1D2PassParticles[] = {"4He", "12Be"};

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

int main(int argc, char **argv) {
	cxxopts::Options options("track_t0", "Track and identify T0 particles.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
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
		return -1;
	}

	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		return -1;
	}
	const int run = result["run"].as<int>();
	if (config.IsJumpRun(run)) {
		std::cout << "Skipping jump run " << run << ".\n";
		return 0;
	}

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.root.workspace, config.paths.ingot);
	const std::string cali_dir = brill::JoinPath(config.root.workspace, config.paths.calibration);

	// load PID cuts
	std::vector<ParticleCutInfo> d1d2_stop_cuts;
	std::vector<ParticleCutInfo> d1d2_tail_cuts;
	for (const auto &particle : kD1D2StopParticles) {
		ParticleCutInfo cut_info;
		cut_info.particle = particle;
		if (!ParseParticleName(cut_info.particle, cut_info.charge, cut_info.mass)) continue;
		if (brill::ParseCutFile(
			config.root.workspace,
			"d1d2",
			particle,
			false,
			cut_info.cut
		)) {
			std::cerr << "Error: Load stopped cut for " << particle << " failed.\n";
			return -2;
		}
		d1d2_stop_cuts.push_back(std::move(cut_info));
	}
	for (const auto &particle : kD1D2PassParticles) {
		ParticleCutInfo cut_info;
		cut_info.particle = particle;
		if (!ParseParticleName(cut_info.particle, cut_info.charge, cut_info.mass)) continue;
		if (brill::ParseCutFile(
			config.root.workspace,
			"d1d2",
			particle,
			true,
			cut_info.cut
		)) {
			std::cerr << "Error: Load tail cut for " << particle << " failed.\n";
			return -2;
		}
		d1d2_tail_cuts.push_back(std::move(cut_info));
	}
	std::map<std::string, GAGGStraightCut> gagg_cuts;
	gagg_cuts.insert(std::make_pair<std::string, GAGGStraightCut>(
		"12Be", {"12Be", {}, 121.08, 126.86, 90.0}
	));
	gagg_cuts["12Be"].param.Read(cali_dir + "/gagg_straight_Be.txt");
	gagg_cuts.insert(std::make_pair<std::string, GAGGStraightCut>(
		"4He", {"4He", {}, 35.2, 37.7, 80.0}
	));
	gagg_cuts["4He"].param.Read(cali_dir + "/gagg_straight_4He.txt");


	// load calibration parameters
	brill::CalibrationParameters t0_cali(2);
	if (t0_cali.Read(cali_dir + "/t0.txt")) {
		std::cerr << "Error: Failed to read t0 calibration parameters.\n";
		return -1;
	}
	brill::t0::GAGGCalibrationParameters gagg_cali(25);
	TString gagg_cali_path = TString::Format(
		"%s/gagg_layer1_%c_Be.txt",
		cali_dir.c_str(),
		run < 1079 ? 'a' : 'b'
	);
	if (gagg_cali.Read(gagg_cali_path.Data())) {
		std::cerr << "Error: Failed to read gagg calibration parameters.\n";
		return -1;
	}

	TChain chain1("tree");
	chain1.Add(TString::Format(
		"%s/t0d1_%04d.root",
		match_dir.c_str(),
		run
	));
	TChain chain2("tree");
	chain2.Add(TString::Format(
		"%s/t0d2_%04d.root",
		match_dir.c_str(),
		run
	));
	TChain chain_gagg("tree");
	chain_gagg.Add(TString::Format(
		"%s/gagg_%04d.root",
		ingot_dir.c_str(),
		run
	));
	chain1.AddFriend(&chain2, "d2");
	chain1.AddFriend(&chain_gagg, "gagg");

	brill::DssdMatchEvent d1_event;
	brill::DssdMatchEvent d2_event;
	brill::GaggEvent gagg_event;
	brill::SetupInput(&chain1, d1_event);
	brill::SetupInput(&chain1, d2_event, "d2.");
	brill::SetupInput(&chain1, gagg_event, "gagg.");

	TString output_path = TString::Format(
		"%s/t0_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.telescope).c_str(),
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
		int used_d1 = 0;
		int used_d2 = 0;
		int used_gagg = 0;

		for (int i = 0; i < d1_event.num; ++i) {
			if (used_d1 & (1 << i)) continue;
			double d1_energy = t0_cali.p0[0] + t0_cali.p1[0] * d1_event.energy[i];
			for (int j = 0; j < d2_event.num; ++j) {
				if (used_d1 & (1 << i)) break;
				if (used_d2 & (1 << j)) continue;
				double d2_energy = t0_cali.p0[1] + t0_cali.p1[1] * d2_event.energy[j];
				// check tracking window
				if (!brill::t0::IsInTrackWindow(
					d1_event.front_strip[i],
					d1_event.back_strip[i],
					d2_event.front_strip[j],
					d2_event.back_strip[j]
				)) continue;
				for (int k = 0; k < gagg_event.num; ++k) {
					if (used_d2 & (1 << j)) break;
					if (used_gagg & (1 << k)) continue;
					if (gagg_event.index[k] >= 24) continue;
					double gagg_energy = gagg_cali.CaliEnergy(gagg_event.index[k], gagg_event.energy[k]);
					if (!brill::t0::IsInTrackWindow(
						d2_event.front_strip[j],
						d2_event.back_strip[j],
						gagg_event.index[k]
					)) continue;
					for (const auto &cut : d1d2_tail_cuts) {
						if (!cut.cut->IsInside(d2_event.energy[j], d1_event.energy[i])) continue;
						GAGGStraightCut &gagg_cut = gagg_cuts[cut.particle];
						double ef = gagg_cut.param.FixedEnergy(d2_energy, gagg_energy);
						if (ef < gagg_cut.min_ef || ef > gagg_cut.max_ef || gagg_energy > gagg_cut.max_gagg_energy) continue;
						t0_event.layer[num] = 3;
						t0_event.flag[num] = 0x7;
						t0_event.charge[num] = cut.charge;
						t0_event.mass[num] = cut.mass;
						t0_event.energy[num][0] = d1_energy;
						t0_event.energy[num][1] = d2_energy;
						t0_event.energy[num][2] = gagg_energy;
						t0_event.time[num][0] = d1_event.time[i];
						t0_event.time[num][1] = d2_event.time[j];
						t0_event.time[num][2] = gagg_event.time[k];
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
						t0_event.scintillator_index[num] = gagg_event.index[k];
						t0_event.last[t0_event.num][0] = i;
						t0_event.last[t0_event.num][1] = j;
						t0_event.last[t0_event.num][2] = k;
						t0_event.num++;
						used_d1 |= (1 << i);
						used_d2 |= (1 << j);
						used_gagg |= (1 << k);
						break;
					}
				}
				if (used_d2 & (1 << j)) continue;
				// check PID cuts
				for (const auto &cut : d1d2_stop_cuts) {
					if (!cut.cut->IsInside(d2_event.energy[j], d1_event.energy[i])) continue;
					t0_event.layer[num] = 2;
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
					used_d1 |= (1 << i);
					used_d2 |= (1 << j);
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