#include <iostream>
#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/t0/dssd.h"
#include "include/event/gagg_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/t0_utils.h"
#include "include/utils.h"

struct Event {
	int run;
	int entry;
	int index;
	double d2_energy;
	double gagg_energy;
};

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("calibrate_t0", "Calibrate T0 energies from tracked events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("p,particle", "Be or He", cxxopts::value<std::string>(), "particle")
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

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;

	std::string particle = "Be";
	if (result.count("particle")) {
		if (result["particle"].as<std::string>() == "He") {
			particle = "4He";
		} else if (result["particle"].as<std::string>() != "Be") {
			std::cerr << "Error: Invalid particle " << result["particle"].as<std::string>() << ".\n";
			return 2;
		}
	}

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.root.workspace, config.paths.ingot);
	const std::string cali_dir = brill::JoinPath(config.root.workspace, config.paths.calibration);
	const std::string estimate_dir = brill::JoinPath(config.root.workspace, config.paths.estimate);

	TChain chain("tree");
	TChain gagg_chain("tree");
	int added_files = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
		std::string path = TString::Format(
			"%s/t0d2_%04d.root",
			match_dir.c_str(),
			current_run
		).Data();
		if (!std::filesystem::exists(path)) {
			std::cerr << "Warning: Skip missing tracked file " << path << ".\n";
			continue;
		}
		chain.Add(path.c_str());

		std::string gagg_path = TString::Format(
			"%s/gagg_%04d.root",
			ingot_dir.c_str(),
			current_run
		).Data();
		if (!std::filesystem::exists(gagg_path)) {
			std::cerr << "Warning: Skip missing ingot file " << gagg_path << ".\n";
			continue;
		}
		gagg_chain.Add(gagg_path.c_str());
		++added_files;
	}
	if (added_files == 0) {
		std::cout << "No tracked files to process.\n";
		return 0;
	}
	chain.AddFriend(&gagg_chain, "gagg");

	brill::DssdMatchEvent event;
	brill::SetupInput(&chain, event);
	brill::GaggEvent gagg_event;
	brill::SetupInput(&chain, gagg_event, "gagg.");

	// load T0 silicon calibration parameters
	brill::CalibrationParameters t0_cali_params(2);
	t0_cali_params.Read(cali_dir + "/t0.txt");

	TString output_path = TString::Format(
		"%s/gagg_calibration_%s_%04d_%04d.root",
		estimate_dir.c_str(),
		particle.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	double gagg_max = particle == "Be" ? 300.0 : 100.0;
	TH2F pid("pid", "D2 energy vs. GAGG energy", 1000, 0, gagg_max, 1000, 0, 80000);
	TTree opt("tree", "Calibrated GAGG events");
	Event out_event;
	opt.Branch("run", &out_event.run, "run/I");
	opt.Branch("entry", &out_event.entry, "entry/I");
	opt.Branch("index", &out_event.index, "index/I");
	opt.Branch("d2_energy", &out_event.d2_energy, "de/D");
	opt.Branch("gagg_energy", &out_event.gagg_energy, "e/D");

	// // load T0 silicon calibration parameters
	// brill::CalibrationParameters t0_cali_params(2);
	// t0_cali_params.Read(cali_dir + "/t0.txt");

	// calibration parameters
	brill::t0::GAGGCalibrationParameters cali_param_a(25);
	TString cali_param_path_a = TString::Format(
		"%s/gagg_layer1_a_%s.txt",
		cali_dir.c_str(),
		particle.c_str()
	);
	cali_param_a.Read(cali_param_path_a.Data());
	brill::t0::GAGGCalibrationParameters cali_param_b(25);
	TString cali_param_path_b = TString::Format(
		"%s/gagg_layer1_b_%s.txt",
		cali_dir.c_str(),
		particle.c_str()
	);
	cali_param_b.Read(cali_param_path_b.Data());

	long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Collecting GAGG calibration samples   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain.GetEntry(entry);
		for (int i = 0; i < event.num; ++i) {
			for (int j = 0; j < gagg_event.num; ++j) {
				if (gagg_event.index[j] >= 24) continue;
				if (!brill::t0::IsInTrackWindow(
					event.front_strip[i],
					event.back_strip[i],
					gagg_event.index[j])
				) continue;
				out_event.run = event.run;
				out_event.entry = entry;
				out_event.index = gagg_event.index[j];
				out_event.d2_energy = event.energy[i];
				// out_event.d2_energy = t0_cali_params.p0[1] + t0_cali_params.p1[1] * event.energy[i];
				out_event.gagg_energy = event.run < 1079
					? cali_param_a.CaliEnergy(out_event.index, gagg_event.amplitude[j])
					: cali_param_b.CaliEnergy(out_event.index, gagg_event.amplitude[j]);
				pid.Fill(out_event.gagg_energy, out_event.d2_energy);
				opt.Fill();
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	// save graph
	opf.cd();
	pid.Write();
	opt.Write();
	// close files
	opf.Close();

	return 0;
}

