#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include <TChain.h>
#include <TString.h>
#include <TFile.h>
#include <TF1.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/energy_calculator/delta_energy_calculator.h"
#include "include/event/t0/t0_event.h"
#include "include/utils.h"


constexpr int kLayerCount = 2;

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

struct ParticlePidInfo {
	// charge number of this particle
	int charge;
	// mass number of this particle
	int mass;
	// left bound of pid
	double left;
	// right bound of pid
	double right;
};

const std::vector<ParticlePidInfo> pid_info {
	// d1d2
	{1, 2, 1400.0, 2800.0},
	{2, 4, 4300.0, 8800.0},
	{3, 7, 8800.0, 16200.0},
	{4, 12, 16200.0, 38000.0},
};

class PidFitFunc {
public:
	PidFitFunc(
		const std::vector<std::pair<int, int>> &projectiles,
		const brill::AppConfig &config
	) {
		for (const auto &p : projectiles) {
			calculators_.insert(
				std::make_pair(
					p.first * 100 + p.second,
					std::make_unique<brill::DeltaEnergyCalculator>(
						config,
						p.first,
						p.second
					)
				)
			);
		}
	}

	double operator()(double *x, double *par) const {
		// identify particle
		for (const auto &info : pid_info) {
			if (x[0] > info.left && x[0] < info.right) {
				double de = par[0] + par[1] * x[0];
				auto calculator = calculators_.at(info.charge*100+info.mass);
				double e = calculator->Energy(0, de);
				return (e - par[2]) / par[3];
			}
		}
		return 0.0;
	}
private:
	std::map<int, std::shared_ptr<brill::DeltaEnergyCalculator>> calculators_;
};

int WriteCalibrationParameters(
	const char *path,
	const double *parameters
) {
	std::filesystem::path file_path(path);
	if (!file_path.parent_path().empty()) {
		std::filesystem::create_directories(file_path.parent_path());
	}
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output parameter file " << path << " failed.\n";
		return -1;
	}

	fout << "# layer p0 p1\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		fout << layer << " " << parameters[layer * 2] << " " << parameters[layer * 2 + 1] << "\n";
	}
	return 0;
}

int main(int argc, char **argv) {
	cxxopts::Options options("calibrate_t0", "Calibrate T0 energies from tracked events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
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
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}

	const brill::SiliconDetectorConfig *detectors[kLayerCount] = {
		config.FindDetector("t0d1"),
		config.FindDetector("t0d2")
	};
	for (int i = 0; i < kLayerCount; ++i) {
		if (!detectors[i]) {
			std::cerr << "Error: Missing detector config for T0D" << i+1 << ".\n";
			return 1;
		}
	}

	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);
	const std::string tele_dir = brill::JoinPath(config.root.workspace, config.paths.telescope);

	TChain chain("tree");
	int added_files = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
		std::string path = TString::Format(
			"%s/t0_%s%04d.root",
			tele_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		).Data();
		if (!std::filesystem::exists(path)) {
			std::cerr << "Warning: Skip missing tracked file " << path << ".\n";
			continue;
		}
		chain.Add(path.c_str());
		++added_files;
	}
	if (added_files == 0) {
		std::cout << "No tracked files to process.\n";
		return 0;
	}

	brill::T0Event event;
	brill::SetupInput(&chain, event);

	TString output_path = TString::Format(
		"%s/t0_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.calibration).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TGraph gcali;

	long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Collecting T0 calibration samples   0%%");
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
			for (const auto &info : pid_info) {
				if (
					event.charge[i] == info.charge
					&& event.mass[i] == info.mass
					&& event.energy[i][0] > info.left
					&& event.energy[i][0] < info.right
				) {
					gcali.AddPoint(
						event.energy[i][0],
						event.energy[i][1]
					);
					break;
				}
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	double initial_calibration_parameters[4] = {
		-0.3, 0.00298,
		4.2, 0.00281,
	};

	std::vector<std::pair<int, int>> projectiles {
		// 2H
		{1, 2},
		// 4He
		{2, 4},
		// 7Li
		{3, 7},
		// 12Be
		{4, 12}
	};
	PidFitFunc pid_fit(projectiles, config);
	// fit function
	TF1 fcali("fcali", pid_fit, 0.0, 40'000.0, 4);
	fcali.SetNpx(1000);
	for (size_t i = 0; i < 4; ++i) {
		fcali.SetParameter(i, initial_calibration_parameters[i]);
	}
	// fcali.SetParLimits(0, -5.0, 5.0);
	// fcali.SetParLimits(1, 0.002, 0.004);
	// fcali.SetParLimits(2, -5.0, 5.0);
	// fcali.SetParLimits(3, 0.002, 0.004);
	// fit
	gcali.Fit(&fcali, "R+");

	double *parameters = fcali.GetParameters();
	std::cout << "Calibration parameters:\n";
	for (int layer = 0; layer < kLayerCount; ++layer) {
		std::cout
			<< "  T0D" << layer + 1
			<< ": p0 = " << parameters[layer * 2]
			<< ", p1 = " << parameters[layer * 2 + 1]
			<< "\n";
	}

	output_path = TString::Format(
		"%s/t0.txt",
		brill::JoinPath(config.root.workspace, config.paths.calibration).c_str()
	).Data();
	WriteCalibrationParameters(output_path, parameters);
	// save graph
	opf.cd();
	gcali.Write("gcali");
	// close files
	opf.Close();
	return 0;
}
