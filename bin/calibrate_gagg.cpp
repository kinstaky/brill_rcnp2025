#include <iostream>
#include <memory>

#include <TChain.h>
#include <TF1.h>
#include <TFile.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/gagg_event.h"
#include "include/utils.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/t0_utils.h"
#include "include/energy_calculator/delta_energy_calculator.h"
#include "include/t0/dssd.h"

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
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return 1;
	}
	char run_letter = run < 1079 && end_run < 1079 ? 'a'
		: run >= 1079 && end_run >= 1079 ? 'b'
		: '?';
	if (run_letter == '?') {
		std::cerr << "Error: Invalid run range " << run << "-" << end_run << ".\n";
		return 1;
	}

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

	TChain chain("tree");
	TChain gagg_chain("tree");
	int added_files = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsSkipRun(current_run)) continue;
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

	// load cuts
	std::vector<std::unique_ptr<TCutG>> be_cuts;
	be_cuts.resize(24);
	for (int i = 0; i < 24; ++i) {
		brill::ParseCutFile(
			config.root.workspace,
			"gagg_" + std::to_string(i) + run_letter,
			particle,
			false,
			be_cuts[i]
		);
	}

	// initialize delta energy calculator
	std::unique_ptr<brill::DeltaEnergyCalculator> particle_calculator = particle == "Be"
		? std::make_unique<brill::DeltaEnergyCalculator>(config, 4, 12, 500)
		: std::make_unique<brill::DeltaEnergyCalculator>(config, 2, 4, 500);

	// load T0 silicon calibration parameters
	brill::CalibrationParameters t0_cali_params(2);
	t0_cali_params.Read(cali_dir + "/t0.txt");

	TString output_path = TString::Format(
		"%s/gagg_%s_%04d_%04d.root",
		cali_dir.c_str(),
		particle.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TGraph gcali[25];

	// output parameters
	brill::t0::GAGGCalibrationParameters cali_param(25);

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
				if (!be_cuts[gagg_event.index[j]]->IsInside(
					gagg_event.amplitude[j],
					event.energy[i]
				)) continue;
				double x = gagg_event.amplitude[j];
				double de = t0_cali_params.p0[1] + t0_cali_params.p1[1] * event.energy[i];
				double y = particle_calculator->Energy(1, de);
				gcali[gagg_event.index[j]].AddPoint(x, y);
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	// fit
	if (particle == "Be") {
		constexpr double fit_range_a[] = {
			10000.0, 10000.0, 10000.0, 10000.0,
			10000.0, 10000.0, 10000.0, 10000.0,
			10000.0, 10000.0, 10000.0, 10000.0,
			10000.0, 10000.0, 10000.0, 10000.0,
			25000.0, 25000.0, 2500.0, 25000.0,
			25000.0, 25000.0, 25000.0, 25000.0,
			25000.0
		};
		constexpr double fit_range_b[] = {
			2000.0, 2000.0, 2000.0, 2000.0,
			2000.0, 2000.0, 2000.0, 2000.0,
			2000.0, 2000.0, 2000.0, 2000.0,
			2000.0, 2000.0, 2000.0, 2000.0,
			5000.0, 5000.0, 600.0, 5000.0,
			5000.0, 5000.0, 5000.0, 5000.0,
			5000.0
		};
		for (int i = 0; i < 25; ++i) {
			double fit_range = run_letter == 'a' ? fit_range_a[i] : fit_range_b[i];
			std::cout << "===========" << "Fitting GAGG " << i << "===========" << std::endl;
			TF1 flinear(TString::Format("fl%d", i), "pol1", 0.0, fit_range);
			flinear.SetLineColor(kBlue);
			gcali[i].Fit(&flinear, "R+");
			TF1 fcali(TString::Format("f%d", i), "[0]+[1]*x+[2]*exp(-x/[3])", 0.0, fit_range);
			fcali.SetParameter(0, flinear.GetParameter(0));
			fcali.SetParameter(1, flinear.GetParameter(1));
			fcali.SetParameter(2, -5.0);
			fcali.SetParameter(3, fit_range/20.0);
			if (run_letter == 'a' && i == 4) {
				fcali.SetParameter(2, -3.8);
				fcali.SetParameter(3, fit_range);
			}
			gcali[i].Fit(&fcali, "R+ ROB=0.7");
			cali_param.p0[i] = fcali.GetParameter(0);
			cali_param.p1[i] = fcali.GetParameter(1);
			cali_param.p2[i] = fcali.GetParameter(2);
			cali_param.p3[i] = fcali.GetParameter(3);
		}
	} else if (particle == "4He") {
		constexpr double fit_range_a[] = {
			16000.0, 16000.0, 16000.0, 16000.0,
			16000.0, 18000.0, 16000.0, 16000.0,
			16000.0, 16000.0, 16000.0, 16000.0,
			16000.0, 16000.0, 16000.0, 16000.0,
			54000.0, 54000.0, 3000.0, 54000.0,
			54000.0, 54000.0, 54000.0, 54000.0,
			54000.0
		};
		constexpr double fit_range_b[] = {
			5500.0, 5500.0, 5500.0, 5500.0,
			5500.0, 5500.0, 3000.0, 5500.0,
			5500.0, 5500.0, 5500.0, 5500.0,
			5500.0, 5500.0, 5500.0, 5500.0,
			16000.0, 16000.0, 1800.0, 18000.0,
			16000.0, 16000.0, 16000.0, 16000.0,
			16000.0
		};
		for (int i = 0; i < 25; ++i) {
			double fit_range = run_letter == 'a' ? fit_range_a[i] : fit_range_b[i];
			std::cout << "===========" << "Fitting GAGG " << i << "===========" << std::endl;
			TF1 flinear(TString::Format("fl%d", i), "pol1", 0.0, fit_range);
			flinear.SetLineColor(kBlue);
			gcali[i].Fit(&flinear, "R+");
			TF1 fcali(TString::Format("f%d", i), "[0]+[1]*x+[2]*exp(-x/[3])", 0.0, fit_range);
			fcali.SetParameter(0, flinear.GetParameter(0));
			fcali.SetParameter(1, flinear.GetParameter(1));
			fcali.SetParameter(2, -5.0);
			fcali.SetParameter(3, fit_range/20.0);
			gcali[i].Fit(&fcali, "R+ ROB=0.9");
			cali_param.p0[i] = fcali.GetParameter(0);
			cali_param.p1[i] = fcali.GetParameter(1);
			cali_param.p2[i] = fcali.GetParameter(2);
			cali_param.p3[i] = fcali.GetParameter(3);
		}
	}

	// save graph
	opf.cd();
	for (int i = 0; i < 25; ++i) {
		gcali[i].Write(TString::Format("g%d", i));
	}
	// close files
	opf.Close();

	// save parameters
	TString cali_output = TString::Format(
		"%s/gagg_layer1_%c_%s.txt",
		cali_dir.c_str(),
		run_letter,
		particle.c_str()
	);
	cali_param.Write(cali_output.Data());
	return 0;
}
