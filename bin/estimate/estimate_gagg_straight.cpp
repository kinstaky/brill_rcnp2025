#include <cmath>
#include <iostream>
#include <memory>

#include <TCutG.h>
#include <TChain.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TF1.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/utils.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/event/gagg_event.h"
#include "include/t0_utils.h"
#include "include/t0/dssd.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double PidFit(double *x, double *par) {
	return (0.5 / par[0]) * (
		std::sqrt(std::pow(x[0], 2.0) + 4.0 * par[0] * std::pow(par[1] * x[0] - par[2], 2.0))
		- x[0]
	);
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_straight", "Estimate T0 straight PID.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
		("p,particle", "Particle type to use, Be, He.", cxxopts::value<std::string>(), "particle")
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
	std::string particle = "Be";
	if (result.count("particle")) {
		if (result["particle"].as<std::string>() == "He") {
			particle = "4He";
		} else if (result["particle"].as<std::string>() != "Be") {
			std::cerr << "Error: Invalid particle " << result["particle"].as<std::string>() << ".\n";
			return -1;
		}
	}

	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		return -1;
	}

	const int run = result["run"].as<int>();
	const int end_run = result.count("end-run") ? result["end-run"].as<int>() : run;
	if (end_run < run) {
		std::cerr << "Error: end run " << end_run << " is smaller than run " << run << ".\n";
		return -1;
	}

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.root.workspace, config.paths.ingot);
	const std::string cali_dir = brill::JoinPath(config.root.workspace, config.paths.calibration);

	TChain chain_d2("tree");
	TChain chain_gagg("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsSkipRun(current_run)) continue;
		++added_runs;
		chain_d2.Add(TString::Format(
			"%s/t0d2_%04d.root",
			match_dir.c_str(),
			current_run
		));
		chain_gagg.Add(TString::Format(
			"%s/gagg_%04d.root",
			ingot_dir.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}

	chain_d2.AddFriend(&chain_gagg, "gagg");
	brill::DssdMatchEvent d2_event;
	brill::GaggEvent gagg_event;
	brill::SetupInput(&chain_d2, d2_event);
	brill::SetupInput(&chain_d2, gagg_event, "gagg.");

	// Load cuts
	std::unique_ptr<TCutG> cuts_a[25], cuts_b[25];
	for (int i = 0; i < 24; ++i) {
		if (brill::ParseCutFile(
			config.root.workspace,
			"gagg_" + std::to_string(i) + "a",
			particle,
			false,
			cuts_a[i]
		)) {
			std::cerr << "Error: Parse GAGG " << i << " 12Be cut a failed.\n";
			return -2;
		}
		if (brill::ParseCutFile(
			config.root.workspace,
			"gagg_" + std::to_string(i) + "b",
			particle,
			false,
			cuts_b[i]
		)) {
			std::cerr << "Error: Parse GAGG " << i << " 12Be cut b failed.\n";
			return -2;
		}
	}

	// Load calibration parameters
	brill::CalibrationParameters t0_cali(2);
	t0_cali.Read(cali_dir + "/t0.txt");
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

	TString output_path = TString::Format(
		"%s/gagg_straight_%s_%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.estimate).c_str(),
		particle.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TGraph gcurve;
	std::vector<TH2F> pids, spids;
	std::vector<TH1F> epids;
	double gagg_max = particle == "Be" ? 300.0 : 100.0;
	for (int i = 0; i < 25; ++i) {
		pids.emplace_back(
			TString::Format("p%d", i),
			TString::Format("PID of GAGG %d", i),
			1000, 0, gagg_max, 1000, 0, 300.0
		);
		spids.emplace_back(
			TString::Format("s%d", i),
			TString::Format("straight PID of GAGG %d", i),
			1000, 0, gagg_max, 1000, 0, 300.0
		);
		epids.emplace_back(
			TString::Format("e%d", i),
			TString::Format("Straight energy of GAGG %d", i),
			1000, 0, 300.0
		);
	}
	TH2F pid("pid", "PID of all GAGG", 1000, 0, gagg_max, 1000, 0, 300.0);
	TH2F spid("spid", "straight PID of all GAGG", 1000, 0, gagg_max, 1000, 0, 300.0);
	TH1F epid("epid", "Straight energy of all GAGG", 1000, 0, 300.0);

	const long long total = chain_d2.GetEntries();
	long long last_percentage = -1;
	std::printf("Filling GAGG straight fit points   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain_d2.GetEntry(entry);
		for (int i = 0; i < d2_event.num; ++i) {
			for (int j = 0; j < gagg_event.num; ++j) {
				if (gagg_event.index[j] >= 25) continue;
				// position check
				if (!brill::t0::IsInTrackWindow(
					d2_event.front_strip[i],
					d2_event.back_strip[i],
					gagg_event.index[j]
				)) continue;
				// PID check
				std::unique_ptr<TCutG> &cut = d2_event.run < 1079
					? cuts_a[gagg_event.index[j]]
					: cuts_b[gagg_event.index[j]];
				if (!cut->IsInside(gagg_event.amplitude[j], d2_event.energy[i])) continue;
				// calibrated energy
				double d2_energy = t0_cali.p0[1] + t0_cali.p1[1] * d2_event.energy[i];
				double gagg_energy = d2_event.run < 1079
					? cali_param_a.CaliEnergy(gagg_event.index[j], gagg_event.amplitude[j])
					: cali_param_b.CaliEnergy(gagg_event.index[j], gagg_event.amplitude[j]);
				gcurve.AddPoint(gagg_energy, d2_energy);
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	// fit parameters
	double a = 0.29;
	double b = -0.266;
	double c = 100.0;

	TF1 *f1 = new TF1("f1", PidFit, 0.0, gagg_max, 3);
	f1->SetParameter(0, a);
	f1->SetParameter(1, b);
	f1->SetParameter(2, c);
	gcurve.Fit(f1, "RQ+ ROB=0.8");
	a = f1->GetParameter(0);
	b = f1->GetParameter(1);
	c = f1->GetParameter(2);
	std::cout
		<< particle
		<< ": A " << a
		<< ", B " << b
		<< ", C " << c
		<< "\n";

	last_percentage = -1;
	std::printf("Filling GAGG straight PID   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain_d2.GetEntry(entry);
		for (int i = 0; i < d2_event.num; ++i) {
			for (int j = 0; j < gagg_event.num; ++j) {
				if (gagg_event.index[j] >= 25) continue;
				// position check
				if (!brill::t0::IsInTrackWindow(
					d2_event.front_strip[i],
					d2_event.back_strip[i],
					gagg_event.index[j]
				)) continue;
				// calibrated energy
				double d2_energy = t0_cali.p0[1] + t0_cali.p1[1] * d2_event.energy[i];
				double gagg_energy = d2_event.run < 1079
					? cali_param_a.CaliEnergy(gagg_event.index[j], gagg_event.amplitude[j])
					: cali_param_b.CaliEnergy(gagg_event.index[j], gagg_event.amplitude[j]);
				double ef = std::sqrt(d2_energy*gagg_energy + a*d2_energy*d2_energy) + b*gagg_energy;
				pid.Fill(gagg_energy, d2_energy);
				spid.Fill(gagg_energy, ef);
				epid.Fill(ef);
				pids[gagg_event.index[j]].Fill(gagg_energy, d2_energy);
				spids[gagg_event.index[j]].Fill(gagg_energy, ef);
				epids[gagg_event.index[j]].Fill(ef);
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	gcurve.Write("gfit");
	for (auto &hist : pids) hist.Write();
	for (auto &hist : spids) hist.Write();
	for (auto &hist : epids) hist.Write();
	pid.Write();
	spid.Write();
	epid.Write();
	opf.Close();

	brill::t0::GAGGStraightParameters straight_param;
	straight_param.a = a;
	straight_param.b = b;
	straight_param.Write(TString::Format(
		"%s/gagg_straight_%s.txt",
		cali_dir.c_str(),
		particle.c_str()
	).Data());

	return 0;
}