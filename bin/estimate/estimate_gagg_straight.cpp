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

struct StraightSlice {
	std::string key;
	TCutG *cut;
	TGraph graph;
	std::unique_ptr<TF1> fit;
	double a = 0.5;
	double b = -0.04;
	double c = 0.0;

	StraightSlice(const std::string &k, TCutG* c) : key(k), cut(c) {}
};


void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

double PidFit(double *x, double *par) {
	return (0.5 / par[0]) * (
		std::sqrt(std::pow(x[0], 2.0) + 4.0 * par[0] * std::pow(par[1] * x[0] - par[2], 2.0))
		- x[0]
	);
}

int FitSlice(StraightSlice &slice) {
	if (slice.graph.GetN() < 5 || !slice.cut) {
		return -1;
	}

	double x_min = slice.cut->GetX()[0];
	double x_max = slice.cut->GetX()[0];
	double y_max = slice.cut->GetY()[0];
	for (int i = 1; i < slice.cut->GetN(); ++i) {
		if (slice.cut->GetX()[i] < x_min) x_min = slice.cut->GetX()[i];
		if (slice.cut->GetX()[i] > x_max) x_max = slice.cut->GetX()[i];
		if (slice.cut->GetY()[i] > y_max) y_max = slice.cut->GetY()[i];
	}

	std::string fit_name = "f" + slice.key;
	slice.fit = std::make_unique<TF1>(fit_name.c_str(), PidFit, x_min, x_max, 3);
	slice.fit->SetParameter(0, slice.a);
	slice.fit->SetParameter(1, slice.b);
	slice.fit->SetParameter(2, y_max);
	slice.fit->SetParLimits(2, 0.0, 1e10);
	slice.graph.Fit(slice.fit.get(), "RQ+ ROB=0.8");
	slice.a = slice.fit->GetParameter(0);
	slice.b = slice.fit->GetParameter(1);
	slice.c = slice.fit->GetParameter(2);
	std::cout
		<< "  " << slice.key
		<< ": A " << slice.a
		<< ", B " << slice.b
		<< ", C " << slice.c
		<< "\n";
	return 0;
}

void FillGraphs(
	const double a,
	const double b,
	const double de,
	const double e,
	TH2F &pid,
	TH2F &spid,
	TH1F &epid
) {
	double ef = std::sqrt(de * e + a * de * de) + b * e;
	pid.Fill(e, de);
	spid.Fill(e, ef);
	epid.Fill(ef);
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
		return 1;
	}
	int use_particle = 0;
	if (result.count("particle")) {
		if (result["particle"].as<std::string>() == "Be") {
			use_particle = 0;
		} else if (result["particle"].as<std::string>() == "He") {
			use_particle = 1;
		}
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

	std::string run_letter = run < 1079 && end_run < 1079 ? "a"
		: run >= 1079 && end_run >= 1079 ? "b"
		: "";
	if (run_letter.empty()) {
		std::cerr << "Error: Invalid run range " << run << "-" << end_run << ".\n";
		return 1;
	}

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.root.workspace, config.paths.ingot);

	TChain chain_d2("tree");
	TChain chain_gagg("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
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

	std::vector<StraightSlice> slices[25];

	std::unique_ptr<TCutG> be_cut_0_15, he_cut_0_15;
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_0_15_"+run_letter,
		"Be",
		false,
		be_cut_0_15
	)) {
		std::cerr << "Error: Parse cut file gagg_0_15_"+run_letter+"_Be failed.\n";
		return 1;
	}
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_0_15_"+run_letter,
		"4He",
		false,
		he_cut_0_15
	)) {
		std::cerr << "Error: Parse cut file gagg_0_15_"+run_letter+"_4He failed.\n";
		return 1;
	}
	std::unique_ptr<TCutG> be_cut_16_24, he_cut_16_24;
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_16_24_"+run_letter,
		"Be",
		false,
		be_cut_16_24
	)) {
		std::cerr << "Error: Parse cut file gagg_16_24_"+run_letter+"_Be failed.\n";
		return 1;
	}
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_16_24_"+run_letter,
		"4He",
		false,
		he_cut_16_24
	)) {
		std::cerr << "Error: Parse cut file gagg_16_24_"+run_letter+"_4He failed.\n";
		return 1;
	}

	std::unique_ptr<TCutG> be_cut_18, he_cut_18;
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_18"+run_letter,
		"Be",
		false,
		be_cut_18
	)) {
		std::cerr << "Error: Parse cut file gagg_18_"+run_letter+"_Be failed.\n";
		return 1;
	}
	if (brill::ParseCutFile(
		config.root.workspace,
		"gagg_18"+run_letter,
		"4He",
		false,
		he_cut_18
	)) {
		std::cerr << "Error: Parse cut file gagg_18_"+run_letter+"_4He failed.\n";
		return 1;
	}

	for (int i = 0; i < 16; ++i) {
		slices[i].push_back(StraightSlice("be"+std::to_string(i), be_cut_0_15.get()));
		slices[i].push_back(StraightSlice("he"+std::to_string(i), he_cut_0_15.get()));
	}
	for (int i = 16; i < 25; ++i) {
		if (i == 18) {
			slices[i].push_back(StraightSlice("be18", be_cut_18.get()));
			slices[i].push_back(StraightSlice("he18", he_cut_18.get()));
		}
		slices[i].push_back(StraightSlice("be"+std::to_string(i), be_cut_16_24.get()));
		slices[i].push_back(StraightSlice("he"+std::to_string(i), he_cut_16_24.get()));
	}

	TString output_path = TString::Format(
		"%s/gagg_straight_%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.estimate).c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	std::vector<TH2F> pids, spids;
	for (int i = 0; i < 16; ++i) {
		pids.emplace_back(
			TString::Format("p%d", i),
			TString::Format("PID of GAGG %d", i),
			1000, 0, run_letter == "a" ? 10000 : 6000,
			1000, 0, 80000
		);
		spids.emplace_back(
			TString::Format("s%d", i),
			TString::Format("straight PID of GAGG %d", i),
			1000, 0, run_letter == "a" ? 10000 : 6000,
			1000, 0, 80000
		);
	}
	for (int i = 16; i < 25; ++i) {
		if (i == 18) {
			pids.emplace_back(
				TString::Format("p%d", i),
				TString::Format("PID of GAGG %d", i),
				1000, 0, run_letter == "a" ? 18000 : 2000,
				1000, 0, 80000
			);
			spids.emplace_back(
				TString::Format("s%d", i),
				TString::Format("straight PID of GAGG %d", i),
				1000, 0, run_letter == "a" ? 18000 : 2000,
				1000, 0, 80000
			);
			continue;
		}
		pids.emplace_back(
			TString::Format("p%d", i),
			TString::Format("PID of GAGG %d", i),
			1000, 0, run_letter == "a" ? 60000 : 25000,
			1000, 0, 80000
		);
		spids.emplace_back(
			TString::Format("s%d", i),
			TString::Format("straight PID of GAGG %d", i),
			1000, 0, run_letter == "a" ? 60000 : 25000,
			1000, 0, 80000
		);
	}
	std::vector<TH1F> epids;
	for (int i = 0; i < 25; ++i) {
		epids.emplace_back(
			TString::Format("e%d", i),
			TString::Format("Straight energy of GAGG %d", i),
			1000, 0, 80000
		);
	}


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
				if (!brill::t0::IsInTrackWindow(
					d2_event.front_strip[i],
					d2_event.back_strip[i],
					gagg_event.index[j]
				)) continue;
				for (auto &slice : slices[gagg_event.index[j]]) {
					if (!slice.cut->IsInside(gagg_event.amplitude[j], d2_event.energy[i])) continue;
					slice.graph.AddPoint(gagg_event.amplitude[j], d2_event.energy[i]);
				}
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	for (int i = 0; i < 16; ++i) {
		slices[i][0].a = 0.5;
		slices[i][0].b = -0.05;
		slices[i][1].a = 0.5;
		slices[i][1].b = -0.05;
	}
	for (int i = 16; i < 25; ++i) {
		slices[i][0].a = 0.5;
		slices[i][0].b = -0.05;
		slices[i][1].a = 0.5;
		slices[i][1].b = -0.05;
	}
	slices[18][0].a = 0.5;
	slices[18][0].b = -0.05;
	slices[18][1].a = 0.5;
	slices[18][1].b = -0.05;

	for (int i = 0; i < 25; ++i) {
		if (FitSlice(slices[i][0])) {
			std::cerr << "Error: Failed to fit slice " << i << " (Be).\n";
		}
		if (FitSlice(slices[i][1])) {
			std::cerr << "Error: Failed to fit slice " << i << " (4He).\n";
		}
	}

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
				if (!brill::t0::IsInTrackWindow(
					d2_event.front_strip[i],
					d2_event.back_strip[i],
					gagg_event.index[j]
				)) continue;
				FillGraphs(
					slices[gagg_event.index[j]][use_particle].a,
					slices[gagg_event.index[j]][1].b,
					d2_event.energy[i],
					gagg_event.amplitude[j],
					pids[gagg_event.index[j]],
					spids[gagg_event.index[j]],
					epids[gagg_event.index[j]]
				);
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	for (int i = 0; i < 25; ++i) {
		slices[i][0].graph.Write(TString::Format("gbe%d", i));
		slices[i][1].graph.Write(TString::Format("ghe%d", i));
	}
	for (auto &hist : pids) hist.Write();
	for (auto &hist : spids) hist.Write();
	for (auto &hist : epids) hist.Write();
	opf.Close();

	return 0;
}