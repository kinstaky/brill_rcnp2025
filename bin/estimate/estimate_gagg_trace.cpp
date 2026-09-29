#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraph.h>
#include <TF1.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_pid", "Estimate T0 PID after DSSD matching.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Start run number.", cxxopts::value<int>(), "run")
		("e,end-run", "End run number.", cxxopts::value<int>(), "run")
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

	const std::string ore_dir = brill::JoinPath(config.root.workspace, config.paths.ore);

	TChain chain("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsSkipRun(current_run)) continue;
		++added_runs;
		chain.Add(TString::Format(
			"%s/gagg_fit_%04d.root",
			ore_dir.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}

	short index, time;
	double amplitude, integral, tau;
	chain.SetBranchAddress("index", &index);
	chain.SetBranchAddress("time", &time);
	chain.SetBranchAddress("amplitude", &amplitude);
	chain.SetBranchAddress("integral", &integral);
	chain.SetBranchAddress("tau", &tau);

	TString output_path = TString::Format(
		"%s/gagg_fit_%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.estimate).c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TGraph g_amp_int[16];
	std::unique_ptr<TH1F> hist_residual[16];
	for (int i = 0; i < 16; ++i) {
		hist_residual[i] = std::make_unique<TH1F>(
			TString::Format("h%d", i),
			TString::Format("Residual for GAGG %d", i),
			100, -200, 200
		);
	}

	const long long total = chain.GetEntries();
	long long last_percentage = -1;
	std::printf("Fitting GAGG parameters   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain.GetEntry(entry);
		if (index >= 16) continue;
		if (time < 20 || time > 236) continue;
		if (integral > 1.2e5) continue;
		g_amp_int[index].AddPoint(integral, amplitude);
	}
	std::printf("\b\b\b\b100%%\n");

	// fit
	double p0[16], p1[16];
	for (int i = 0; i < 16; ++i) {
		TF1 *f1 = new TF1(TString::Format("f%d", i), "pol1", 8.2e4, 1.2e5);
		g_amp_int[i].Fit(f1, "QR+ ROB=0.8");
		p0[i] = f1->GetParameter(0);
		p1[i] = f1->GetParameter(1);
	}
	for (int i = 0; i < 16; ++i) std::cout << p0[i] << ",\n"[i == 15];
	for (int i = 0; i < 16; ++i) std::cout << p1[i] << ",\n"[i == 15];

	for (int i = 0; i < 16; ++i) {
		int n = g_amp_int[i].GetN();
		double *x = g_amp_int[i].GetX();
		double *y = g_amp_int[i].GetY();
		for (int j = 0; j < n; ++j) {
			if (x[j] < 8.2e4 || x[j] > 1.2e5) continue;
			double residual = y[j] - (p0[i] + p1[i] * x[j]);
			hist_residual[i]->Fill(residual);
		}
	}

	opf.cd();
	for (int i = 0; i < 16; ++i) g_amp_int[i].Write(TString::Format("g%d", i));
	for (int i = 0; i < 16; ++i) {
		hist_residual[i]->Write();
		hist_residual[i] = nullptr;
	}
	opf.Close();

	return 0;
}
