#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TString.h>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/event/ingot/silicon_event.h"
#include "include/event/t0/dssd_match_event.h"
#include "include/utils.h"
#include "include/t0_utils.h"

void PrintUsage(const cxxopts::Options &options) {
	std::cout << options.help() << "\n";
}

int main(int argc, char **argv) {
	cxxopts::Options options("estimate_t0_pid", "Estimate T0 PID after DSSD matching.");
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

	const std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);
	const std::string ingot_dir = brill::JoinPath(config.root.workspace, config.paths.ingot);
	const std::string trigger_infix = brill::TriggerInfix(config.root.trigger);

	TChain chain1("tree");
	TChain chain2("tree");
	int added_runs = 0;
	for (int current_run = run; current_run <= end_run; ++current_run) {
		if (config.IsJumpRun(current_run)) continue;
		++added_runs;
		chain1.Add(TString::Format(
			"%s/t0d1_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
		chain2.Add(TString::Format(
			"%s/t0d2_%s%04d.root",
			match_dir.c_str(),
			trigger_infix.c_str(),
			current_run
		));
	}
	if (added_runs == 0) {
		std::cout << "No runs to process after applying jump_run.\n";
		return 0;
	}
	chain1.AddFriend(&chain2, "d2");

	brill::DssdMatchEvent event1;
	int special_flag1[8];
	brill::DssdMatchEvent event2;
	int special_flag2[8];
	brill::SetupInput(&chain1, event1);
	chain1.SetBranchAddress("special_flag", special_flag1);
	brill::SetupInput(&chain1, event2, "d2.");
	chain1.SetBranchAddress("d2.special_flag", special_flag2);

	TString output_path = TString::Format(
		"%s/t0_pid_%s%04d_%04d.root",
		brill::JoinPath(config.root.workspace, config.paths.estimate).c_str(),
		trigger_infix.c_str(),
		run,
		end_run
	);
	TFile opf(output_path, "recreate");
	TH2F d1d2_pid("d1d2", "D1-D2 PID", 5000, 0.0, 80000.0, 5000, 0.0, 60000.0);
	TH2F d1_strip("d1s", "D1 strips", 128, 0, 128, 128, 0, 128);
	TH2F d2_strip("d2s", "D2 strips", 128, 0, 128, 128, 0, 128);
	TH2F d1d2_pid_es("d1d2es", "D1-D2 PID exclude special strips", 5000, 0.0, 80000.0, 5000, 0.0, 60000.0);
	TH2F d1_strip_es("d1ses", "D1 strips exclude special strips", 128, 0, 128, 128, 0, 128);
	TH2F d2_strip_es("d2ses", "D2 strips exclude special strips", 128, 0, 128, 128, 0, 128);

	const long long total = chain1.GetEntries();
	long long last_percentage = -1;
	std::printf("Filling T0 PID   0%%");
	std::fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		long long percentage = total > 0 ? entry * 100ll / total : 100ll;
		if (percentage > last_percentage) {
			last_percentage = percentage;
			std::printf("\b\b\b\b%3lld%%", percentage);
			std::fflush(stdout);
		}
		chain1.GetEntry(entry);
		for (int i = 0; i < event1.num; ++i) {
			for (int j = 0; j < event2.num; ++j) {
				int d1fs = event1.front_strip[i];
				int d1bs = event1.back_strip[i];
				int d2fs = event2.front_strip[j];
				int d2bs = event2.back_strip[j];
				if (brill::t0::IsInTrackWindow(d1fs, d1bs, d2fs, d2bs)) {
					d1d2_pid.Fill(event2.energy[j], event1.energy[i]);
					d1_strip.Fill(d1bs, d1fs);
					d2_strip.Fill(d2bs, d2fs);
					if (!special_flag1[i] && !special_flag2[j]) {
						d1d2_pid_es.Fill(event2.energy[j], event1.energy[i]);
						d1_strip_es.Fill(event1.back_strip[i], event1.front_strip[i]);
						d2_strip_es.Fill(event2.back_strip[j], event2.front_strip[j]);
					}
				}
			}
		}
	}
	std::printf("\b\b\b\b100%%\n");

	opf.cd();
	d1d2_pid.Write();
	d1_strip.Write();
	d2_strip.Write();
	d1d2_pid_es.Write();
	d1_strip_es.Write();
	d2_strip_es.Write();
	opf.Close();

	return 0;
}
