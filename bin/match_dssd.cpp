#include <iostream>
#include <set>
#include <vector>
#include <memory>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>

#include <Eigen/Dense>

#include "external/cxxopts.hpp"
#include "include/config.h"
#include "include/utils.h"
#include "include/t0/dssd.h"
#include "include/normalize_extra.h"
#include "include/strip_combination.h"
#include "include/t0_utils.h"

bool IsNormalStrip(const brill::EnergyGuess &guess) {
	return guess.type == brill::StripType::Normal
		|| guess.type == brill::StripType::Piecewise;
}

bool IsNormalSharedStrip(const brill::EnergyGuess &guess) {
	return guess.type == brill::StripType::ChargeShared
		|| guess.type == (brill::StripType::Piecewise | brill::StripType::ChargeShared);
}

bool IsLineStrip(const brill::EnergyGuess &guess) {
	return guess.type == brill::StripType::Short
		|| guess.type == brill::StripType::BrokenAdjacent;
}

bool IsPlaneStrip(const brill::EnergyGuess &guess) {
	return guess.type == (brill::StripType::BrokenAdjacent | brill::StripType::ChargeShared);
}

bool IsShortSharedStrip(const brill::EnergyGuess &guess) {
	return guess.type == (brill::StripType::ChargeShared | brill::StripType::Short);
}

void MatchT0D1WithSpecialStrips(
	const brill::SiliconDetectorConfig* detector,
	TTree *ipt,
	brill::DssdEvent &raw,
	const brill::DssdNormalizeParameters &parameters,
	const brill::T0D1ExtraNormalizeParameters &extra,
	const std::map<std::string, brill::FullPCAParameter> &pca_parameters,
	const long long single_entry
) {
	if (single_entry != -1) {
		ipt->GetEntry(single_entry);
		for (int i = 0; i < raw.front_num; ++i) {
			raw.front_energy[i] = raw.front_integral[i];
		}
		for (int i = 0; i < raw.back_num; ++i) {
			raw.back_energy[i] = raw.back_integral[i];
		}

		// search for possible single side combination
		std::vector<std::unique_ptr<brill::StripCombination>> front_comb, back_comb;
		brill::t0::GetT0D1FrontCombinations(
			detector,
			raw,
			parameters,
			extra,
			pca_parameters,
			front_comb
		);
		brill::t0::GetT0D1BackCombinations(
			detector,
			raw,
			parameters,
			pca_parameters,
			back_comb
		);
		// match
		for (size_t i = 0; i < front_comb.size(); ++i) {
			for (size_t j = 0; j < back_comb.size(); ++j) {
				std::cout << "Trying " << i << ", " << j << std::endl;
				brill::MatchResult result = front_comb[i]->Match(
					*back_comb[j],
					detector->match_tolerance,
					true
				);
				if (!result.valid) continue;
				std::cout << "Valid result: "
					<< result.energy << ", "
					<< result.distance << ", "
					<< result.flag << std::endl;
			}
		}
		return;
	}
	// estimate guess distance
	TH1F hist_distance("hd", "Distance of strips", 500, 0, 5000);
	TH1F hist_distance_nn("hdnn", "Distance of normal-normal strips", 500, 0, 5000);
	TH1F hist_distance_nc("hdnc", "Distance of normal-shared strips", 500, 0, 5000);
	TH1F hist_distance_cc("hdcc", "Distance of shared-shared strips", 500, 0, 5000);
	TH1F hist_distance_np("hdnp", "Distance of normal-plane strips", 500, 0, 5000);
	TH1F hist_distance_nl("hdnl", "Distance of normal-line strips", 500, 0, 5000);
	TH1F hist_distance_cp("hdcp", "Distance of shared-plane strips", 500, 0, 5000);
	TH1F hist_distance_cl("hdcl", "Distance of shared-line strips", 500, 0, 5000);
	TH1F hist_distance_nsc("hdnsc", "Distance of normal-short_shared strips", 500, 0, 5000);
	TH1F hist_distance_csc("hdcsc", "Distance of shared-short_shared strips", 500, 0, 5000);


	// residual events
	TH1F hist_front_num("hrfn", "Front residual hit distribution", 8, 0, 8);
	TH1F hist_back_num("hrbn", "Back residual hit distribution", 8, 0, 8);
	TH1F hist_front_strip("hrfs", "Front residual strip distribution", 128, 0, 128);
	TH1F hist_back_strip("hrbs", "Back residual strip distribution", 128, 0, 128);
	TH2F hist_front_strip_pair("hrfsp", "Front residual strip pair distribution", 128, 0, 128, 128, 0, 128);
	TH2F hist_back_strip_pair("hrbsp", "Back residual strip pair distribution", 128, 0, 128, 128, 0, 128);

	// data tree
	TTree opt("tree", "Match T0D1");
	brill::DssdMatchEvent match;
	int special_flag[8] = {0};
	brill::SetupOutput(&opt, match);
	opt.Branch("special_flag", special_flag, "sf[num]/I");

	// residual tree
	TTree rtree("rtree", "Residual T0D1 events");
	brill::DssdEvent residual_event;
	bool front_has_strip[128], back_has_strip[128];
	brill::SetupOutput(&rtree, residual_event);
	rtree.Branch("front_has_strip", front_has_strip, "fhs[128]/O");
	rtree.Branch("back_has_strip", back_has_strip, "bhs[128]/O");

	// loop events
	long long total = ipt->GetEntries();
	long long last_percentage = 0;
	printf("Matching   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		for (int i = 0; i < raw.front_num; ++i) {
			raw.front_energy[i] = raw.front_integral[i];
		}
		for (int i = 0; i < raw.back_num; ++i) {
			raw.back_energy[i] = raw.back_integral[i];
		}

		// search for possible single side combination
		std::vector<std::unique_ptr<brill::StripCombination>> front_comb, back_comb;
		brill::t0::GetT0D1FrontCombinations(
			detector,
			raw,
			parameters,
			extra,
			pca_parameters,
			front_comb
		);
		brill::t0::GetT0D1BackCombinations(
			detector,
			raw,
			parameters,
			pca_parameters,
			back_comb
		);
		// estimate guess distances
		for (auto &fc : front_comb) {
			for (auto &bc : back_comb) {
				for (const auto &fg : fc->Guesses()) {
					for (const auto &bg : bc->Guesses()) {
						double distance = fabs(fg.value - bg.value);
						hist_distance.Fill(distance);
						if (IsNormalStrip(fg) && IsNormalStrip(bg)) {
							hist_distance_nn.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsNormalSharedStrip(bg))
							|| (IsNormalSharedStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nc.Fill(distance);
						} else if (IsNormalSharedStrip(fg) && IsNormalSharedStrip(bg)) {
							hist_distance_cc.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsLineStrip(bg))
							|| (IsLineStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nl.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsPlaneStrip(bg))
							|| (IsPlaneStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_np.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsLineStrip(bg))
							|| (IsLineStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_cl.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsPlaneStrip(bg))
							|| (IsPlaneStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_cp.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsShortSharedStrip(bg))
							|| (IsShortSharedStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nsc.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsShortSharedStrip(bg))
							|| (IsShortSharedStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_csc.Fill(distance);
						}
					}
				}
			}
		}

		// match
		std::vector<brill::MatchResult> match_results;
		for (size_t i = 0; i < front_comb.size(); ++i) {
			for (size_t j = 0; j < back_comb.size(); ++j) {
				brill::MatchResult result = front_comb[i]->Match(
					*back_comb[j],
					detector->match_tolerance,
					false
				);
				if (!result.valid) continue;
				match_results.push_back(result);
			}
		}

		// sort
		std::sort(
			match_results.begin(),
			match_results.end(),
			[](const brill::MatchResult &a, const brill::MatchResult &b) {
				return a.distance < b.distance;
			}
		);

		// pick
		int used_flag = 0;
		match.num = 0;
		for (auto &m : match_results) {
			if (used_flag & m.flag) continue;
			used_flag |= m.flag;
			match.flag[match.num] = m.flag;
			special_flag[match.num] = 0;
			special_flag[match.num] |= (m.type[0] > 1 ? 1 : 0);
			special_flag[match.num] |= (m.type[1] > 1 ? 2 : 0);
			match.front_strip[match.num] = m.strip[0];
			match.back_strip[match.num] = m.strip[1];
			match.energy[match.num] = m.energy;
			match.time[match.num] = m.time;
			brill::t0::GetPixelPosition(
				detector,
				m.strip[0],
				m.strip[1],
				match.x[match.num],
				match.y[match.num],
				match.z[match.num]
			);
			++match.num;
		}
		match.run = raw.run;
		match.entry = int(entry);

		// have a look at the distribution of residual events
		residual_event.front_num = 0;
		residual_event.back_num = 0;
		memset(front_has_strip, 0, sizeof(front_has_strip));
		memset(back_has_strip, 0, sizeof(back_has_strip));
		for (int i = 0; i < raw.front_num; ++i) {
			if (used_flag & (0x1 << i)) continue;
			residual_event.front_strip[residual_event.front_num] = raw.front_strip[i];
			residual_event.front_energy[residual_event.front_num] = parameters.NormEnergy(
				0, raw.front_strip[i], raw.front_integral[i]
			);
			front_has_strip[raw.front_strip[i]] = true;
			++residual_event.front_num;
		}
		for (int i = 0; i < raw.back_num; ++i) {
			if (used_flag & (0x100 << i)) continue;
			residual_event.back_strip[residual_event.back_num] = raw.back_strip[i];
			residual_event.back_energy[residual_event.back_num] = parameters.NormEnergy(
				1, raw.back_strip[i], raw.back_integral[i]
			);
			back_has_strip[raw.back_strip[i]] = true;
			++residual_event.back_num;
		}
		// fill to histogram
		if (residual_event.front_num) hist_front_num.Fill(residual_event.front_num);
		for (int i = 0; i < residual_event.front_num; ++i) {
			hist_front_strip.Fill(residual_event.front_strip[i]);
			for (int j = i+1; j < residual_event.front_num; ++j) {
				hist_front_strip_pair.Fill(
					residual_event.front_strip[i],
					residual_event.front_strip[j]
				);
			}
		}
		if (residual_event.back_num) hist_back_num.Fill(residual_event.back_num);
		for (int i = 0; i < residual_event.back_num; ++i) {
			hist_back_strip.Fill(residual_event.back_strip[i]);
			for (int j = i+1; j < residual_event.back_num; ++j) {
				hist_back_strip_pair.Fill(
					residual_event.back_strip[i],
					residual_event.back_strip[j]
				);
			}
		}
		residual_event.run = raw.run;
		residual_event.entry = raw.entry;

		opt.Fill();
		rtree.Fill();
	}
	printf("\b\b\b\b100%%\n");

	// save
	hist_distance.Write();
	hist_distance_nn.Write();
	hist_distance_nc.Write();
	hist_distance_cc.Write();
	hist_distance_nl.Write();
	hist_distance_np.Write();
	hist_distance_cl.Write();
	hist_distance_cp.Write();
	hist_distance_nsc.Write();
	hist_distance_csc.Write();
	hist_front_num.Write();
	hist_front_strip.Write();
	hist_front_strip_pair.Write();
	hist_back_num.Write();
	hist_back_strip.Write();
	hist_back_strip_pair.Write();
	opt.Write();
	rtree.Write();
}

void MatchT0D2WithSpecialStrips(
	const brill::SiliconDetectorConfig* detector,
	TTree *ipt,
	brill::DssdEvent &raw,
	const brill::DssdNormalizeParameters &parameters,
	const brill::T0D2ExtraNormalizeParameters &extra,
	const std::map<std::string, brill::FullPCAParameter> &pca_parameters,
	const long long single_entry
) {
	if (single_entry != -1) {
		ipt->GetEntry(single_entry);
		for (int i = 0; i < raw.front_num; ++i) {
			raw.front_energy[i] = raw.front_integral[i];
		}
		for (int i = 0; i < raw.back_num; ++i) {
			raw.back_energy[i] = raw.back_integral[i];
		}

		// search for possible single side combination
		std::vector<std::unique_ptr<brill::StripCombination>> front_comb, back_comb;
		brill::t0::GetT0D2FrontCombinations(
			detector,
			raw,
			parameters,
			extra,
			pca_parameters,
			front_comb
		);
		brill::t0::GetT0D2BackCombinations(
			detector,
			raw,
			parameters,
			pca_parameters,
			back_comb
		);

		// match
		for (size_t i = 0; i < front_comb.size(); ++i) {
			for (size_t j = 0; j < back_comb.size(); ++j) {
				std::cout << "Trying " << i << ", " << j << std::endl;

				brill::MatchResult result = front_comb[i]->Match(
					*back_comb[j],
					detector->match_tolerance,
					true
				);
				if (!result.valid) continue;
				std::cout << "Valid result: "
					<< result.energy << ", "
					<< result.distance << ", "
					<< result.flag << std::endl;
			}
		}
		return;
	}

	// estimate guess distance
	TH1F hist_distance("hd", "Distance of strips", 500, 0, 5000);
	TH1F hist_distance_nn("hdnn", "Distance of normal-normal strips", 500, 0, 5000);
	TH1F hist_distance_nc("hdnc", "Distance of normal-shared strips", 500, 0, 5000);
	TH1F hist_distance_cc("hdcc", "Distance of shared-shared strips", 500, 0, 5000);
	TH1F hist_distance_np("hdnp", "Distance of normal-plane strips", 500, 0, 5000);
	TH1F hist_distance_nl("hdnl", "Distance of normal-line strips", 500, 0, 5000);
	TH1F hist_distance_cp("hdcp", "Distance of shared-plane strips", 500, 0, 5000);
	TH1F hist_distance_cl("hdcl", "Distance of shared-line strips", 500, 0, 5000);
	TH1F hist_distance_nsc("hdnsc", "Distance of normal-short_shared strips", 500, 0, 5000);
	TH1F hist_distance_csc("hdcsc", "Distance of shared-short_shared strips", 500, 0, 5000);


	// residual events
	TH1F hist_front_num("hrfn", "Front residual hit distribution", 8, 0, 8);
	TH1F hist_back_num("hrbn", "Back residual hit distribution", 8, 0, 8);
	TH1F hist_front_strip("hrfs", "Front residual strip distribution", 128, 0, 128);
	TH1F hist_back_strip("hrbs", "Back residual strip distribution", 128, 0, 128);
	TH2F hist_front_strip_pair("hrfsp", "Front residual strip pair distribution", 128, 0, 128, 128, 0, 128);
	TH2F hist_back_strip_pair("hrbsp", "Back residual strip pair distribution", 128, 0, 128, 128, 0, 128);

	// data tree
	TTree opt("tree", "Match T0D2");
	brill::DssdMatchEvent match;
	int special_flag[8] = {0};
	brill::SetupOutput(&opt, match);
	opt.Branch("special_flag", special_flag, "sf[num]/I");

	// residual tree
	TTree rtree("rtree", "Residual T0D2 events");
	brill::DssdEvent residual_event;
	bool front_has_strip[128], back_has_strip[128];
	brill::SetupOutput(&rtree, residual_event);
	rtree.Branch("front_has_strip", front_has_strip, "fhs[128]/O");
	rtree.Branch("back_has_strip", back_has_strip, "bhs[128]/O");

	// loop events
	long long total = ipt->GetEntries();
	long long last_percentage = 0;
	printf("Matching   0%%");
	fflush(stdout);
	for (long long entry = 0; entry < total; ++entry) {
		if (entry * 100 / total > last_percentage) {
			last_percentage = entry * 100 / total;
			printf("\b\b\b\b%3lld%%", last_percentage);
			fflush(stdout);
		}
		ipt->GetEntry(entry);
		for (int i = 0; i < raw.front_num; ++i) {
			raw.front_energy[i] = raw.front_integral[i];
		}
		for (int i = 0; i < raw.back_num; ++i) {
			raw.back_energy[i] = raw.back_integral[i];
		}

		// search for possible single side combination
		std::vector<std::unique_ptr<brill::StripCombination>> front_comb, back_comb;
		brill::t0::GetT0D2FrontCombinations(
			detector,
			raw,
			parameters,
			extra,
			pca_parameters,
			front_comb
		);
		brill::t0::GetT0D2BackCombinations(
			detector,
			raw,
			parameters,
			pca_parameters,
			back_comb
		);
		// estimate guess distances
		for (auto &fc : front_comb) {
			for (auto &bc : back_comb) {
				for (const auto &fg : fc->Guesses()) {
					for (const auto &bg : bc->Guesses()) {
						double distance = fabs(fg.value - bg.value);
						hist_distance.Fill(distance);
						if (IsNormalStrip(fg) && IsNormalStrip(bg)) {
							hist_distance_nn.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsNormalSharedStrip(bg))
							|| (IsNormalSharedStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nc.Fill(distance);
						} else if (IsNormalSharedStrip(fg) && IsNormalSharedStrip(bg)) {
							hist_distance_cc.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsLineStrip(bg))
							|| (IsLineStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nl.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsPlaneStrip(bg))
							|| (IsPlaneStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_np.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsLineStrip(bg))
							|| (IsLineStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_cl.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsPlaneStrip(bg))
							|| (IsPlaneStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_cp.Fill(distance);
						} else if (
							(IsNormalStrip(fg) && IsShortSharedStrip(bg))
							|| (IsShortSharedStrip(fg) && IsNormalStrip(bg))
						) {
							hist_distance_nsc.Fill(distance);
						} else if (
							(IsNormalSharedStrip(fg) && IsShortSharedStrip(bg))
							|| (IsShortSharedStrip(fg) && IsNormalSharedStrip(bg))
						) {
							hist_distance_csc.Fill(distance);
						}
					}
				}
			}
		}

		// match
		std::vector<brill::MatchResult> match_results;
		for (size_t i = 0; i < front_comb.size(); ++i) {
			for (size_t j = 0; j < back_comb.size(); ++j) {
				brill::MatchResult result = front_comb[i]->Match(
					*back_comb[j],
					detector->match_tolerance,
					false
				);
				if (!result.valid) continue;
				match_results.push_back(result);
			}
		}

		// sort
		std::sort(
			match_results.begin(),
			match_results.end(),
			[](const brill::MatchResult &a, const brill::MatchResult &b) {
				return a.distance < b.distance;
			}
		);

		// pick
		int used_flag = 0;
		match.num = 0;
		for (auto &m : match_results) {
			if (used_flag & m.flag) continue;
			used_flag |= m.flag;
			match.flag[match.num] = m.flag;
			special_flag[match.num] = 0;
			special_flag[match.num] |= (m.type[0] > 1 ? 1 : 0);
			special_flag[match.num] |= (m.type[1] > 1 ? 2 : 0);
			match.front_strip[match.num] = m.strip[0];
			match.back_strip[match.num] = m.strip[1];
			match.energy[match.num] = m.energy;
			match.time[match.num] = m.time;
			brill::t0::GetPixelPosition(
				detector,
				m.strip[0],
				m.strip[1],
				match.x[match.num],
				match.y[match.num],
				match.z[match.num]
			);
			++match.num;
		}
		match.run = raw.run;
		match.entry = raw.entry;

		// have a look at the distribution of residual events
		residual_event.front_num = 0;
		residual_event.back_num = 0;
		memset(front_has_strip, 0, sizeof(front_has_strip));
		memset(back_has_strip, 0, sizeof(back_has_strip));
		for (int i = 0; i < raw.front_num; ++i) {
			if (used_flag & (0x1 << i)) continue;
			residual_event.front_strip[residual_event.front_num] = raw.front_strip[i];
			residual_event.front_energy[residual_event.front_num] = parameters.NormEnergy(
				0, raw.front_strip[i], raw.front_integral[i]
			);
			front_has_strip[raw.front_strip[i]] = true;
			++residual_event.front_num;
		}
		for (int i = 0; i < raw.back_num; ++i) {
			if (used_flag & (0x100 << i)) continue;
			residual_event.back_strip[residual_event.back_num] = raw.back_strip[i];
			residual_event.back_energy[residual_event.back_num] = parameters.NormEnergy(
				1, raw.back_strip[i], raw.back_integral[i]
			);
			back_has_strip[raw.back_strip[i]] = true;
			++residual_event.back_num;
		}
		// fill to histogram
		if (residual_event.front_num) hist_front_num.Fill(residual_event.front_num);
		for (int i = 0; i < residual_event.front_num; ++i) {
			hist_front_strip.Fill(residual_event.front_strip[i]);
			for (int j = i+1; j < residual_event.front_num; ++j) {
				hist_front_strip_pair.Fill(
					residual_event.front_strip[i],
					residual_event.front_strip[j]
				);
			}
		}
		if (residual_event.back_num) hist_back_num.Fill(residual_event.back_num);
		for (int i = 0; i < residual_event.back_num; ++i) {
			hist_back_strip.Fill(residual_event.back_strip[i]);
			for (int j = i+1; j < residual_event.back_num; ++j) {
				hist_back_strip_pair.Fill(
					residual_event.back_strip[i],
					residual_event.back_strip[j]
				);
			}
		}
		residual_event.run = raw.run;
		residual_event.entry = raw.entry;

		opt.Fill();
		rtree.Fill();
	}
	printf("\b\b\b\b100%%\n");

	// save
	hist_distance.Write();
	hist_distance_nn.Write();
	hist_distance_nc.Write();
	hist_distance_cc.Write();
	hist_distance_nl.Write();
	hist_distance_np.Write();
	hist_distance_cl.Write();
	hist_distance_cp.Write();
	hist_distance_nsc.Write();
	hist_distance_csc.Write();
	hist_front_num.Write();
	hist_front_strip.Write();
	hist_front_strip_pair.Write();
	hist_back_num.Write();
	hist_back_strip.Write();
	hist_back_strip_pair.Write();
	opt.Write();
	rtree.Write();
}

int main(int argc, char **argv) {
	// Parse arguments
	cxxopts::Options options("match", "Match DSSD front-back correlated events.");
	options.add_options()
		("h,help", "Print help information.")
		("r,run", "Run number.", cxxopts::value<int>(), "run")
		(
			"c,config",
			"Config file path.",
			cxxopts::value<std::string>()->default_value("config.toml"),
			"file"
		)
		(
			"s,single",
			"Single entry for debug.",
			cxxopts::value<long long>(),
			"entry"
		)
		(
			"detector",
			"Detectors to match.",
			cxxopts::value<std::vector<std::string>>(),
			"detector"
		);
	options.parse_positional({"detector"});

	auto result = options.parse(argc, argv);
	if (result.count("help")) {
		std::cout << options.help() << std::endl;
		return 0;
	}
	if (!result.count("run")) {
		std::cerr << "Error: Missing required option --run.\n";
		std::cout << options.help() << std::endl;
		return -1;
	}
	if (result.count("run") > 1) {
		std::cerr << "Error: Multiple --run options, only single run is allowed.\n";
		std::cout << options.help() << std::endl;
		return -1;
	}
	if (!result.count("detector")) {
		std::cerr << "Error: Missing detector positional argument.\n";
		std::cout << options.help() << std::endl;
		return -2;
	}
	long long single_entry = result.count("single")
		? result["single"].as<long long>()
		: -1;

	// check detectors
	const std::set<std::string> allowed_detectors = {"t0d1", "t0d2"};
	std::vector<std::string> detectors = result["detector"].as<std::vector<std::string>>();
	for (const auto &detector : detectors) {
		if (allowed_detectors.count(detector) == 0) {
			std::cerr << "Warning: Unsupported detector " << detector << ".\n";
			std::cout << "  Supported detectors: t0d1, t0d2." << std::endl;
			return -2;
		}
	}

	// load config
	brill::AppConfig config;
	if (config.Load(result["config"].as<std::string>())) {
		std::cerr << "Error: Load config failed." << std::endl;
		return -3;
	}
	const int run = result["run"].as<int>();
	if (config.IsSkipRun(run)) {
		std::cerr << "Info: Run " << run << " is jump run." << std::endl;
		return -1;
	}

	// names
	const std::map<std::string, std::vector<std::string>> strip_names = {
		{"t0d1", {"fs50", "fs94", "fs120"}},
		{"t0d2", {"fs2", "fs1819", "fs70", "bs32a", "bs32b", "bs32c", "bs36", "bs48", "bs68"}},
	};

	// apply normalize
	for (const auto &detector : detectors) {
		const brill::SiliconDetectorConfig *detector_config =
			config.FindDetector(detector);
		if (!detector_config) {
			std::cerr << "Error: Find detector " << detector << " failed.\n";
			return -2;
		}

		// get normalize run
		int normalize_run = run;
		for (const auto &[start, use] : config.normalize.runs) {
			if (run >= start) normalize_run = use;
		}

		// prepare paths
		std::string normalize_dir = brill::JoinPath(config.root.workspace, config.paths.normalize);
		std::string match_dir = brill::JoinPath(config.root.workspace, config.paths.match);

		// setup input
		TString input_filename = TString::Format(
			"%s/%s_%04d.root",
			brill::JoinPath(config.root.workspace, config.paths.ingot).c_str(),
			detector.c_str(),
			run
		);
		TFile ipf(input_filename, "read");
		TTree *ipt = (TTree*)ipf.Get("tree");
		if (!ipt) {
			std::cerr << "Error: Get tree from " << input_filename << " failed.\n";
			return -4;
		}
		brill::DssdEvent raw_event;
		brill::SetupInput(ipt, raw_event);

		// read parameters
		brill::DssdNormalizeParameters parameters(
			detector_config->front_strips, detector_config->back_strips
		);
		std::string parameter_path = TString::Format(
			"%s/%s_%04d.txt",
			normalize_dir.c_str(),
			detector.c_str(),
			normalize_run
		).Data();
		if (parameters.Read(parameter_path)) {
			std::cerr << "Error: Read " << detector << " normalize parameters failed.\n";
			return -4;
		}

		// read extra parameters
		std::string extra_parameter_path = TString::Format(
			"%s/%s_extra_%04d.txt",
			normalize_dir.c_str(),
			detector.c_str(),
			normalize_run
		).Data();

		TFile *opf = nullptr;
		if (single_entry == -1) {
			// setup output
			TString output_filename = TString::Format(
				"%s/%s_%04d.root",
				match_dir.c_str(),
				detector.c_str(),
				run
			);
			opf = new TFile(output_filename, "recreate");
		}

		if (detector == "t0d1") {
			brill::T0D1ExtraNormalizeParameters extra_parameters;
			if (extra_parameters.Read(extra_parameter_path)) {
				std::cerr << "Error: Read T0D1 extra normalize parameters failed.\n";
				return -4;
			}
			std::map<std::string, brill::FullPCAParameter> pca_parameters = {
				{"fs50", {0, {49, 51}, 50, false, {}, {{}, {}}, {400.0, 400.0, 400.0}}},
				{"fs94", {0, {93, 95}, 94, false, {}, {{}, {}}, {600.0, 400.0, 400.0}}},
				{"fs120", {0, {119, 121}, 120, false, {}, {{}, {}}, {400.0, 400.0, 400.0}}},
			};
			for (auto &[name, pca] : pca_parameters) {
				const auto &line_pca = extra_parameters.pca.at(name);
				const auto &plane0_pca = extra_parameters.pca.at(name + "p0");
				const auto &plane1_pca = extra_parameters.pca.at(name + "p1");
				pca.line = brill::Line3D {line_pca.mean, line_pca.direction};
				pca.plane[0] = brill::Plane {plane0_pca.mean, plane0_pca.direction};
				pca.plane[1] = brill::Plane {plane1_pca.mean, plane1_pca.direction};
			}
			const auto &bs104_pca = extra_parameters.pca.at("bs104");
			pca_parameters.emplace("bs104", brill::FullPCAParameter {
				1, {104, 109}, 104, false,
				{bs104_pca.mean, bs104_pca.direction},
				{{}, {}},
				{500.0, 0.0, 0.0}
			});
			MatchT0D1WithSpecialStrips(
				detector_config,
				ipt,
				raw_event,
				parameters,
				extra_parameters,
				pca_parameters,
				single_entry
			);
		} else if (detector == "t0d2") {
			brill::T0D2ExtraNormalizeParameters extra_parameters;
			if (extra_parameters.Read(extra_parameter_path)) {
				std::cerr << "Error: Read T0D2 extra normalize parameters failed.\n";
				return -4;
			}
			std::map<std::string, brill::FullPCAParameter> pca_parameters = {
				{"fs2", {0, {1, 3}, 2, false, {}, {{}, {}}, {200.0, 600.0, 600.0}}},
				{"fs1819", {0, {17, 20}, 18, false, {}, {{}, {}}, {300.0, 300.0, 200.0}}},
				{"fs70", {0, {69, 71}, 70, false, {}, {{}, {}}, {300.0, 300.0, 300.0}}},
				{"bs32a", {1, {31, 33}, 32, true, {}, {{}, {}}, {400.0, 500.0, 600.0}}},
				{"bs32b", {1, {31, 34}, 32, true, {}, {{}, {}}, {500.0, 500.0, 500.0}}},
				{"bs32c", {1, {34, 31}, 32, true, {}, {{}, {}}, {500.0, 500.0, 600.0}}},
				{"bs36", {1, {35, 37}, 36, false, {}, {{}, {}}, {600.0, 600.0, 600.0}}},
				{"bs48", {1, {47, 49}, 48, false, {}, {{}, {}}, {600.0, 700.0, 700.0}}},
				{"bs68", {1, {67, 69}, 68, false, {}, {{}, {}}, {500.0, 600.0, 600.0}}},
			};
			for (auto &[name, pca] : pca_parameters) {
				const auto &line_pca = extra_parameters.pca.at(name);
				const auto &plane0_pca = extra_parameters.pca.at(name + "p0");
				const auto &plane1_pca = extra_parameters.pca.at(name + "p1");
				pca.line = brill::Line3D {line_pca.mean, line_pca.direction};
				pca.plane[0] = brill::Plane {plane0_pca.mean, plane0_pca.direction};
				pca.plane[1] = brill::Plane {plane1_pca.mean, plane1_pca.direction};
			}
			MatchT0D2WithSpecialStrips(
				detector_config,
				ipt,
				raw_event,
				parameters,
				extra_parameters,
				pca_parameters,
				single_entry
			);
		}

		// close files
		ipf.Close();
		if (single_entry == -1) {
			opf->Close();
			delete opf;
		}
	}

	return 0;
}
