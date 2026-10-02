#include "include/t0_utils.h"

#include <cmath>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace brill::t0 {

void GetT0D1FrontCombinations(
	const SiliconDetectorConfig *detector,
	const DssdEvent &raw,
	const DssdNormalizeParameters &parameters,
	const T0D1ExtraNormalizeParameters &extra,
	const std::map<std::string, FullPCAParameter> &pca_parameters,
	std::vector<std::unique_ptr<StripCombination>> &combinations
) {
	combinations.clear();
	for (int i = 0; i < raw.front_num; ++i) {
		brill::StripInfo strip1 {
			raw.front_strip[i],
			raw.front_energy[i],
			raw.front_time[i],
			parameters.front_p0[raw.front_strip[i]],
			parameters.front_p1[raw.front_strip[i]],
			parameters.front_p2[raw.front_strip[i]],
			parameters.front_p3[raw.front_strip[i]]
		};
		// piecewise
		if (raw.front_strip[i] == 101) {
			combinations.push_back(
				std::make_unique<brill::PiecewiseStripCombination>(
					1 << i,
					strip1,
					extra.piecewise,
					detector->match_tolerance
				)
			);
		} else {
			// normal
			combinations.push_back(
				std::make_unique<brill::NormalStripCombination>(
					1 << i,
					strip1,
					detector->match_tolerance
				)
			);
		}
		for (int j = i+1; j < raw.front_num; ++j) {
			brill::StripInfo strip2 {
				raw.front_strip[j],
				raw.front_energy[j],
				raw.front_time[j],
				parameters.front_p0[raw.front_strip[j]],
				parameters.front_p1[raw.front_strip[j]],
				parameters.front_p2[raw.front_strip[j]],
				parameters.front_p3[raw.front_strip[j]]
			};
			if (abs(raw.front_strip[i] - raw.front_strip[j]) == 1) {
				// piecewise shared
				if (raw.front_strip[i] == 101 || raw.front_strip[j] == 101) {
					combinations.push_back(
						std::make_unique<brill::PiecewiseSharedStripCombination>(
							(1 << i) | (1 << j),
							raw.front_strip[i] == 101 ? strip1 : strip2,
							raw.front_strip[i] == 101 ? strip2 : strip1,
							extra.piecewise,
							detector->match_tolerance
						)
					);
				} else {
					// normal shared
					combinations.push_back(
						std::make_unique<brill::NormalSharedStripCombination>(
							(1 << i) | (1 << j),
							strip1,
							strip2,
							detector->match_tolerance
						)
					);
				}
				continue;
			}
			// broken adjacent
			for (const auto &[name, pca] : pca_parameters) {
				if (pca.side == 1) continue;
				if (!AreStrips(
					raw.front_strip[i], raw.front_strip[j],
					pca.strips[0], pca.strips[1],
					pca.has_order
				)) continue;
				combinations.push_back(
					std::make_unique<brill::BrokenAdjacentStripCombination>(
						(1 << i) | (1 << j),
						raw.front_strip[i] == pca.strips[0] ? strip1 : strip2,
						raw.front_strip[i] == pca.strips[0] ? strip2 : strip1,
						pca.broken_strip,
						pca.line,
						pca.plane[0],
						pca.plane[1],
						pca.threshold[0],
						pca.threshold[1],
						pca.threshold[2]
					)
				);
			}
		}
	}
}

/// @brief Get back strips combination
/// @param[in] raw raw event
/// @param[out] combinations found combinations
void GetT0D1BackCombinations(
	const SiliconDetectorConfig *detector,
	const DssdEvent &raw,
	const DssdNormalizeParameters &parameters,
	const std::map<std::string, FullPCAParameter> &pca_parameters,
	std::vector<std::unique_ptr<StripCombination>> &combinations
) {
	combinations.clear();
	// short strip
	const FullPCAParameter &short_pca = pca_parameters.at("bs104");
	for (int i = 0; i < raw.back_num; ++i) {
		brill::StripInfo strip1 {
			raw.back_strip[i],
			raw.back_energy[i],
			raw.back_time[i],
			parameters.back_p0[raw.back_strip[i]],
			parameters.back_p1[raw.back_strip[i]],
			parameters.back_p2[raw.back_strip[i]],
			parameters.back_p3[raw.back_strip[i]]
		};
		// normal
		if (raw.back_strip[i] != 104 && raw.back_strip[i] != 109) {
			combinations.push_back(
				std::make_unique<brill::NormalStripCombination>(
					0x100 << i,
					strip1,
					detector->match_tolerance
				)
			);
		}
		for (int j = i+1; j < raw.back_num; ++j) {
			brill::StripInfo strip2 {
				raw.back_strip[j],
				raw.back_energy[j],
				raw.back_time[j],
				parameters.back_p0[raw.back_strip[j]],
				parameters.back_p1[raw.back_strip[j]],
				parameters.back_p2[raw.back_strip[j]],
				parameters.back_p3[raw.back_strip[j]]
			};
			if (abs(raw.back_strip[i] - raw.back_strip[j]) == 1) {
				// normal shared
				combinations.push_back(
					std::make_unique<brill::NormalSharedStripCombination>(
						(0x100 << i) | (0x100 << j),
						strip1,
						strip2,
						detector->match_tolerance
					)
				);
			}
			if (AreStrips(
				raw.back_strip[i], raw.back_strip[j],
				short_pca.strips[0], short_pca.strips[1],
				short_pca.has_order
			)) {
				combinations.push_back(
					std::make_unique<brill::ShortStripCombination>(
						(0x100 << i) | (0x100 << j),
						raw.back_strip[i] == short_pca.strips[0] ? strip1 : strip2,
						raw.back_strip[i] == short_pca.strips[1] ? strip2 : strip1,
						short_pca.line,
						short_pca.threshold[0]
					)
				);
			}
			for (int k = j+1; k < raw.back_num; ++k) {
				brill::StripInfo strip3 {
					raw.back_strip[k],
					raw.back_energy[k],
					raw.back_time[k],
					parameters.back_p0[raw.back_strip[k]],
					parameters.back_p1[raw.back_strip[k]],
					parameters.back_p2[raw.back_strip[k]],
					parameters.back_p3[raw.back_strip[k]]
				};
				if (
					AreStrips(
						raw.back_strip[i], raw.back_strip[j],
						short_pca.strips[0], short_pca.strips[1],
						short_pca.has_order
					) && (
						abs(raw.back_strip[k] - raw.back_strip[i]) == 1
						|| abs(raw.back_strip[k] - raw.back_strip[j]) == 1
					)
				) {
					combinations.push_back(
						std::make_unique<brill::ShortSharedStripCombination>(
							(0x100 << i) | (0x100 << j) | (0x100 << k),
							raw.back_strip[i] == short_pca.strips[0] ? strip1 : strip2,
							raw.back_strip[i] == short_pca.strips[0] ? strip2 : strip1,
							strip3,
							short_pca.line,
							short_pca.threshold[0]
						)
					);
				} else if (
					AreStrips(
						raw.back_strip[i], raw.back_strip[k],
						short_pca.strips[0], short_pca.strips[1],
						short_pca.has_order
					) && (
						abs(raw.back_strip[j] - raw.back_strip[i]) == 1
						|| abs(raw.back_strip[j] - raw.back_strip[k]) == 1
					)
				) {
					combinations.push_back(
						std::make_unique<brill::ShortSharedStripCombination>(
							(0x100 << i) | (0x100 << j) | (0x100 << k),
							raw.back_strip[i] == short_pca.strips[0] ? strip1 : strip3,
							raw.back_strip[i] == short_pca.strips[0] ? strip3 : strip1,
							strip2,
							short_pca.line,
							short_pca.threshold[0]
						)
					);
				} else if (
					AreStrips(
						raw.back_strip[j], raw.back_strip[k],
						short_pca.strips[0], short_pca.strips[1],
						short_pca.has_order
					) && (
						abs(raw.back_strip[i] - raw.back_strip[j]) == 1
						|| abs(raw.back_strip[i] - raw.back_strip[k]) == 1
					)
				) {
					combinations.push_back(
						std::make_unique<brill::ShortSharedStripCombination>(
							(0x100 << i) | (0x100 << j) | (0x100 << k),
							raw.back_strip[j] == short_pca.strips[0] ? strip2 : strip3,
							raw.back_strip[j] == short_pca.strips[0] ? strip3 : strip2,
							strip1,
							short_pca.line,
							short_pca.threshold[0]
						)
					);
				}
			}
		}
	}
}

void GetT0D2FrontCombinations(
	const SiliconDetectorConfig *detector,
	const DssdEvent &raw,
	const DssdNormalizeParameters &parameters,
	const T0D2ExtraNormalizeParameters &extra,
	const std::map<std::string, FullPCAParameter> &pca_parameters,
	std::vector<std::unique_ptr<StripCombination>> &combinations
) {
	combinations.clear();
	for (int i = 0; i < raw.front_num; ++i) {
		StripInfo strip1 {
			raw.front_strip[i],
			raw.front_energy[i],
			raw.front_time[i],
			parameters.front_p0[raw.front_strip[i]],
			parameters.front_p1[raw.front_strip[i]],
			parameters.front_p2[raw.front_strip[i]],
			parameters.front_p3[raw.front_strip[i]]
		};
		if (raw.front_strip[i] == 75) {
			// piecewise
			combinations.push_back(
				std::make_unique<PiecewiseStripCombination>(
					1 << i,
					strip1,
					extra.pfs75,
					detector->match_tolerance
				)
			);
		} else if (raw.front_strip[i] == 99) {
			// piecewise
			combinations.push_back(
				std::make_unique<brill::PiecewiseStripCombination>(
					1 << i,
					strip1,
					extra.pfs99,
					detector->match_tolerance
				)
			);
		} else {
			// normal
			combinations.push_back(
				std::make_unique<brill::NormalStripCombination>(
					1 << i,
					strip1,
					detector->match_tolerance
				)
			);
		}
		for (int j = i+1; j < raw.front_num; ++j) {
			brill::StripInfo strip2 {
				raw.front_strip[j],
				raw.front_energy[j],
				raw.front_time[j],
				parameters.front_p0[raw.front_strip[j]],
				parameters.front_p1[raw.front_strip[j]],
				parameters.front_p2[raw.front_strip[j]],
				parameters.front_p3[raw.front_strip[j]]
			};
			if (abs(raw.front_strip[i] - raw.front_strip[j]) == 1) {
				if (raw.front_strip[i] == 75 || raw.front_strip[j] == 75) {
					combinations.push_back(
						std::make_unique<brill::PiecewiseSharedStripCombination>(
							(1 << i) | (1 << j),
							raw.front_strip[i] == 75 ? strip1 : strip2,
							raw.front_strip[i] == 75 ? strip2 : strip1,
							extra.pfs75,
							detector->match_tolerance
						)
					);
				} else if (raw.front_strip[i] == 99 || raw.front_strip[j] == 99) {
					combinations.push_back(
						std::make_unique<brill::PiecewiseSharedStripCombination>(
							(1 << i) | (1 << j),
							raw.front_strip[i] == 99 ? strip1 : strip2,
							raw.front_strip[i] == 99 ? strip2 : strip1,
							extra.pfs99,
							detector->match_tolerance
						)
					);
				} else {
					// normal shared
					combinations.push_back(
						std::make_unique<brill::NormalSharedStripCombination>(
							(1 << i) | (1 << j),
							strip1,
							strip2,
							detector->match_tolerance
						)
					);
				}
				continue;
			}
			// broken adjacent
			for (const auto &[name, pca] : pca_parameters) {
				if (pca.side == 1) continue;
				if (!AreStrips(
					raw.front_strip[i], raw.front_strip[j],
					pca.strips[0], pca.strips[1],
					pca.has_order
				)) continue;
				combinations.push_back(
					std::make_unique<brill::BrokenAdjacentStripCombination>(
						(1 << i) | (1 << j),
						raw.front_strip[i] == pca.strips[0] ? strip1 : strip2,
						raw.front_strip[i] == pca.strips[0] ? strip2 : strip1,
						pca.broken_strip,
						pca.line,
						pca.plane[0],
						pca.plane[1],
						pca.threshold[0],
						pca.threshold[1],
						pca.threshold[2]
					)
				);
			}
		}
	}
}

void GetT0D2BackCombinations(
	const SiliconDetectorConfig *detector,
	const DssdEvent &raw,
	const DssdNormalizeParameters &parameters,
	const std::map<std::string, FullPCAParameter> &pca_parameters,
	std::vector<std::unique_ptr<StripCombination>> &combinations
) {
	combinations.clear();
	for (int i = 0; i < raw.back_num; ++i) {
		brill::StripInfo strip1 {
			raw.back_strip[i],
			raw.back_energy[i],
			raw.back_time[i],
			parameters.back_p0[raw.back_strip[i]],
			parameters.back_p1[raw.back_strip[i]],
			parameters.back_p2[raw.back_strip[i]],
			parameters.back_p3[raw.back_strip[i]]
		};
		// normal
		if (raw.back_strip[i] != 104 && raw.back_strip[i] != 109) {
			combinations.push_back(
				std::make_unique<brill::NormalStripCombination>(
					0x100 << i,
					strip1,
					detector->match_tolerance
				)
			);
		}
		for (int j = i+1; j < raw.back_num; ++j) {
			brill::StripInfo strip2 {
				raw.back_strip[j],
				raw.back_energy[j],
				raw.back_time[j],
				parameters.back_p0[raw.back_strip[j]],
				parameters.back_p1[raw.back_strip[j]],
				parameters.back_p2[raw.back_strip[j]],
				parameters.back_p3[raw.back_strip[j]]
			};
			if (abs(raw.back_strip[i] - raw.back_strip[j]) == 1) {
				// normal shared
				combinations.push_back(
					std::make_unique<brill::NormalSharedStripCombination>(
						(0x100 << i) | (0x100 << j),
						strip1,
						strip2,
						detector->match_tolerance
					)
				);
			}
			// broken adjacent
			for (const auto &[name, pca] : pca_parameters) {
				if (pca.side == 0) continue;
				if (!AreStrips(
					raw.back_strip[i], raw.back_strip[j],
					pca.strips[0], pca.strips[1],
					pca.has_order
				)) continue;
				combinations.push_back(
					std::make_unique<brill::BrokenAdjacentStripCombination>(
						(0x100 << i) | (0x100 << j),
						raw.back_strip[i] == pca.strips[0] ? strip1 : strip2,
						raw.back_strip[i] == pca.strips[0] ? strip2 : strip1,
						pca.broken_strip,
						pca.line,
						pca.plane[0],
						pca.plane[1],
						pca.threshold[0],
						pca.threshold[1],
						pca.threshold[2]
					)
				);
			}
		}
	}
}

bool IsInTrackWindow(
	int &d1fs,
	int &d1bs,
	int &d2fs,
	int &d2bs
) {
	int fs_diff = abs(d1fs - d2fs);
	int bs_diff = abs(d1bs - d2bs);
	if (fs_diff <= 1 && bs_diff <= 1) {
		if (d2fs != 18 && d2fs != 19) return true;
		if (d1fs == 18 && d2fs == 19) {
			d2fs = 18;
			return true;
		}
		if (d1fs == 19 && d2fs == 18) {
			d2fs = 19;
			return true;
		}
	}
	if (fs_diff > 1) {
		if (d2fs == 18 && abs(d1fs-19) <= 1) {
			d2fs = 19;
			fs_diff = abs(d1fs-19);
		} else if (d2fs == 19 && abs(d1fs-18) <= 1) {
			d2fs = 18;
			fs_diff = abs(d1fs-18);
		}
	}
	if (bs_diff > 1) {
		if (d1bs == 104 && abs(109-d2bs) <= 1) {
			d1bs = 109;
			bs_diff = abs(109-d2bs);
		} else if (d1bs == 109 && abs(104-d2bs) <= 1) {
			d1bs = 104;
			bs_diff = abs(104-d2bs);
		}
	}
	if (fs_diff <= 1 && bs_diff <= 1) return true;
	return false;
}

void GetPixelPosition(
	const SiliconDetectorConfig *detector,
	const int front_strip,
	const int back_strip,
	double &x,
	double &y,
	double &z
) {
	double tmp_x = detector->center_x_mm
		+ detector->size_x_mm * ((front_strip+0.5)/double(detector->front_strips) - 0.5);
	double tmp_y = detector->center_y_mm
		+ detector->size_y_mm * ((back_strip+0.5)/double(detector->back_strips) - 0.5);
	// rotate 135 degrees in anticlockwise
	x = -0.5*sqrt(2.0)*(tmp_x + tmp_y);
	y = -0.5*sqrt(2.0)*(tmp_x - tmp_y);
	z = detector->z_mm;
}


inline bool CheckStripIndex(const int strip, const int index) {
	if (index == 0 && strip >= 0 && strip <= 27) return true;
	else if (index == 1 && strip >= 26 && strip <= 51) return true;
	else if (index == 2 && strip >= 50 && strip <= 75) return true;
	else if (index == 3 && strip >= 74 && strip <= 99) return true;
	else if (index == 4 && strip >= 98 && strip <= 127) return true;
	return false;
}

bool IsInTrackWindow(const int d2fs, const int d2bs, const int index) {
	int findex = 4 - (index / 5);
	int bindex = index % 5;
	return CheckStripIndex(d2fs, findex) && CheckStripIndex(d2bs, bindex);
}

GAGGCalibrationParameters::GAGGCalibrationParameters(const int counts)
: size(counts) {
	for (int i = 0; i < counts; ++i) {
		p0[i] = 0.0;
		p1[i] = 1.0;
		p2[i] = 0.0;
		p3[i] = 1.0;
	}
}

double GAGGCalibrationParameters::CaliEnergy(const int index, const double raw_energy) const {
	return p0[index] + p1[index]*raw_energy + p2[index]*exp(-raw_energy/p3[index]);
}

int GAGGCalibrationParameters::Write(const std::string &path) const {
	std::filesystem::path file_path(path);
	if (!file_path.parent_path().empty()) {
		std::filesystem::create_directories(file_path.parent_path());
	}
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Write GAGG calibration parameters to " << path << " failed.\n";
		return -1;
	}
	fout << "index p0 p1 p2 p3\n";
	for (int i = 0; i < size; ++i) {
		fout << i << " " << p0[i] << " " << p1[i] << " " << p2[i] << " " << p3[i] << "\n";
	}
	fout.close();
	return 0;
}

int GAGGCalibrationParameters::Read(const std::string &path) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Read GAGG calibration parameters from " << path << " failed.\n";
		return -1;
	}

	std::string line;
	if (!std::getline(fin, line)) {
		std::cerr << "Error: Read GAGG calibration parameters header from " << path << " failed.\n";
		return -2;
	}

	while (std::getline(fin, line)) {
		if (line.empty()) continue;
		std::istringstream iss(line);
		int index = -1;
		double value[4];
		iss >> index;
		for (int i = 0; i < 4; ++i) iss >> value[i];
		if (index < 0 || index >= size) continue;
		p0[index] = value[0];
		p1[index] = value[1];
		p2[index] = value[2];
		p3[index] = value[3];
	}
	fin.close();
	return 0;
}

double GAGGStraightParameters::FixedEnergy(const double d2_energy, const double gagg_energy) const {
	return std::sqrt(d2_energy*gagg_energy + a*d2_energy*d2_energy) + b*gagg_energy;
}

int GAGGStraightParameters::Write(const std::string &path) const {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Errro: Write GAGG straight parameters to " << path << " failed.\n";
		return -1;
	}
	fout << a << " " << b << "\n";
	fout.close();
	return 0;
}

int GAGGStraightParameters::Read(const std::string &path) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Read GAGG straight parameters from " << path << " failed.\n";
		return -1;
	}
	fin >> a >> b;
	fin.close();
	return 0;
}


}
