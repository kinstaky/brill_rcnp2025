#include "include/t0_utils.h"

#include <cmath>
#include <vector>
#include <filesystem>
#include <fstream>
#include <iostream>

namespace brill::t0 {

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
		if (!(iss >> index >> value[0] >> value[1] >> value[2] >> value[3])) continue;
		if (index < 0 || index >= size) continue;
		p0[index] = value[0];
		p1[index] = value[1];
		p2[index] = value[2];
		p3[index] = value[3];
	}
	fin.close();
	return 0;
}

}