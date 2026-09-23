#pragma once

#include <cmath>
#include "include/config.h"

namespace brill::t0 {

constexpr int kMaxGAGG = 64;

/// @brief Check whether strips in D1D2 are in window
/// @param[inout] d1fs D1 front strip, output corrected result
/// @param[inout] d1bs D1 back strip, output corrected result
/// @param[inout] d2fs D2 front strip, output corrected result
/// @param[inout] d2bs D2 back strip, output corrected result
/// @returns is in window or not
bool IsInTrackWindow(
	int &d1fs,
	int &d1bs,
	int &d2fs,
	int &d2bs
);

/// @brief Get DSSD pixel position
/// @param[in] detector detector config
/// @param[in] front_strip front strip index
/// @param[in] back_strip back strip index
/// @param[out] x x position
/// @param[out] y y position
/// @param[out] z z position
void GetPixelPosition(
	const brill::SiliconDetectorConfig *detector,
	const int front_strip,
	const int back_strip,
	double &x,
	double &y,
	double &z
);

/// @brief Check whether T0D2 and GAGG are in window
/// @param[in] d2fs T0D2 front strip
/// @param[in] d2bs T0D2 back strip
/// @param[in] index GAGG index
/// @returns is in window or not
bool IsInTrackWindow(const int d2fs, const int d2bs, const int index);

class GAGGCalibrationParameters {
public:
	GAGGCalibrationParameters(const int counts);

	double CaliEnergy(const int index, const double raw_energy) const;
	int Write(const std::string &path) const;
	int Read(const std::string &path);

	int size;
	double p0[kMaxGAGG];
	double p1[kMaxGAGG];
	double p2[kMaxGAGG];
	double p3[kMaxGAGG];
};

class GAGGStraightParameters {
public:
	GAGGStraightParameters() = default;

	double FixedEnergy(const double d2_energy, const double gagg_energy) const;
	int Write(const std::string &path) const;
	int Read(const std::string &path);

	double a;
	double b;
};

}