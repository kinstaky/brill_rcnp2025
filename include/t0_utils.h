#pragma once

#include "include/config.h"

namespace brill::t0 {

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

}