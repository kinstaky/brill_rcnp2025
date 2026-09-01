#pragma once

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

}