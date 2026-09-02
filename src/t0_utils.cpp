#include "include/t0_utils.h"

#include <cmath>
#include <vector>

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

}