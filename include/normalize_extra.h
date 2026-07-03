#ifndef __NORMALIZE_EXTRA_H__
#define __NORMALIZE_EXTRA_H__

#include <map>
#include <string>

#include "include/pca.h"

namespace brill {

struct NonLinearParameter {
	double p0[4], p1[4], p2[4];
	double xmin[4], xmax[4], ymin[4], ymax[4];
};

std::istream &operator>>(std::istream &in, NonLinearParameter &non_linear);
std::ostream &operator<<(std::ostream &out, const NonLinearParameter &non_linear);

struct ExtraNormalizeParameters {
	std::map<std::string, PCAParameter> pca;
};

struct T0D1ExtraNormalizeParameters: public ExtraNormalizeParameters {
	NonLinearParameter non_linear;
};

struct T0D2ExtraNormalizeParameters: public ExtraNormalizeParameters {};

int ReadExtraNormalizeParameters(
	const std::string &path,
	T0D1ExtraNormalizeParameters &parameters
);

int WriteExtraNormalizeParameters(
	const std::string &path,
	const T0D1ExtraNormalizeParameters &parameters
);

int ReadExtraNormalizeParameters(
	const std::string &path,
	T0D2ExtraNormalizeParameters &parameters
);

int WriteExtraNormalizeParameters(
	const std::string &path,
	const T0D2ExtraNormalizeParameters &parameters
);


void PCAPrint(
	const std::vector<Eigen::Vector3d> &points,
	const std::string &description,
	const std::string &name,
	const bool assume_line,
	brill::ExtraNormalizeParameters &extra_parameters
);

}
#endif