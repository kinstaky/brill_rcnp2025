#pragma once

#include <map>
#include <string>

#include "include/pca.h"

namespace brill {


struct Pol2Parameter {
	double p0, p1, p2;
};

struct Range {
	double xmin, xmax, ymin, ymax;
};

std::istream& operator>>(std::istream &is, Pol2Parameter &parameter);
std::ostream& operator<<(std::ostream &os, const Pol2Parameter &parameter);

struct PiecewiseParameter {
	std::vector<Pol2Parameter> pol2;
	std::vector<Range> range;
};

std::istream& operator>>(std::istream &is, PiecewiseParameter &parameter);
std::ostream& operator<<(std::ostream &os, const PiecewiseParameter &parameter);


class ExtraNormalizeParameters {
public:
	virtual ~ExtraNormalizeParameters() = default;

	virtual int Read(const std::string &path) = 0;
	virtual int Write(const std::string &path) const = 0;

	std::map<std::string, PCAParameter> pca;
	std::map<std::string, Pol2Parameter> pol2;
};

class T0D1ExtraNormalizeParameters: public ExtraNormalizeParameters {
public:
	virtual ~T0D1ExtraNormalizeParameters() = default;

	virtual int Read(const std::string &path) override;
	virtual int Write(const std::string &path) const override;

	PiecewiseParameter piecewise;
};

class T0D2ExtraNormalizeParameters: public ExtraNormalizeParameters {
public:
	virtual ~T0D2ExtraNormalizeParameters() = default;

	virtual int Read(const std::string &path) override;
	virtual int Write(const std::string &path) const override;

	PiecewiseParameter pfs75;
	PiecewiseParameter pfs99;
	double rfs17_param[2];
	double rfs20_param[2];
};

void PCAPrint(
	const std::vector<Eigen::Vector3d> &points,
	const std::string &description,
	const std::string &name,
	const bool assume_line,
	brill::ExtraNormalizeParameters &extra_parameters
);

}