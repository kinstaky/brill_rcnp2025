#include "include/normalize_extra.h"

#include <fstream>

namespace brill {

std::istream &operator>>(std::istream &is, NonLinearParameter &non_linear) {
	for (int i = 0; i < 4; ++i) {
		is >> non_linear.xmin[i]
			>> non_linear.xmax[i]
			>> non_linear.ymin[i]
			>> non_linear.ymax[i]
			>> non_linear.p0[i]
			>> non_linear.p1[i]
			>> non_linear.p2[i];
	}
	return is;
}

std::ofstream &operator<<(std::ofstream &os, const NonLinearParameter &non_linear) {
	for (int i = 0; i < 4; ++i) {
		os << non_linear.xmin[i] << " "
			<< non_linear.xmax[i] << " "
			<< non_linear.ymin[i] << " "
			<< non_linear.ymax[i] << " "
			<< non_linear.p0[i] << " "
			<< non_linear.p1[i] << " "
			<< non_linear.p2[i] << "\n";
	}
	return os;
}


int ReadExtraNormalizeParameters(
	const std::string &path,
	T0D1ExtraNormalizeParameters &parameters
) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	fin >> parameters.non_linear;
	while (fin.good()) {
		std::string name;
		PCAParameter par;
		fin >> name >> par;
		parameters.pca.insert(std::make_pair(name, par));
	}
	return 0;
}


int WriteExtraNormalizeParameters(
	const std::string &path,
	const T0D1ExtraNormalizeParameters &parameters
) {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	fout << parameters.non_linear;
	for (const auto &p : parameters.pca) {
		fout << p.first << " " << p.second << "\n";
	}
	return 0;
}


int ReadExtraNormalizeParameters(
	const std::string &path,
	T0D2ExtraNormalizeParameters &parameters
) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	while (fin.good()) {
		std::string name;
		PCAParameter par;
		fin >> name >> par;
		parameters.pca.insert(std::make_pair(name, par));
	}
	return 0;
}


int WriteExtraNormalizeParameters(
	const std::string &path,
	const T0D2ExtraNormalizeParameters &parameters
) {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	for (const auto &p : parameters.pca) {
		fout << p.first << " " << p.second << "\n";
	}
	return 0;
}

void PCAPrint(
	const std::vector<Eigen::Vector3d> &points,
	const std::string &description,
	const std::string &name,
	const bool assume_line,
	brill::ExtraNormalizeParameters &extra_parameters
) {
	brill::PCAParameter par;
	Eigen::Vector3d eigen_values;
	Eigen::Matrix3d eigen_vectors;
	int result = PrincipalComponentAnalysis(
		points,
		eigen_values,
		par.mean,
		eigen_vectors
	);
	if (result == -1) {
		std::cerr << "Error: too few points: " << points.size() << std::endl;
	} else if (result == -2) {
		std::cerr << "Error: Solve covariance matrix failed." << std::endl;
	} else {
		par.linearity = eigen_values(2) / eigen_values.sum();
		par.planarity = eigen_values(0) / eigen_values.sum();
		if (assume_line) par.direction = eigen_vectors.col(2);
		else par.direction = eigen_vectors.col(0);
		std::cout
			<< description << ":\n"
			<< "  Eigen values: " << eigen_values.transpose() << "\n"
			<< "  Linearity: " << par.linearity << "\n"
			<< "  Planarity: " << par.planarity << "\n"
			<< "  Mean: " << par.mean.transpose() << "\n"
			<< "  Direction: " << par.direction.transpose() << "\n";
		extra_parameters.pca.insert(std::make_pair(name, par));
	}
}

}