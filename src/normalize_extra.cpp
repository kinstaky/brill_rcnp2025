#include "include/normalize_extra.h"

#include <fstream>

namespace brill {

std::istream& operator>>(std::istream &is, Pol2Parameter &parameter) {
	is >> parameter.p0 >> parameter.p1 >> parameter.p2;
	return is;
}

std::ostream& operator<<(std::ostream &os, const Pol2Parameter &parameter) {
	os << parameter.p0 << " " << parameter.p1 << " " << parameter.p2;
	return os;
}

std::istream& operator>>(std::istream &is, PiecewiseParameter &parameter) {
	int count;
	is >> count;

	for (int i = 0; i < count; ++i) {
		Range range;
		Pol2Parameter pol2;
		is >> range.xmin >> range.xmax
			>> range.ymin >> range.ymax
			>> pol2;
		parameter.range.push_back(range);
		parameter.pol2.push_back(pol2);
	}
	return is;
}

std::ostream& operator<<(std::ostream &os, const PiecewiseParameter &parameter) {
	os << parameter.range.size() << "\n";
	for (size_t i = 0; i < parameter.range.size(); ++i) {
		os << parameter.range[i].xmin << " "
			<< parameter.range[i].xmax << " "
			<< parameter.range[i].ymin << " "
			<< parameter.range[i].ymax << " "
			<< parameter.pol2[i] << "\n";
	}
	return os;
}


int T0D1ExtraNormalizeParameters::Read(const std::string &path) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	char c;
	std::string type, name;
	int count;
	fin >> c >> type >> count;
	while (c == '#' && fin.good()) {
		if (type == "piecewise") {
			fin >> piecewise;
		} else if (type == "pol2") {
			for (int i = 0; i < count; ++i) {
				Pol2Parameter par;
				fin >> name >> par;
				pol2.insert(std::make_pair(name, par));
			}
		} else if (type == "pca") {
			for (int i = 0; i < count; ++i) {
				PCAParameter par;
				fin >> name >> par;
				pca.insert(std::make_pair(name, par));
			}
		}
		fin >> c >> type >> count;
	}
	return 0;
}


int T0D1ExtraNormalizeParameters::Write(const std::string &path) const {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}

	fout << "# piecewise 1\n";
	fout << piecewise;

	fout << "# pca " << pca.size() << "\n";
	for (const auto &p : pca) {
		fout << p.first << " " << p.second << "\n";
	}

	fout << "# pol2 " << pol2.size() << "\n";
	for (const auto &p : pol2) {
		fout << p.first << " " << p.second << "\n";
	}
	return 0;
}


int T0D2ExtraNormalizeParameters::Read(const std::string &path) {
	std::ifstream fin(path);
	if (!fin.good()) {
		std::cerr << "Error: Open normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}
	char c;
	std::string type, name;
	int count;
	fin >> c >> type >> count;
	while (c == '#' && fin.good()) {
		if (type == "pol2") {
			for (int i = 0; i < count; ++i) {
				Pol2Parameter par;
				fin >> name >> par;
				pol2.insert(std::make_pair(name, par));
			}
		} else if (type == "pca") {
			for (int i = 0; i < count; ++i) {
				PCAParameter par;
				fin >> name >> par;
				pca.insert(std::make_pair(name, par));
			}
		} else if (type == "rs") {
			fin >> rfs17_param[0] >> rfs17_param[1]
				>> rfs20_param[0] >> rfs20_param[1];
		}
		fin >> c >> type >> count;
	}
	return 0;
}


int T0D2ExtraNormalizeParameters::Write(const std::string &path) const {
	std::ofstream fout(path);
	if (!fout.good()) {
		std::cerr << "Error: Open output normalize parameter file "
			<< path << " failed.\n";
		return -1;
	}

	fout << "# rs 1 "
		<< rfs17_param[0] << " " << rfs17_param[1] << " "
		<< rfs20_param[0] << " " << rfs20_param[1] << "\n";

	fout << "# pca " << pca.size() << "\n";
	for (const auto &p : pca) {
		fout << p.first << " " << p.second << "\n";
	}

	fout << "# pol2 " << pol2.size() << "\n";
	for (const auto &p : pol2) {
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