#include "include/pca.h"

#include <iostream>

namespace brill {

std::ostream &operator<<(std::ostream &os, const PCAParameter &par) {
	os << par.mean(0) << " " << par.mean(1) << " " << par.mean(2) << " "
		<< par.direction(0) << " " << par.direction(1) << " "
		<< par.direction(2) << " " << par.linearity << " " << par.planarity;
	return os;
}

std::istream &operator>>(std::istream &is, PCAParameter &par) {
	is >> par.mean(0) >> par.mean(1) >> par.mean(2)
		>> par.direction(0) >> par.direction(1) >> par.direction(2)
		>> par.linearity >> par.planarity;
	return is;
}

int PrincipalComponentAnalysis(
	const std::vector<Eigen::Vector3d> &points,
	Eigen::Vector3d &eigen_values,
	Eigen::Vector3d &centroid,
	Eigen::Matrix3d &directions
) {
	// PCA of share energy strips
	int n = points.size();
	if (n < 3) return -1;
	Eigen::MatrixXd mat(n, 3);
	// fill matrix
	for (int i = 0; i < n; ++i) {
		mat(i, 0) = points[i](0);
		mat(i, 1) = points[i](1);
		mat(i, 2) = points[i](2);
	}
	// center data
	Eigen::RowVector3d center = mat.colwise().mean();
	Eigen::MatrixXd centered = mat.rowwise() - center;
	centroid= center.transpose();
	// covariance matrix
	Eigen::MatrixXd cov = (centered.adjoint() * centered) / double(n-1);
	// eigen decomposition
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver;
	solver.compute(cov);
	if (solver.info() != Eigen::Success) return -2;
	eigen_values = solver.eigenvalues();
	directions = solver.eigenvectors();
	return 0;
}

double PointVerticalAngle(
	const Eigen::Vector3d &point,
	const Eigen::Vector3d &centroid,
	const Eigen::Vector3d &direction
) {
	Eigen::Vector3d origin(0.0, 0.0, 0.0);
	origin -= centroid;
	Eigen::Vector3d origin_vertical = origin - origin.dot(direction) * direction;
	origin_vertical /= origin_vertical.norm();
	Eigen::Vector3d rpoint = point - centroid;
	Eigen::Vector3d point_vertical = rpoint - rpoint.dot(direction) * direction;
	point_vertical /= point_vertical.norm();
	return point_vertical.dot(origin_vertical);
}

}