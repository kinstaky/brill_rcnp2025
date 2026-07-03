#ifndef __PCA_H__
#define __PCA_H__

#include <iostream>
#include <Eigen/Dense>

namespace brill {

struct PCAParameter {
	double linearity, planarity;
	Eigen::Vector3d mean;
	Eigen::Vector3d direction;
};

std::istream &operator>>(std::istream &in, PCAParameter &parameter);
std::ostream &operator<<(std::ostream &out, const PCAParameter &parameter);

int PrincipalComponentAnalysis(
	const std::vector<Eigen::Vector3d> &points,
	Eigen::Vector3d &eigen_values,
	Eigen::Vector3d &centroid,
	Eigen::Matrix3d &directions
);

inline double PointLineDistance(
	const Eigen::Vector3d &point,
	const Eigen::Vector3d &centroid,
	const Eigen::Vector3d &direction
) {
	Eigen::Vector3d t = point - centroid;
	return (t - t.dot(direction) * direction).norm();
}

inline double PointPlaneDistance(
	const Eigen::Vector3d &point,
	const Eigen::Vector3d &centroid,
	const Eigen::Vector3d &direction
) {
	return fabs((point-centroid).dot(direction));
}

double PointVerticalAngle(
	const Eigen::Vector3d &point,
	const Eigen::Vector3d &centroid,
	const Eigen::Vector3d &direction
);


}
#endif