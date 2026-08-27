#pragma once

#include <Eigen/Dense>

namespace brill {

struct Line3D {
	Eigen::Vector3d centroid;
	Eigen::Vector3d direction;

	inline Eigen::Vector3d ProjectedByPoint(const Eigen::Vector3d &point) const {
		return centroid + direction.dot(point - centroid) * direction;
	}

	inline double DistanceFromPoint(const Eigen::Vector3d &point) const {
		Eigen::Vector3d t = point - centroid;
		return (t - t.dot(direction) * direction).norm();
	}
};

struct Plane {
	Eigen::Vector3d centroid;
	Eigen::Vector3d normal;

	inline Eigen::Vector3d ProjectedByPoint(const Eigen::Vector3d &point) const {
		Eigen::Vector3d t = point - centroid;
		return centroid + (t - t.dot(normal) * normal);
	}

	inline double DistanceFromPoint(const Eigen::Vector3d &point) const {
		return fabs((point-centroid).dot(normal));
	}
};

} // namespace brill
