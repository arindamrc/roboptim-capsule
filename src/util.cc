// Copyright (C) 2012 by Antonio El Khoury, CNRS.
//
// This file is part of the roboptim-capsule.
//
// roboptim-capsule is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// roboptim-capsule is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with roboptim-capsule.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ROBOPTIM_CAPSULE_UTIL_CC_
# define ROBOPTIM_CAPSULE_UTIL_CC_

# include <iostream>
# include <set>
# include <limits>
# include <cstdio>

# include <boost/foreach.hpp>

# include <roboptim/capsule/util.hh>
# include "pugixml.cpp"

namespace roboptim {
namespace capsule {

polyhedron_t convexHullFromPoints(const std::vector<point_t>& points) {
	polyhedron_t convexPolyhedron;

# ifdef HAVE_QHULL
	int numpoints = static_cast<int> (points.size ());
	int dim = 3;

	std::vector<value_type> rboxpoints (dim * numpoints);
	size_t iter = 0;

	for (std::vector<point_t>::const_iterator
			point = points.begin (); point != points.end (); ++point)
	{
		rboxpoints[iter++] = (*point)[0];
		rboxpoints[iter++] = (*point)[1];
		rboxpoints[iter++] = (*point)[2];
	}

	// Compute the convex hull with qhull
	char flags[25];
	sprintf (flags, "qhull Qc Qt Qi");
	// Note: using stderr instead of NULL to avoid a bug in older versions of
	// qhull:
	// > QH6232 Qhull internal error (userprintf.c): fp is 0.  Wrong
	// > qh_fprintf called.
	int exitcode = qh_new_qhull (dim, numpoints,
			rboxpoints.data (), 0,
			flags, NULL, stderr);

	if (exitcode != 0)
	{
		// Return empty polyhedron
		return convexPolyhedron;
	}

	// Get the list of points
	convexPolyhedron.clear ();
	vertexT* vertex = qh vertex_list;

	// TODO: find how to get the number of vertices on the convex
	// hull with qhull.

	FORALLvertices {
		// qh_pointid (vertex->point) is the point id of the vertex
		//int id = qh_pointid (vertex->point);
		// vertex->point is the coordinates of the vertex
		convexPolyhedron.push_back(point_t (vertex->point[0],
						vertex->point[1],
						vertex->point[2]));
	}

	qh_freeqhull (!qh_ALL);
# else
	std::cerr << "Qhull not found, cannot compute the convex hull."
			<< std::endl;
# endif //! HAVE_QHULL
	// Return the convex hull as a polyhedron
	return convexPolyhedron;
}

value_type distancePointToSegment(const point_t& p, const point_t& a,
		const point_t& b) {
	value_type d_ab = (b - a).norm();

	// If the segment is a point, i.e. a = b
	if (d_ab < 1e-6)
		return (a - p).norm();

	return (p - projectionOnSegment(p, a, b)).norm();
}

point_t projectionOnSegment(const point_t& p, const point_t& a,
		const point_t& b) {
	value_type d_ab = (b - a).norm();

	// If the segment is a point, i.e. a = b
	if (d_ab < 1e-6)
		return a;

	// We note q the projection of p on the line (a,b)
	value_type d_aq = (p - a).dot(b - a) / d_ab;
	if (d_aq > d_ab)
		return b;
	else if (d_aq < 0.)
		return a;
	else
		return a + d_aq * (b - a).normalized();
}

value_type distancePointToLine(const point_t& point, const point_t& linePoint,
		const vector3_t& dir) {
	assert(dir.norm() > 0.);
	return (dir.cross(linePoint - point)).norm() / dir.norm();
}

Eigen::Matrix3d covarianceMatrix(const std::vector<point_t>& points) {
	value_type oon = 1.0 / (value_type) points.size();
	point_t c(0., 0., 0.);
	value_type e00, e11, e22, e01, e02, e12;

	// compute the center of mass of the points
	for (size_t i = 0; i < points.size(); ++i)
		c += points[i];
	c *= oon;

	// compute covariance elements
	e00 = e11 = e22 = e01 = e02 = e12 = 0.0;
	for (size_t i = 0; i < points.size(); ++i) {
		// translate points so center of mass is at origin
		point_t p = points[i] - c;
		// compute covariance of translated points
		e00 += p[0] * p[0];
		e11 += p[1] * p[1];
		e22 += p[2] * p[2];
		e01 += p[0] * p[1];
		e02 += p[0] * p[2];
		e12 += p[1] * p[2];
	}

	// fill in the covariance matrix elements
	Eigen::Matrix3d cov;

	cov(0, 0) = e00 * oon;
	cov(1, 1) = e11 * oon;
	cov(2, 2) = e22 * oon;
	cov(0, 1) = cov(1, 0) = e01 * oon;
	cov(0, 2) = cov(2, 0) = e02 * oon;
	cov(1, 2) = cov(2, 1) = e12 * oon;

	return cov;
}

// Returns indices imin and imax into pt[] array of the least and
// most, respectively, distant points along the direction dir
void extremePointsAlongDirection(vector3_t dir,
		const std::vector<point_t>& points, int& imin, int& imax) {
	double minproj = std::numeric_limits<double>::max();
	double maxproj = -minproj;

	for (size_t i = 0; i < points.size(); ++i) {
		// Project vector from origin to point onto direction vector
		double proj = points[i].dot(dir);
		// Keep track of least distant point along direction vector
		if (proj < minproj) {
			minproj = proj;
			imin = static_cast<int>(i);
		}
		// Keep track of most distant point along direction vector
		if (proj > maxproj) {
			maxproj = proj;
			imax = static_cast<int>(i);
		}
	}
}

Capsule capsuleFromPoints(const std::vector<point_t>& points) {
	assert(points.size() > 0 && "Cannot compute capsule for empty polyhedron.");

	// Create the covariance matrix for PCA
	Eigen::Matrix3d covariance = covarianceMatrix(points);

	// Compute eigenvectors and eigenvalues
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
	es.compute(covariance, Eigen::ComputeEigenvectors);

	Eigen::Matrix3d eigenVectors = es.eigenvectors();
	Eigen::Matrix3d eigenValues = es.eigenvalues().asDiagonal();

	// Find the largest eigenvalue and the corresponding direction
	// (largest spread)
	unsigned maxc, minc;
	maxc = minc = 0;
	value_type absev, maxev, minev;
	maxev = minev = std::fabs(eigenValues(0, 0));

	if ((absev = std::fabs(eigenValues(1, 1))) > maxev) {
		maxc = 1;
		maxev = absev;
	} else {
		minc = 1;
		minev = absev;
	}

	if ((absev = std::fabs(eigenValues(2, 2))) > maxev) {
		maxc = 2;
		maxev = absev;
	} else if (minev > absev) {
		minc = 2;
		minev = absev;
	}

	vector3_t dirLargestSpread = eigenVectors.col(maxc);
	dirLargestSpread.normalize();

	// Find the most extreme points along the largest spread direction.
	// Those points will help to find the length of the capsule.
	int iminLargestSpread = 0;
	int imaxLargestSpread = 0;
	extremePointsAlongDirection(dirLargestSpread, points, iminLargestSpread,
			imaxLargestSpread);
	point_t minptLargestSpread = points[iminLargestSpread];
	point_t maxptLargestSpread = points[imaxLargestSpread];

	// Compute the start point
	// The cylinder axis will be (average point, largest spread direction).
	// However, a better point could be found with a more complicated
	// algorithm, thus reducing the volume of the capsule.
	point_t average(0., 0., 0.);
	for (size_t i = 0; i < points.size(); ++i) {
		average += points[i];
	}
	average /= static_cast<value_type>(points.size());

	// Find the correct radius for the capsule.
	value_type radius = 0;
	for (size_t i = 0; i < points.size(); ++i) {
		value_type dist = distancePointToLine(points[i], average,
				dirLargestSpread);
		if (dist > radius)
			radius = dist;
	}

	// Find the correct length for the capsule (cylinder part)
	value_type length = (maxptLargestSpread - minptLargestSpread).norm();

	// Length used to find the correct center position on the direction axis
	value_type maxLengthFromAverage = std::fabs(
			(maxptLargestSpread - average).dot(dirLargestSpread));
	point_t center = average
			+ (maxLengthFromAverage - 0.5 * length) * dirLargestSpread;

	// Optimization of the volume
	// - We determine the points located at
	//   both extremities (+/-)(0.5 * length - radius)
	// - For all of those points, we look for the start/endpoint position
	//   that will minimize the capsule volume.
	std::vector<point_t> nearStartPoints;
	std::vector<point_t> nearEndPoints;
	point_t start = center - (0.5 * length - radius) * dirLargestSpread;
	point_t end = center + (0.5 * length - radius) * dirLargestSpread;

	for (size_t i = 0; i < points.size(); ++i) {
		// if located near the start boundary
		value_type dirDist = dirLargestSpread.dot(points[i] - center);
		if (-dirDist > 0.5 * length - radius)
			nearStartPoints.push_back(points[i]);
		// else if located near the end boundary
		else if (dirDist > 0.5 * length - radius)
			nearEndPoints.push_back(points[i]);
	}

	// we move the position of the start point to include all points in its
	// vicinity
	for (size_t i = 0; i < nearStartPoints.size(); ++i) {
		if ((nearStartPoints[i] - start).norm() > radius) {
			// using pythagore theorem
			value_type h = distancePointToLine(nearStartPoints[i], center,
					dirLargestSpread);
			value_type l = (nearStartPoints[i] - start).dot(-dirLargestSpread);
			if (l - sqrt(radius * radius - h * h) > 0)
				start -= (l - sqrt(radius * radius - h * h)) * dirLargestSpread;
		}
	}

	// we move the position of the end point to include all points in its
	// vicinity
	for (size_t i = 0; i < nearEndPoints.size(); ++i) {
		if ((nearEndPoints[i] - end).norm() > radius) {
			// using pythagore theorem
			value_type l = (nearEndPoints[i] - end).dot(dirLargestSpread);
			value_type h = distancePointToLine(nearEndPoints[i], center,
					dirLargestSpread);
			if (l - std::sqrt(radius * radius - h * h) > 0)
				end += (l - std::sqrt(radius * radius - h * h))
						* dirLargestSpread;
		}
	}

	Capsule capsule;
	capsule.P0 = start;
	capsule.P1 = end;
	capsule.radius = radius;

	return capsule;
}

void convertCapsuleToSolverParam(argument_ref dst, const point_t& endPoint1,
		const point_t& endPoint2, const value_type& radius) {
	dst.resize(7);

	dst[0] = endPoint1[0];
	dst[1] = endPoint1[1];
	dst[2] = endPoint1[2];
	dst[3] = endPoint2[0];
	dst[4] = endPoint2[1];
	dst[5] = endPoint2[2];
	dst[6] = radius;
}

void convertSolverParamToCapsule(point_t& endPoint1, point_t& endPoint2,
		value_type& radius, const argument_t src) {
	assert(src.size() == 7 && "Incorrect src size, expected 7.");
	assert(
			src[6] > 0
					&& "Invalid value for radius, expected non-negative value.");

	endPoint1[0] = src[0];
	endPoint1[1] = src[1];
	endPoint1[2] = src[2];
	endPoint2[0] = src[3];
	endPoint2[1] = src[4];
	endPoint2[2] = src[5];
	radius = src[6];
}

void convertPolyhedronVectorToPolyhedron(polyhedron_t& polyhedron,
		const polyhedrons_t& polyhedrons) {
	assert(polyhedrons.size() != 0 && "Empty polyhedron vector.");
	assert(polyhedron.size() == 0 && "Union polyhedron must be empty.");

	BOOST_FOREACH (polyhedron_t poly, polyhedrons) {
		BOOST_FOREACH (point_t point, poly)
			polyhedron.push_back(point);
	}
}

void computeBoundingCapsulePolyhedron(const polyhedrons_t& polyhedrons,
		point_t& endPoint1, point_t& endPoint2, value_type& radius) {
	assert(polyhedrons.size() != 0 && "Empty polyhedron vector.");

	// Retrieve vector of points from polyhedrons.
	size_type nbPoints = 0;
	BOOST_FOREACH (polyhedron_t polyhedron, polyhedrons) {
		nbPoints += polyhedron.size();
	}

	std::vector<point_t> points(nbPoints);
	size_type k = 0;
	BOOST_FOREACH (polyhedron_t polyhedron, polyhedrons) {
		BOOST_FOREACH (point_t point, polyhedron) {
			points[k++] = point;
		}
	}

	// Compute bounding capsule of points.
	Capsule capsule = capsuleFromPoints(points);

	// Get capsule parameters.
	endPoint1 = capsule.P0;
	endPoint2 = capsule.P1;
	radius = capsule.radius;
}

void computeConvexPolyhedron(const polyhedrons_t& polyhedrons,
		polyhedrons_t& convexPolyhedrons) {
	assert(polyhedrons.size() != 0 && "Empty polyhedron vector.");
	assert(
			convexPolyhedrons.size() == 0
					&& "Convex polyhedron vector must be empty.");

	polyhedron_t polyhedron;
	convertPolyhedronVectorToPolyhedron(polyhedron, polyhedrons);

	assert(polyhedron.size() > 0 && "Polyhedron merging failed.");

	// Build convex polyhedron that contains unique points.
	polyhedron_t convexPolyhedron = convexHullFromPoints(polyhedron);

	assert(
			convexPolyhedron.size() > 0
					&& "Convex polyhedron computation failed.");

	convexPolyhedrons.push_back(convexPolyhedron);
}

} // end of namespace capsule.

namespace geometry {
Eigen::Matrix4d create_transformation_matrix(Eigen::Vector3d rpy,
		Eigen::Vector3d xyz, Eigen::Vector3d scale) {
	double ax, ay, az, tx, ty, tz, sx, sy, sz;
	ax = rpy(0);
	ay = rpy(1);
	az = rpy(2);
	tx = xyz(0);
	ty = xyz(1);
	tz = xyz(2);
	sx = scale(0);
	sy = scale(1);
	sz = scale(2);
	Eigen::Affine3d rx = Eigen::Affine3d(
			Eigen::AngleAxisd(ax, Eigen::Vector3d::UnitX()));
	Eigen::Affine3d ry = Eigen::Affine3d(
			Eigen::AngleAxisd(ay, Eigen::Vector3d::UnitY()));
	Eigen::Affine3d rz = Eigen::Affine3d(
			Eigen::AngleAxisd(az, Eigen::Vector3d::UnitZ()));
	Eigen::Affine3d r = rz * ry * rx;
	Eigen::Affine3d t(Eigen::Translation3d(Eigen::Vector3d(tx, ty, tz)));
	Eigen::Affine3d s(Eigen::Scaling(Eigen::Vector3d(sx, sy, sz)));
	Eigen::Matrix4d m = (t * r * s).matrix();
	return m;
}

Eigen::Vector3d axisAngles(Eigen::Vector3d rot_i) {
	Eigen::Vector3d z(Eigen::Vector3d::UnitZ());

	Eigen::Vector3d a_n = rot_i.normalized();

	Eigen::Vector3d N = z.cross(a_n); // 	 default alignment of the cylinder is with the z-axis

	double sine = N.norm();
	double cosine = z.dot(a_n);

	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d Nx = Eigen::Matrix3d::Zero();

	Eigen::Vector3d x = N.normalized();

	Nx << 0, -x(2), x(1), x(2), 0, -x(0), -x(1), x(0), 0;

//	Eigen::Matrix3d R = I + Nx + ((1.0 / (1 + cosine)) * (Nx * Nx));
	Eigen::Matrix3d R = I + sine * Nx + (1.0 - cosine) * (Nx * Nx);

	double yaw = atan2((double) R(1, 0), (double) R(0, 0));
	double pitch = atan2((double) -R(2, 0),
			(double) pow((double) (R(2, 1) * R(2, 1) + R(2, 2) * R(2, 2)),
					0.5));
	double roll = atan2((double) R(2, 1), (double) R(2, 2));

	Eigen::Vector3d rpy = Eigen::Vector3d(roll, pitch, yaw);
	return rpy;
}

/* Get a vector from a space separated value list*/
void getAsVector(const char* str, std::vector<double>* vec) {
	std::vector<double> retVal;
	std::istringstream vertexStream(str);
	double val;
	while (vertexStream >> val) {
		vec->push_back(val);
	}
}

} // end of namespace geometry.

} // end of namespace roboptim.

#endif //! ROBOPTIM_CAPSULE_UTIL_CC_
