// Copyright (C) 2014 by Benjamin Chrétien, CNRS-LIRMM.
//
// This file is part of the roboptim-capsule.
//
// roboptim-capsule is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// roboptim-capsule is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with roboptim-capsule.  If not, see
// <http://www.gnu.org/licenses/>.

/**
 * \file src/capsule-generator.cc
 *
 * \brief CLI capsule generator.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <regex>
#include <math.h>

#include <boost/program_options.hpp>
#include <boost/regex.hpp>

#include <roboptim/capsule/fitter.hh>
#include <roboptim/capsule/util.hh>

#include <roboptim/capsule/pugiconfig.hpp>
#include <roboptim/capsule/pugixml.hpp>

#include <eigen3/Eigen/Geometry>

#include "pugixml.cpp"

using namespace roboptim;
using namespace roboptim::capsule;

const double side = 1;

Eigen::Vector3d axisAngles(Eigen::Vector3d a) {
	Eigen::Vector3d z(Eigen::Vector3d::UnitZ());

	Eigen::Vector3d a_n = a.normalized();

	Eigen::Vector3d N = a_n.cross(z); // default alignment of the cylinder is with the z-axis
	std::cout << "a: " << a << ", z: " << z << ", N: " << N << std::endl;

	double sine = N.norm();
	double cosine = a_n.dot(z);

	std::cout << "sin: " << sine << ", cos: " << cosine << std::endl;

	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	Eigen::Matrix3d Nx = Eigen::Matrix3d::Zero();

	Nx << 0, -N(2), N(1), N(2), 0, -N(0), -N(1), N(0), 0;

	Eigen::Matrix3d R = I + Nx + ((1.0 / (1 + cosine)) * (Nx * Nx));

	std::cout << "Nx: " << Nx << std::endl;
	std::cout << "Nx2: " << (Nx * Nx) << std::endl;
	std::cout << "scalar: " << (1.0 / (1 + cosine)) << std::endl;
	std::cout << "R: " << R << std::endl;

	double yaw = atan2((double) R(1, 0), (double) R(0, 0));
	double pitch = atan2((double) -R(2, 0),
			(double) pow((double) (R(2, 1) * R(2, 1) + R(2, 2) * R(2, 2)),
					0.5));
	double roll = atan2((double) R(2, 1), (double) R(2, 2));

	Eigen::Vector3d rpy = Eigen::Vector3d(roll, pitch, yaw);
	std::cout << "rpy: " << rpy << std::endl;
	return rpy;
}

Eigen::Matrix4d create_transformation_matrix(double ax, double ay, double az,
		double tx, double ty, double tz, double sx, double sy, double sz) {
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

std::vector<double> getPointsFromDAE(const char* filename) {
	std::vector<double> points;
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(filename);
	std::cout << "Load meshfile: " << result.description()
			<< ", meshfile child name: " << doc.first_child().name()
			<< std::endl;

	const std::string path =
			"/COLLADA/library_geometries/geometry/mesh/source[*]/float_array[contains(@id,'positions')]";
	const pugi::xpath_node_set nodes = doc.select_nodes(path.c_str());

	for (const pugi::xpath_node* it = nodes.begin(); it != nodes.end(); ++it) {
		const char* vertexString = it->node().text().get();
		std::istringstream vertexStream(vertexString);
		double point;
		while (vertexStream >> point) {
			points.push_back(point);
		}
	}
	return points;
}

Eigen::Matrix4d getTransformationMatriyFromDAE(const char* filename) {
	Eigen::Matrix4d t;
	std::vector<double> els;
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(filename);
	std::cout << "Load meshfile: " << result.description() << std::endl;
	const std::string path =
			"/COLLADA/library_visual_scenes/visual_scene[@id=\"Scene\"]/node/matrix[@sid=\"transform\"]";
	const pugi::xpath_node_set nodes = doc.select_nodes(path.c_str());

	for (const pugi::xpath_node* it = nodes.begin(); it != nodes.end(); ++it) {
		const char* elString = it->node().text().get();
		std::cout << "matrix string: " << elString << std::endl;
		std::istringstream elStream(elString);
		double el;
		while (elStream >> el) {
			els.push_back(el);
		}
	}
	t = Eigen::Matrix4d(els.data());
	std::cout << "collada transform matrix: " << t << std::endl;
	return t;
}

Fitter getCapsuleParams(std::vector<double> points) {
	std::string solver = "ipopt";

	if (points.size() % 3 != 0) {
		throw std::invalid_argument(
				"Error: points should be an array of 3D points, e.g. x0 y0 z0 x1 y1 z1 etc.");
	}

	// Load polyhedron
	polyhedron_t polyhedron;
	for (size_t i = 0; i < points.size(); i += 3) {
		point_t p(points[i], points[i + 1], points[i + 2]);
		polyhedron.push_back(p);
	}

	std::cout << "polyhedron constructed" << std::endl;

	// Fitter expects a vector of polyhedrons
	polyhedrons_t polyhedrons;
	polyhedrons.push_back(polyhedron);

	// Create fitter
	Fitter fitter(polyhedrons, solver);

	fitter.logDirectory() = "log";
	// Compute initial guess
	point_t P0;
	point_t P1;
	value_type r = 0.;
	argument_t initParam(7);

	polyhedrons_t convexPolyhedrons;
	computeConvexPolyhedron(polyhedrons, convexPolyhedrons);
	std::cout << "convex polyhedron computed" << std::endl;
	computeBoundingCapsulePolyhedron(convexPolyhedrons, P0, P1, r);
	std::cout << "initial bounding capsule computed" << std::endl;
	convertCapsuleToSolverParam(initParam, P0, P1, r);

	// Compute optimal capsule
	fitter.computeBestFitCapsule(initParam);
	std::cout << "best capsule computed" << std::endl;

	return fitter;
}

std::vector<double> parseComponents(const char* str) {
	std::cout << "parsing components from: " << str << std::endl;
	std::vector<double> components;
	std::istringstream stream(str);
	std::string ci;

	boost::smatch mXacro, mValue;
	boost::regex eXacro("\\${side\\*\\(([-+]?[0-9]*\.[0-9]+|[0-9]+)\\)}");
	boost::regex eValue("[-+]?([0-9]*\.[0-9]+|[0-9]+)");

	while (stream >> ci) {
		std::cout << "ci: " << ci << std::endl;
		if (boost::regex_search(ci, mXacro, eXacro)) {
//			std::cout << "match: " << mXacro[0].str() << std::endl;
//			std::cout << "g1: " << mXacro[1].str() << std::endl;
			components.push_back(side * atof(mXacro[1].str().c_str()));
		} else {
			components.push_back(atof(ci.c_str()));
		}
	}

//	std::cout << components << std::endl;
	return components;
}

void makeCapsuleFromXml(pugi::xml_node linkNode, std::vector<double> e1,
		std::vector<double> e2, double r, const char* rpyStr,
		const char* xyzStr, const char* scaleStr, Eigen::Matrix4d tDAE) {

	printf(
			"raw: e1[0]: %lf, e1[1]: %lf, e1[2]: %lf, e2[0]: %lf, e2[1]: %lf, e2[2]: %lf, r: %lf\n",
			e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], r);

	std::vector<double> rpy(3);
	std::vector<double> xyz(3);
	std::vector<double> scale(3);
//	double scale = 0.012;

//	endpoint1 *= scale;
//	endpoint2 *= scale;
//	r *= scale;

	if (strlen(scaleStr) > 0) {
		scale = parseComponents(scaleStr);
	} else {
		scale[0] = 1;
		scale[1] = 1;
		scale[2] = 1;
	}

	if (strlen(rpyStr) > 0) {
		rpy = parseComponents(rpyStr);
	} else {
		rpy[0] = 0;
		rpy[1] = 0;
		rpy[2] = 0;
	}

	if (strlen(xyzStr) > 0) {
		xyz = parseComponents(xyzStr);
	} else {
		xyz[0] = 0;
		xyz[1] = 0;
		xyz[2] = 0;
	}

	std::cout << "xyz: " << xyz << std::endl;
	std::cout << "rpy: " << rpy << std::endl;
	std::cout << "scale: " << scale << std::endl;

	Eigen::Matrix4d transform = create_transformation_matrix(rpy[0], rpy[1],
			rpy[2], xyz[0], xyz[1], xyz[2], scale[0], scale[1], scale[2]);

	std::cout << "transformation matrix: " << transform << std::endl;

//	e1[0] *= scale[0];
//	e1[1] *= scale[1];
//	e1[2] *= scale[2];
//
//	e2[0] *= scale[0];
//	e2[1] *= scale[1];
//	e2[2] *= scale[2];

	Eigen::Vector3d endpoint1(e1[0], e1[1], e1[2]);
	Eigen::Vector3d endpoint2(e2[0], e2[1], e2[2]);

	Eigen::Vector4d homoe1 = endpoint1.homogeneous();
	Eigen::Vector4d homoe2 = endpoint2.homogeneous();
	std::cout << "h1: " << homoe1 << " | " << "h2: " << homoe2 << std::endl;

	homoe1 = transform * tDAE * homoe1;
	homoe2 = transform * tDAE * homoe2;
	r = r * tDAE(0, 0);

	std::cout << "h1: " << homoe1 << " | " << "h2: " << homoe2 << std::endl;
	std::cout << "e1: " << endpoint1 << " | " << "e2: " << endpoint2
			<< std::endl;

	endpoint1 = homoe1.hnormalized();
	endpoint2 = homoe2.hnormalized();

	Eigen::Vector3d capAxis = endpoint2 - endpoint1;
	std::cout << "e1: " << endpoint1 << " | " << "e2: " << endpoint2
			<< std::endl;

	double caplen = (endpoint1 - endpoint2).norm();
	Eigen::Vector3d midpoint = (endpoint1 + endpoint2) * 0.5;
	std::cout << "capsule len: " << caplen << " | " << "midpoint: " << midpoint
			<< std::endl;

	Eigen::Vector3d rot = axisAngles(capAxis);
	std::cout << "new rotation: " << rot << std::endl;

//	std::cout << "scale: " << scale << std::endl;
//	std::cout << "cyl len mag: " << cyl_l << std::endl;
//	std::cout << "cyl len: " << len << std::endl;

//	printf(
//			"scaled: e1[0]: %lf, e1[1]: %lf, e1[2]: %lf, e2[0]: %lf, e2[1]: %lf, e2[2]: %lf, r: %lf\n",
//			e1[0], e1[1], e1[2], e2[0], e2[1], e2[2], r);

	std::cout << "constructing s1" << std::endl;
	// The first sphere
	pugi::xml_node s1Collision = linkNode.append_child("collision");
	pugi::xml_node s1Origin = s1Collision.append_child("origin");
	std::stringstream s1rpySS;
	s1rpySS << rot(0) << " " << rot(1) << " " << rot(2);
	std::cout << "s1rpySS: " << s1rpySS.str() << std::endl;
	s1Origin.append_attribute("rpy").set_value(s1rpySS.str().c_str());
	std::stringstream s1xyzSS;
	s1xyzSS << (endpoint1(0)) << " " << (endpoint1(1)) << " " << (endpoint1(2));
	std::cout << "s1xyzSS: " << s1xyzSS.str() << std::endl;
	s1Origin.append_attribute("xyz").set_value(s1xyzSS.str().c_str());
	pugi::xml_node s1Geom = s1Collision.append_child("geometry");
	pugi::xml_node s1 = s1Geom.append_child("sphere");
	pugi::xml_attribute rs1 = s1.append_attribute("radius");
	rs1.set_value(r);

	std::cout << "constructing cyl" << std::endl;
	// The middle cylinder
	pugi::xml_node crCollision = linkNode.append_child("collision");
	pugi::xml_node crOrigin = crCollision.append_child("origin");
	std::stringstream crrpySS;
	crrpySS << rot(0) << " " << rot(1) << " " << rot(2);
	crOrigin.append_attribute("rpy").set_value(crrpySS.str().c_str());
	std::stringstream crxyzSS;
	crxyzSS << (midpoint(0)) << " " << (midpoint(1)) << " " << (midpoint(2));
	std::cout << "crxyzSS: " << crxyzSS.str() << std::endl;
	crOrigin.append_attribute("xyz").set_value(crxyzSS.str().c_str());
	pugi::xml_node crGeom = crCollision.append_child("geometry");
	pugi::xml_node cyl = crGeom.append_child("cylinder");
	pugi::xml_attribute crl = cyl.append_attribute("length");
	crl.set_value(caplen);
	pugi::xml_attribute crr = cyl.append_attribute("radius");
	crr.set_value(r);

	std::cout << "constructing s2" << std::endl;
	// The second sphere
	pugi::xml_node s2Collision = linkNode.append_child("collision");
	pugi::xml_node s2Origin = s2Collision.append_child("origin");
	std::stringstream s2rpySS;
	s2rpySS << rot(0) << " " << rot(1) << " " << rot(2);
	s2Origin.append_attribute("rpy").set_value(s1rpySS.str().c_str());
	std::stringstream s2xyzSS;
	s2xyzSS << (endpoint2(0)) << " " << (endpoint2(1)) << " " << (endpoint2(2));
	s2Origin.append_attribute("xyz").set_value(s2xyzSS.str().c_str());
	pugi::xml_node s2Geom = s2Collision.append_child("geometry");
	pugi::xml_node s2 = s2Geom.append_child("sphere");
	pugi::xml_attribute s2r = s2.append_attribute("radius");
	s2r.set_value(r);

}

int main_old(int argc, char **argv) {
	const char* filename = "./meshes/untitled.stl";
	std::vector<double> points; //= getPoints(filename);
	std::ifstream infile(filename);
	std::string line;
	std::getline(infile, line);
	std::istringstream vertexStream(line);
	double point;
	while (vertexStream >> point) {
		points.push_back(point);
	}
	std::cout << points.size() << " points pushed!" << std::endl;

	Fitter fitter = getCapsuleParams(points);
	std::vector<double> e1(3), e2(3);
	double radius;
	std::cout << "Initial: " << fitter.initParam() << std::endl;
	std::cout << "Solution: " << fitter.solutionParam() << std::endl;

	// assign values to e1,e2 and radius
	e1[0] = fitter.solutionParam()(0, 0);
	e1[1] = fitter.solutionParam()(1, 0);
	e1[2] = fitter.solutionParam()(2, 0);

	std::cout << "e1: 0,1,2: " << e1[0] << "," << e1[1] << "," << e1[2]
			<< std::endl;

	e2[0] = fitter.solutionParam()(3, 0);
	e2[1] = fitter.solutionParam()(4, 0);
	e2[2] = fitter.solutionParam()(5, 0);

	std::cout << "e2: 0,1,2: " << e2[0] << "," << e2[1] << "," << e2[2]
			<< std::endl;

	radius = fitter.solutionParam()(6, 0);

	std::cout << "rad: " << radius << std::endl;
	pugi::xml_node tmp;

//	makeCapsuleFromXml(tmp, e1, e2, radius, "", "", "");
	return 1;
}

int main(int argc, char** argv) {

	try {
		pugi::xml_document doc;
//		pugi::xml_parse_result result = doc.load_file("./src/svh.urdf.xacro");
		pugi::xml_parse_result result = doc.load_file("./src/svh.urdf");
//		pugi::xml_parse_result result = doc.load_file("./src/urdf/test.urdf");
		std::cout << "Load result: " << result.description()
				<< ", 1st child name: " << doc.first_child().name()
				<< std::endl;

		const std::string path = "/robot/link[*]/collision/geometry/mesh";
		const pugi::xpath_node_set nodes = doc.select_nodes(path.c_str());

		for (const pugi::xpath_node* it = nodes.begin(); it != nodes.end();
				++it) {
			pugi::xml_node meshNode = it->node();
			pugi::xml_node collisionNode = meshNode.parent().parent();
			pugi::xml_node linkNode = collisionNode.parent();

			const char *filename = meshNode.attribute("filename").value();

			std::cout << " Meshfile: " << filename << std::endl;

			std::vector<double> points = getPointsFromDAE(filename);
			std::cout << points.size() << " points pushed!" << std::endl;

			Eigen::Matrix4d tDAE = getTransformationMatriyFromDAE(filename);

			Fitter fitter = getCapsuleParams(points);
			std::vector<double> e1(3), e2(3);
			double radius;
			std::cout << "Initial: " << fitter.initParam() << std::endl;
			std::cout << "Solution: " << fitter.solutionParam() << std::endl;

			// assign values to e1,e2 and radius
			e1[0] = fitter.solutionParam()(0, 0);
			e1[1] = fitter.solutionParam()(1, 0);
			e1[2] = fitter.solutionParam()(2, 0);

			std::cout << "e1: 0,1,2: " << e1[0] << "," << e1[1] << "," << e1[2]
					<< std::endl;

			e2[0] = fitter.solutionParam()(3, 0);
			e2[1] = fitter.solutionParam()(4, 0);
			e2[2] = fitter.solutionParam()(5, 0);

			std::cout << "e2: 0,1,2: " << e2[0] << "," << e2[1] << "," << e2[2]
					<< std::endl;

			radius = fitter.solutionParam()(6, 0);

			std::cout << "rad: " << radius << std::endl;

//			std::cout << "link node: " << linkNode.name() << std::endl;
//			std::cout << "collision node: " << collisionNode.name() << std::endl;
			pugi::xml_node originNode = collisionNode.child("origin");
			const char* rpyStr = originNode.attribute("rpy").value();
			const char* xyzStr = originNode.attribute("xyz").value();
			const char* scaleStr = "";
			if (meshNode.attribute("scale")) {
				scaleStr = meshNode.attribute("scale").value();
			}

			if (linkNode.remove_child("collision")) {
				std::cout << "deleted old collision node!!" << std::endl;

				makeCapsuleFromXml(linkNode, e1, e2, radius, rpyStr, xyzStr,
						scaleStr, tDAE);
			}

			std::cout << std::endl;
		}
//		doc.save_file("./src/svh_new.urdf.xacro");
		doc.save_file("./src/svh_new.urdf");
//		doc.save_file("./src/test_new.urdf");
	} catch (std::exception& e) {
		std::cerr << "Unhandled Exception reached the top of main: " << e.what()
				<< ", application will now exit" << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

