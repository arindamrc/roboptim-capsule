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
#include "util.cc"

using namespace roboptim;
using namespace roboptim::capsule;

/*******************************************************************************************************************************************************/

class BestFitCapsule {
public:
	Eigen::Vector3d endpoint1;
	Eigen::Vector3d endpoint2;
	double radius;
	double getLength();
	Eigen::Vector3d getMidpoint();
	Eigen::Vector3d getRotation();
	BestFitCapsule(std::vector<double> points, Eigen::Vector3d translation,
			Eigen::Vector3d rotation, Eigen::Vector3d scale,
			Eigen::Matrix4d transformation);
private:
	void makeCapsule();
	Fitter getCapsuleParams();
	void transformParameters();
	Eigen::Vector3d xyz_i, rpy_i, s_i;
	Eigen::Matrix4d M_i;
	std::vector<double> points;
	double len;
	Eigen::Vector3d midpoint;
	Eigen::Vector3d rot;
};

double BestFitCapsule::getLength() {
	return len;
}

Eigen::Vector3d BestFitCapsule::getMidpoint() {
	return midpoint;
}

Eigen::Vector3d BestFitCapsule::getRotation() {
	return rot;
}

BestFitCapsule::BestFitCapsule(std::vector<double> points,
		Eigen::Vector3d translation, Eigen::Vector3d rotation,
		Eigen::Vector3d scale, Eigen::Matrix4d transformation) {
	radius = 0;
	xyz_i = translation;
	rpy_i = rotation;
	s_i = scale;
	M_i = transformation;
	this->points = points;
	makeCapsule();
}

void BestFitCapsule::makeCapsule() {
	Fitter fitter = getCapsuleParams();

	// assign values to e1,e2 and radius
	this->endpoint1 = Eigen::Vector3d(fitter.solutionParam()(0, 0),
			fitter.solutionParam()(1, 0), fitter.solutionParam()(2, 0));
	this->endpoint2 = Eigen::Vector3d(fitter.solutionParam()(3, 0),
			fitter.solutionParam()(4, 0), fitter.solutionParam()(5, 0));
	this->radius = fitter.solutionParam()(6, 0);

	transformParameters();

}

void BestFitCapsule::transformParameters() {
	Eigen::Vector4d homoe1 = endpoint1.homogeneous();
	Eigen::Vector4d homoe2 = endpoint2.homogeneous();

	Eigen::Matrix4d transform =
			roboptim::geometry::create_transformation_matrix(rpy_i, xyz_i, s_i);

	homoe1 = transform * M_i * homoe1;
	homoe2 = transform * M_i * homoe2;
	radius = radius * M_i(0, 0);

	endpoint1 = homoe1.hnormalized();
	endpoint2 = homoe2.hnormalized();

	Eigen::Vector3d capAxis = endpoint2 - endpoint1;

	len = (endpoint1 - endpoint2).norm();
	midpoint = (endpoint1 + endpoint2) * 0.5;

	rot = roboptim::geometry::axisAngles(capAxis);
}

/* Get the capsule endpoints and radius using roboptim capsule*/
Fitter BestFitCapsule::getCapsuleParams() {
	assert(points.size() > 0 && "Cannot compute capsule for point set.");

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

/*******************************************************************************************************************************************************/

class Urdf {
public:
	const char* getNextMesh();
	bool hasMoreMesh();
	Urdf(const char* filepath, bool xacro);
	void replaceMesh(BestFitCapsule capsule);
	void replaceMeshWithXacro(BestFitCapsule capsule);
	Eigen::Vector3d getMeshRotation();
	Eigen::Vector3d getMeshTranslation();
	Eigen::Vector3d getMeshScale();
	bool save(const char* path);
	bool isXacro();
private:
	const char* filepath;
	const std::string collision_path = "/robot/link[*]/collision/geometry/mesh";
	pugi::xpath_node_set nodes;
	const pugi::xpath_node* it;
	pugi::xml_node lastLinkNode;
	pugi::xml_node lastCollisionNode;
	pugi::xml_node lastMeshNode;
	pugi::xml_document doc;
	void addXacroDef();
	bool xacro;
};

Urdf::Urdf(const char* filepath, bool xacro) {
	this->filepath = filepath;
	this->xacro = xacro;
	doc.load_file(filepath);
	nodes = doc.select_nodes(collision_path.c_str());
	it = nodes.begin();
//	pugi::xml_node meshNode = it->node();
	if (xacro) {
		addXacroDef();
	}
}

bool Urdf::save(const char* path) {
	return doc.save_file(path);
}

void Urdf::addXacroDef() {
	pugi::xml_node robotNode = doc.child("robot");

	pugi::xml_node xacroNode = robotNode.append_child("xacro:macro");
	xacroNode.append_attribute("name").set_value("capsule-collision");
	xacroNode.append_attribute("params").set_value("e1 e2 rad len rpy mid");

	// The first sphere
	pugi::xml_node s1Collision = xacroNode.append_child("collision");
	pugi::xml_node s1Origin = s1Collision.append_child("origin");
	s1Origin.append_attribute("rpy").set_value("${rpy}");
	s1Origin.append_attribute("xyz").set_value("${e1}");

	pugi::xml_node s1Geom = s1Collision.append_child("geometry");
	pugi::xml_node s1 = s1Geom.append_child("sphere");
	pugi::xml_attribute rs1 = s1.append_attribute("radius");
	rs1.set_value("${rad}");

	// The middle cylinder
	pugi::xml_node crCollision = xacroNode.append_child("collision");
	pugi::xml_node crOrigin = crCollision.append_child("origin");
	crOrigin.append_attribute("rpy").set_value("${rpy}");
	crOrigin.append_attribute("xyz").set_value("${mid}");

	pugi::xml_node crGeom = crCollision.append_child("geometry");
	pugi::xml_node cyl = crGeom.append_child("cylinder");
	pugi::xml_attribute crl = cyl.append_attribute("length");
	crl.set_value("${len}");
	pugi::xml_attribute crr = cyl.append_attribute("radius");
	crr.set_value("${rad}");

	// The second sphere
	pugi::xml_node s2Collision = xacroNode.append_child("collision");
	pugi::xml_node s2Origin = s2Collision.append_child("origin");
	s2Origin.append_attribute("rpy").set_value("${rpy}");
	s2Origin.append_attribute("xyz").set_value("${e2}");

	pugi::xml_node s2Geom = s2Collision.append_child("geometry");
	pugi::xml_node s2 = s2Geom.append_child("sphere");
	pugi::xml_attribute s2r = s2.append_attribute("radius");
	s2r.set_value("${rad}");
}

bool Urdf::isXacro() {
	return this->xacro;
}

const char* Urdf::getNextMesh() {
	lastMeshNode = it->node();
	lastCollisionNode = lastMeshNode.parent().parent();
	lastLinkNode = lastCollisionNode.parent();

	const char *filename = lastMeshNode.attribute("filename").value();
	it++;
	return filename;
}

bool Urdf::hasMoreMesh() {
	return it != nodes.end();
}

void Urdf::replaceMesh(BestFitCapsule capsule) {
	//remove old collision element
	if (!lastLinkNode.remove_child("collision")) {
		std::cout << "Can't remove collision element. Skipping..." << std::endl;
	}
	// The first sphere
	pugi::xml_node s1Collision = lastLinkNode.append_child("collision");
	pugi::xml_node s1Origin = s1Collision.append_child("origin");
	Eigen::Vector3d rot = capsule.getRotation();
	Eigen::Vector3d endpoint1 = capsule.endpoint1;
	Eigen::Vector3d endpoint2 = capsule.endpoint2;
	Eigen::Vector3d midpoint = capsule.getMidpoint();

	std::stringstream s1rpySS;
	s1rpySS << rot(0) << " " << rot(1) << " " << rot(2);
	s1Origin.append_attribute("rpy").set_value(s1rpySS.str().c_str());

	std::stringstream s1xyzSS;
	s1xyzSS << (endpoint1(0)) << " " << (endpoint1(1)) << " " << (endpoint1(2));
	s1Origin.append_attribute("xyz").set_value(s1xyzSS.str().c_str());

	pugi::xml_node s1Geom = s1Collision.append_child("geometry");
	pugi::xml_node s1 = s1Geom.append_child("sphere");
	pugi::xml_attribute rs1 = s1.append_attribute("radius");
	rs1.set_value(capsule.radius);

	// The middle cylinder
	pugi::xml_node crCollision = lastLinkNode.append_child("collision");
	pugi::xml_node crOrigin = crCollision.append_child("origin");

	std::stringstream crrpySS;
	crrpySS << rot(0) << " " << rot(1) << " " << rot(2);
	crOrigin.append_attribute("rpy").set_value(crrpySS.str().c_str());

	std::stringstream crxyzSS;
	crxyzSS << (midpoint(0)) << " " << (midpoint(1)) << " " << (midpoint(2));
	crOrigin.append_attribute("xyz").set_value(crxyzSS.str().c_str());

	pugi::xml_node crGeom = crCollision.append_child("geometry");
	pugi::xml_node cyl = crGeom.append_child("cylinder");
	pugi::xml_attribute crl = cyl.append_attribute("length");
	crl.set_value(capsule.getLength());
	pugi::xml_attribute crr = cyl.append_attribute("radius");
	crr.set_value(capsule.radius);

	// The second sphere
	pugi::xml_node s2Collision = lastLinkNode.append_child("collision");
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
	s2r.set_value(capsule.radius);
}

void Urdf::replaceMeshWithXacro(BestFitCapsule capsule) {
	//remove old collision element
	if (!lastLinkNode.remove_child("collision")) {
		std::cout << "Can't remove collision element. Skipping..." << std::endl;
	}
	pugi::xml_node xacroCollision = lastLinkNode.append_child(
			"xacro:capsule-collision");
	Eigen::Vector3d rot = capsule.getRotation();
	Eigen::Vector3d endpoint1 = capsule.endpoint1;
	Eigen::Vector3d endpoint2 = capsule.endpoint2;
	Eigen::Vector3d midpoint = capsule.getMidpoint();

	std::stringstream e1SS;
	e1SS << (endpoint1(0)) << " " << (endpoint1(1)) << " " << (endpoint1(2));
	xacroCollision.append_attribute("e1").set_value(e1SS.str().c_str());

	std::stringstream e2SS;
	e2SS << (endpoint2(0)) << " " << (endpoint2(1)) << " " << (endpoint2(2));
	xacroCollision.append_attribute("e2").set_value(e2SS.str().c_str());

	xacroCollision.append_attribute("rad").set_value(capsule.radius);

	xacroCollision.append_attribute("len").set_value(capsule.getLength());

	std::stringstream rpySS;
	rpySS << (rot(0)) << " " << (rot(1)) << " " << (rot(2));
	xacroCollision.append_attribute("rpy").set_value(rpySS.str().c_str());

	std::stringstream midSS;
	midSS << (midpoint(0)) << " " << (midpoint(1)) << " " << (midpoint(2));
	xacroCollision.append_attribute("mid").set_value(midSS.str().c_str());
}

Eigen::Vector3d Urdf::getMeshRotation() {
	pugi::xml_node originNode = lastCollisionNode.child("origin");
	std::vector<double> rpy;
	if (!originNode.attribute("rpy")) {
		return Eigen::Vector3d(0, 0, 0);
	}
	const char* rpyStr = originNode.attribute("rpy").value();
	roboptim::geometry::getAsVector(rpyStr, &rpy);
	return Eigen::Vector3d(rpy.data());
}

Eigen::Vector3d Urdf::getMeshTranslation() {
	pugi::xml_node originNode = lastCollisionNode.child("origin");
	std::vector<double> xyz;
	if (!originNode.attribute("xyz")) {
		return Eigen::Vector3d(0, 0, 0);
	}
	const char* xyzStr = originNode.attribute("xyz").value();
	roboptim::geometry::getAsVector(xyzStr, &xyz);
	return Eigen::Vector3d(xyz[0], xyz[1], xyz[2]);
}

Eigen::Vector3d Urdf::getMeshScale() {
	std::vector<double> scaleVec;
	const char* scaleStr = "";
	if (!lastMeshNode.attribute("scale")) {
		return Eigen::Vector3d(1, 1, 1);
	}
	scaleStr = lastMeshNode.attribute("scale").value();
	std::cout << "scaleStr: " << scaleStr << std::endl;
	roboptim::geometry::getAsVector(scaleStr, &scaleVec);
	return Eigen::Vector3d(scaleVec[0], scaleVec[1], scaleVec[2]);
}

/*******************************************************************************************************************************************************/

class Dae {
private:
	const char* filename;
	pugi::xml_document doc;
public:
	Dae(const char* filename);
	std::vector<double> getPoints();
	Eigen::Matrix4d getTransformationMatrix();
};

Dae::Dae(const char* filename) {
	this->filename = filename;
	this->doc.load_file(filename);
}

/* Extract all polygon vertices from DAE file*/
std::vector<double> Dae::getPoints() {
	std::vector<double> points;
	const std::string path =
			"/COLLADA/library_geometries/geometry/mesh/source[*]/float_array[contains(@id,'positions')]";
	const pugi::xpath_node_set nodes = doc.select_nodes(path.c_str());

	for (const pugi::xpath_node* it = nodes.begin(); it != nodes.end(); ++it) {
		const char* vertexString = it->node().text().get();
		roboptim::geometry::getAsVector(vertexString, &points);
	}
	return points;
}

/* Get transformation matrix from DAE file*/
Eigen::Matrix4d Dae::getTransformationMatrix() {
	Eigen::Matrix4d t;
	std::vector<double> els;
	const std::string path =
			"/COLLADA/library_visual_scenes/visual_scene[@id=\"Scene\"]/node/matrix[@sid=\"transform\"]";
	const pugi::xpath_node_set nodes = doc.select_nodes(path.c_str());

	for (const pugi::xpath_node* it = nodes.begin(); it != nodes.end(); ++it) {
		const char* elString = it->node().text().get();
		roboptim::geometry::getAsVector(elString, &els);
	}
	t = Eigen::Matrix4d(els.data());
	return t;
}

int main(int argc, char **argv) {
	try {
		const char* filename = "./src/svh.urdf";
		Urdf urdf(filename, true);
		while (urdf.hasMoreMesh()) {
			const char* meshfile = urdf.getNextMesh();
			std::cout << "collada file: " << meshfile << std::endl;
			Dae dae(meshfile);
			std::vector<double> points = dae.getPoints();
			Eigen::Matrix4d M = dae.getTransformationMatrix();
			std::cout << "transformation matrix retrieved." << std::endl;
			Eigen::Vector3d t = urdf.getMeshTranslation();
			std::cout << "t: " << t << std::endl;
			std::cout << "mesh translation retrieved." << std::endl;
			Eigen::Vector3d r = urdf.getMeshRotation();
			std::cout << "r: " << r << std::endl;
			std::cout << "mesh rotation retrieved." << std::endl;
			Eigen::Vector3d s = urdf.getMeshScale();
			std::cout << "scaleVec: " << s << std::endl;
			BestFitCapsule capsule = BestFitCapsule(points, t, r, s, M);
			if (urdf.isXacro()) {
				urdf.replaceMeshWithXacro(capsule);
			} else {
				urdf.replaceMesh(capsule);
			}
		}
		if (urdf.save("./src/svh_new_2.urdf")) {
			std::cout << "Capsules generated successfully!";
		}
	} catch (std::exception& e) {
		std::cerr << "Unhandled Exception reached the top of main: " << e.what()
				<< ", application will now exit" << std::endl;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

