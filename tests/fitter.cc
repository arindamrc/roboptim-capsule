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

#define BOOST_TEST_MODULE fitter

#include <boost/test/unit_test.hpp>
#include <boost/test/output_test_stream.hpp>

#include <roboptim/capsule/util.hh>
#include <roboptim/capsule/fitter.hh>

#define BOOST_CHECK_SMALL_OR_CLOSE(EXP, OBS, TOL) \
    if (std::fabs (EXP) < TOL) { \
        BOOST_CHECK_SMALL(OBS, TOL); \
    } else { \
        BOOST_CHECK_CLOSE(EXP, OBS, TOL); \
    }

using boost::test_tools::output_test_stream;

BOOST_AUTO_TEST_CASE (fitter)
{
  using namespace roboptim::capsule;

  // Build a cubic polyhedron centered in (0,0,0).
  polyhedron_t polyhedron;
  value_type halfLength = 0.5;

  polyhedron.push_back (point_t (-halfLength, -halfLength, -halfLength));
  polyhedron.push_back (point_t (-halfLength, -halfLength, halfLength));
  polyhedron.push_back (point_t (-halfLength, halfLength, -halfLength));
  polyhedron.push_back (point_t (-halfLength, halfLength, halfLength));
  polyhedron.push_back (point_t (halfLength, -halfLength, -halfLength));
  polyhedron.push_back (point_t (halfLength, -halfLength, halfLength));
  polyhedron.push_back (point_t (halfLength, halfLength, -halfLength));
  polyhedron.push_back (point_t (halfLength, halfLength, halfLength));

  polyhedrons_t polyhedrons;
  polyhedrons.push_back (polyhedron);

  // Define initial capsule parameters. The segment must be inside the
  // polyhedron, and the capsule must contain the polyhedron.
  //
  //
  // To do so, compute initial guess by finding a bounding capsule
  // (not the minimum one).
  //
  // If needed, the convex hull of the polyhedron can be first computed
  // to reduce the number of constraints and accelerate the optimization
  // phase.
  polyhedrons_t convexPolyhedrons;
  computeConvexPolyhedron (polyhedrons, convexPolyhedrons);

  // Create fitter. It is used to find the best fitting capsule on the
  // polyhedron vector.
  Fitter fitter_cube (convexPolyhedrons);

  point_t endPoint1 = point_t (0., 0., 0.);
  point_t endPoint2 = point_t (0., 0., 0.);
  value_type radius = 0.;
  computeBoundingCapsulePolyhedron (convexPolyhedrons,
				    endPoint1, endPoint2, radius);

  argument_t initParam (7);
  convertCapsuleToSolverParam (initParam, endPoint1, endPoint2, radius);

  // Compute best fitting capsule.
  fitter_cube.computeBestFitCapsule (initParam);
  argument_t solutionParam = fitter_cube.solutionParam ();
  std::cout << fitter_cube << std::endl;

  double epsilon = 1e-1;
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[0], 0.190393847,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[1], 0.,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[2], 0.,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[3], -0.190393847,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[4], 0.,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[5], 0.,epsilon);
  BOOST_CHECK_SMALL_OR_CLOSE(solutionParam[6], 0.77191705555821011, epsilon)

  polyhedrons.clear ();
  convexPolyhedrons.clear ();

  // Enlarge the cube to get a rectangular box
  int n = 4;
  polyhedron[0][0] -= n * halfLength;
  polyhedron[1][0] -= n * halfLength;
  polyhedron[2][0] -= n * halfLength;
  polyhedron[3][0] -= n * halfLength;
  polyhedron[4][0] += n * halfLength;
  polyhedron[5][0] += n * halfLength;
  polyhedron[6][0] += n * halfLength;
  polyhedron[7][0] += n * halfLength;
  polyhedrons.push_back (polyhedron);
  computeConvexPolyhedron (polyhedrons, convexPolyhedrons);

  Fitter fitter_rect (convexPolyhedrons);
  computeBoundingCapsulePolyhedron (convexPolyhedrons,
				    endPoint1, endPoint2, radius);
  convertCapsuleToSolverParam (initParam, endPoint1, endPoint2, radius);
  fitter_rect.computeBestFitCapsule (initParam);
  std::cout << fitter_rect << std::endl;
}
