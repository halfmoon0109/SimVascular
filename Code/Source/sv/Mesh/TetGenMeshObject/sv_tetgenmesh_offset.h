/* Copyright (c) Stanford University, The Regents of the University of
 *               California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/**
 * @file sv_tetgenmesh_offset.h
 * @brief The outer surface of the wall as the true offset of the interface at
 * the requested thickness, contoured from a distance field.
 *
 * Standard library only: no VTK, no TetGen. The one piece that needs a
 * volume mesher - the Delaunay tetrahedralization of the point cloud the
 * field is sampled on - is supplied by the caller.
 *
 * Why an offset and not an extrusion. Pushing every interface point out along
 * its normal gives every outer point one inner point, and that correspondence
 * cannot represent the wall at a concave junction: the wall of one vessel
 * fills the crotch up to the other vessel's wall, the outer surface there is
 * the crease where the two vessels' offsets meet, and the interface points in
 * the crotch have no outer point at all. Cutting the extrusion along the
 * creases where its sheets cross (sv_tetgenmesh_envelope.h) gets the crease
 * right where the sheets do cross, but where the extrusion folds or dips -
 * every crotch whose fillet is tighter than the wall, or whose two faces ask
 * for different thicknesses - the sheets stand a fraction of the wall above
 * the interface and nothing in the extrusion holds the surface that belongs
 * there. Measured on a carotid model, the wall over the interface at the four
 * junctions was 0.2 to 0.5 of the thickness asked for, whatever was done to
 * the thickness field or to the classification of the pieces.
 *
 * The offset has no correspondence to keep. The wall is the union over the
 * interface of the balls of radius t(y) about its points y (with t
 * interpolated over the triangles), and its outer surface is the zero level
 * of
 *
 *     d(x) = min over interface triangles T of (dist(x, T) - t at the closest point of T)
 *
 * on the wall side of the interface, taken negative in the lumen. Where two
 * vessels' walls meet, the minimum is the thicker one's and the septum is
 * solid; the crease and the rounding of a convex edge come out of the level
 * set rather than being aimed at, and every point of the surface is a wall
 * away from the interface by construction.
 *
 * Why not a grid. The first build of this contoured the field on a uniform
 * grid, and a grid that fits in memory is too coarse for the thin vessels of
 * a model whose walls range over an order of magnitude (measured: spacing 4.7
 * times the thinnest wall, at which the offset of the thin vessels collapsed
 * onto the interface). The resolution has to follow the interface. So the
 * field is sampled on a point cloud made of the interface points and their
 * offsets along the normals (see Options), that cloud
 * is tetrahedralized (Delaunay, by the caller), and the zero level is taken
 * by marching tetrahedra over that mesh: fine where the interface is fine,
 * coarse where it is coarse, with the contour points then put on the exact
 * zero level along their edges. The extruded layers run through each other at
 * a crotch just as the extrusion did, and it does not matter: they are
 * samples, not a surface.
 *
 * The interface is open at the vessel ends. The field alone would round the
 * offset off around each rim; a collar continuing the interface past each
 * cap rim keeps the offset a straight tube through the cap plane, and the
 * caller trims it back to the plane and closes the wall end with an annulus.
 */

#ifndef __SV_TETGENMESH_OFFSET_H__
#define __SV_TETGENMESH_OFFSET_H__

#include <string>
#include <vector>

namespace svoffset
{

/**
 * @brief The interface: the fluid/wall boundary, open at the vessel ends.
 */
struct Interface
{
  std::vector<double> points;        // three per point
  std::vector<double> normals;       // three per point, out of the lumen into the wall; normalized here
  std::vector<double> thickness;     // one per point, the wall thickness asked for there, positive
  std::vector<long long> triangles;  // three point ids per triangle, wound so that the right-hand normal points into the wall
};

/**
 * @brief What can be tuned.
 */
struct Options
{
  // The cloud is the interface points (value -t) and, along each normal, a
  // layer inside the wall at innerLayer times the local thickness, a layer
  // outside at outerLayer times it, and a far layer at farLayer times it but
  // no nearer than farSpacing times the local edge, so that the far points
  // of neighbouring rays are nearer each other than the surface and roof the
  // band: the Delaunay hull is then made of far points, all outside, and no
  // tetrahedron joins an inside point to the far side of the model. The
  // inner layer keeps the tetrahedra local (without it the interface points
  // of a thin vessel join the far side across the lumen and the contour gets
  // points in the lumen), and the outer layer keeps the zero level closely
  // bracketed (with the far layer alone as the outside, a thin coarse tube
  // came out with triangles crossing).
  double innerLayer = 0.5;
  double outerLayer = 1.5;
  double farLayer = 3.0;
  double farSpacing = 1.5;
  int snapIterations = 1;       // secant steps putting a contour point on the exact zero level (0 leaves the linear interpolation)
  double collapseRatio = 0.8;   // edges shorter than this times the local interface size are collapsed; zero or less leaves the contour as marched
  double collapseTurnCosine = 0.7;   // a collapse may not turn a triangle's normal by more than this cosine (slivers, whose normal means little, excepted)
  double collapseFieldTolerance = 0.1;   // no triangle a collapse makes may have its centre farther off the zero level than this fraction of the local thickness
};

/**
 * @brief The result: the offset surface, closed everywhere but where the
 * caller trims it.
 */
struct Surface
{
  std::vector<double> points;        // three per point
  std::vector<long long> triangles;  // three point ids per triangle, facing out of the wall
};

/**
 * @brief What happened, for the log and for judging the result.
 */
struct Report
{
  long long numInterfacePoints = 0;
  long long numInterfaceTriangles = 0;
  long long numRims = 0;                 // cap rims found on the interface
  long long numCollarTriangles = 0;      // added past the rims for the field
  double smallestThickness = 0.0;
  double largestThickness = 0.0;
  long long numCloudPoints = 0;
  long long numTetrahedra = 0;
  long long numTetrahedraCut = 0;        // by the zero level
  long long numFieldEvaluations = 0;
  long long numContourPoints = 0;        // as marched
  long long numContourTriangles = 0;
  long long numCollapsed = 0;            // edges collapsed by the decimation
  long long numPoints = 0;               // of the result
  long long numTriangles = 0;
  long long numBoundaryEdges = 0;        // of the result; none expected before the caller trims it
  long long numNonManifoldEdges = 0;     // edges on more than two triangles
  long long numMiswoundEdges = 0;        // edges traversed the same way by both triangles
  double secondsField = 0.0;
  double secondsDelaunay = 0.0;
  double secondsContour = 0.0;
  double secondsDecimate = 0.0;
  std::string firstFault;
  double firstFaultAt[3] = {0.0, 0.0, 0.0};
};

/**
 * @brief The Delaunay tetrahedralization the caller supplies.
 * @param points The cloud, three per point.
 * @param tetrahedra Set to four point ids per tetrahedron, in the cloud's
 * numbering; every cloud point should appear (a point the mesher dropped is a
 * gap in the sampling, not an error).
 * @param context Whatever the caller passed along.
 * @param error Set to why, on failure.
 * @return true if the tetrahedra were built.
 */
typedef bool (*DelaunayFunction)(const std::vector<double> &points,
    std::vector<long long> &tetrahedra, void *context, std::string &error);

/**
 * @brief Told what stage the build is at, so that a caller can log it as it
 * happens: a build that dies leaves the stage it died in.
 */
typedef void (*ProgressFunction)(const char *stage, void *context);

/**
 * @brief Builds the offset surface.
 * @param input The interface; see Interface.
 * @param options See Options.
 * @param delaunay The tetrahedralization to sample the field on.
 * @param context Passed to delaunay and progress.
 * @param surface Set to the offset surface.
 * @param report Set to what happened.
 * @param error Set to why, when the build fails outright.
 * @param progress Called at the start of each stage, if given.
 * @return 0 if the surface was built (it may still carry faults, counted in
 * the report), 1 if it could not be.
 */
int BuildOffsetSurface(const Interface &input, const Options &options,
    DelaunayFunction delaunay, void *context, Surface &surface, Report &report,
    std::string &error, ProgressFunction progress = nullptr);

/**
 * @brief The field the surface is the zero level of, for measuring: the
 * signed distance to the interface (with its collars) less the thickness
 * there. Built once and evaluated many times.
 */
class OffsetField
{
public:
  OffsetField();
  ~OffsetField();
  OffsetField(const OffsetField &) = delete;
  OffsetField &operator=(const OffsetField &) = delete;

  /**
   * @brief Builds the field's surface: the interface with a collar past each
   * cap rim, and the search structure.
   * @return 0 on success, 1 with error set otherwise.
   */
  int Build(const Interface &input, Report &report, std::string &error);

  /// The value at x; positive outside the wall, negative inside it and in the lumen.
  double Evaluate(const double x[3]) const;

  /// The interface where x is: the target edge length of the offset surface
  /// there (the interface's own edge, shortened where the curvature of the
  /// offset would otherwise put a chord more than a tenth of the wall off it)
  /// and the wall thickness there.
  void Local(const double x[3], double &size, double &thickness) const;

  /// How far from the field's surface a point is surely outside the wall.
  double Reach() const;

  /// The cap rims found: one loop of interface point ids each, in the triangles' winding.
  const std::vector<std::vector<long long> > &Rims() const;

  /// The field's surface: the interface followed by the collars.
  const std::vector<double> &Points() const;
  const std::vector<double> &Normals() const;
  const std::vector<double> &Thickness() const;
  const std::vector<long long> &Triangles() const;

  long long NumEvaluations() const;

private:
  struct Data;
  Data *data_;
};

/**
 * @brief Counts the edges of a triangle surface by how many triangles use
 * them, for the report and the tests.
 */
void CountEdges(const std::vector<long long> &triangles, long long &numBoundary,
    long long &numNonManifold, long long &numMiswound);

}  // namespace svoffset

#endif  // __SV_TETGENMESH_OFFSET_H__
