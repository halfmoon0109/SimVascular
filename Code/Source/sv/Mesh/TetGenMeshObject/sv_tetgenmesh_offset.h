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
  // The three layer points of one ray are moved off the ray sideways, each
  // by this fraction of its own distance out, in three directions 120
  // degrees apart (chosen by a hash of the point, so a run repeats). Four
  // points on one line are a degenerate input to the Delaunay: TetGen
  // returned tetrahedra of no volume on them (1e-17 against a median of
  // 1e-4 on a patch of the user's model), and the level's cut through such
  // a tetrahedron is a triangle lying in its plane, overlapping its
  // neighbours' and wound either way; every crossing of the raw contour on
  // that patch had one of them, and the snap and decimation then folded the
  // surface there. Zero leaves the points on the ray.
  double layerJitter = 0.1;
  int snapIterations = 16;      // at most this many secant steps (bisection when the secant leaves the bracket) putting a contour point on the zero level; 0 leaves the linear interpolation
  double snapTolerance = 1.0e-3;     // a contour point is on the level once |field| is within this fraction of the larger of its edge's end values (about the thickness)
  double collapseRatio = 0.8;   // edges shorter than this times the local interface size are collapsed; zero or less leaves the contour as marched
  double collapseTurnCosine = 0.7;   // a collapse may not turn a triangle's normal by more than this cosine (slivers, whose normal means little, excepted)
  double collapseFieldTolerance = 0.1;   // no triangle a collapse makes may have its centre farther off the zero level than this fraction of the local thickness
  // The target edge of the offset surface is the interface's own edge, or
  // shorter where a chord of that length would stand off a surface of the
  // interface's curvature plus the wall by more than this fraction of the
  // wall thickness. A layer surface inside the wall (a fraction of the
  // thickness) is decimated against the whole wall, so its caller raises this
  // by the reciprocal of the fraction.
  double chordTolerance = 0.05;
  // Every decimation operation (collapse, flip, valence-three and flat-apex
  // removal) is refused when a triangle it would make passes through a live
  // triangle near it, by the judgement of svenvelope::TrianglesCross. The
  // contour as marched is free of crossings, so the surface stays so.
  bool guardCrossings = true;
};

/**
 * @brief The result: the offset surface, closed everywhere but where the
 * caller trims it.
 */
struct Surface
{
  std::vector<double> points;        // three per point
  std::vector<long long> triangles;  // three point ids per triangle, facing out of the wall
  // Which vessel end each point belongs to, for trimming: the index into
  // rims of the cap whose collar, or whose interface within a collar's
  // length of its rim, is the nearest part of the field's surface to the
  // point; -1 for every other point. The part of the surface past a cap
  // plane that belongs to that cap is exactly the part its points own, so
  // the trim can leave a neighbouring vessel's surface alone even where it
  // runs past the plane within a rim radius - a spherical window around the
  // rim cannot (measured: two vessel ends 4.7 apart, radius 1.2 each, and
  // the window of one cut a loop out of the other's offset).
  std::vector<long long> pointRim;
  std::vector<std::vector<long long> > rims;   // the cap rims as loops of interface point ids, in the interface triangles' winding
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
  long long numZeroCloudPoints = 0;      // cloud points on the zero level itself (within rounding)
  long long numTetrahedra = 0;
  long long numTetrahedraCut = 0;        // by the zero level
  long long numFieldEvaluations = 0;
  long long numContourPoints = 0;        // as marched
  long long numContourTriangles = 0;
  long long numContourPointsOffLevel = 0;  // whose last evaluated position was still off the level by more than snapTolerance (the field jumps along the edge, or the steps ran out)
  long long numDegenerateContourTriangles = 0; // left out of the contour: naming a point twice, or emitted twice around points on the level
  long long numCollapsed = 0;            // edges collapsed by the decimation
  long long numRefusedForCrossing = 0;   // decimation operations refused because a triangle they would make passes through a live triangle near it
  long long numFoldedEdges = 0;          // of the result: edges whose two triangles still lie on each other within foldDegrees, which the volume mesher refuses
  double smallestFoldDegrees = 180.0;    // of the result: the smallest angle between the two triangles on an edge (180 is flat, 0 is folded)
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
 * @brief Builds the offset surfaces at several fractions of the wall
 * thickness at once: the surface at fraction f is the zero level of the
 * field with the thickness scaled by f (see OffsetField::Evaluate), so that
 * the surfaces of different fractions are nested. They are marched from
 * one point cloud and one field, so the marched contours cannot cross each
 * other either (the marched value falls with the fraction at every cloud
 * point); each is then decimated on its own. Built for the layered fill of
 * the wall, whose layer surfaces at 1/N, 2/N, ... of the thickness were
 * offsets of separately scaled interfaces before, each on its own cloud,
 * and crossed one another where the walls are thin (measured 2026-09-23 on
 * the user's 178k model: 17 crossings of the two-thirds surface through the
 * outer one at a branch root, and the mesher refused the shell).
 * @param fractions Each in (0, 1]; 1 is the outer surface.
 * @param surfaces One per fraction, in order.
 * @param reports One per fraction; the field and cloud figures are the same
 * in all of them.
 * @return 0 if every surface was built, 1 if not.
 */
int BuildOffsetSurfaces(const Interface &input, const Options &options,
    const std::vector<double> &fractions, DelaunayFunction delaunay, void *context,
    std::vector<Surface> &surfaces, std::vector<Report> &reports,
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
  int Build(const Interface &input, Report &report, std::string &error, double chordTolerance = 0.05);

  /// The value at x; positive outside the wall, negative inside it and in
  /// the lumen. With a thickness scale below one, the wall is that fraction
  /// of its thickness: the zero level is then the layer surface at that
  /// fraction, and the levels of different fractions never cross (the value
  /// falls as the fraction rises, at every point).
  double Evaluate(const double x[3], double thicknessScale = 1.0) const;

  /// The interface where x is: the target edge length of the offset surface
  /// there (the interface's own edge, shortened where the curvature of the
  /// offset would otherwise put a chord more than a twentieth of the wall
  /// off it), the wall thickness there, and the cap rim the field triangle
  /// x's offset stands on belongs to (a collar triangle, or an interface
  /// triangle within two collar lengths of the rim along the interface),
  /// -1 for none. The size and the thickness come from the nearest triangle
  /// by distance; the rim from the triangle the field takes its value from
  /// at x, the least distance less the wall, since that is the piece of the
  /// interface x's offset stands on.
  void Local(const double x[3], double &size, double &thickness, long long &rim) const;


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
 * @brief A cap plane the offset surface is trimmed along: the plane through
 * origin with the outward direction of the vessel end, and the rim (an
 * entry of Surface::pointRim) whose points past the plane the cut takes.
 */
struct CapPlane
{
  double origin[3];
  double outward[3];
  long long rim;
};

/**
 * @brief What TrimSurfaceAtCaps did, for the log.
 */
struct TrimReport
{
  long long numPointsBefore = 0, numTrianglesBefore = 0;
  long long numPointsAfter = 0, numTrianglesAfter = 0;
  std::vector<long long> numCut;      // per plane: points past it that came off
  std::vector<long long> numSnapped;  // per plane: points moved onto it
  long long numRemoved = 0;           // triangles wholly past a plane
  long long numDropped = 0;           // triangles left with nothing on the kept side but a point or an edge in the plane, or lying flat in it
  long long numSplit = 0;             // triangles the plane passed through
  long long numRimEdges = 0;
  long long numRimEdgesMerged = 0;   // rim edges shorter than a quarter of the mean edge, collapsed so the rim has no teeth the annulus stitching would fold on
  double shortestRimEdgeRatio = 0.0;  // shortest rim edge over the mean edge at its ends
  double smallestRimAngleDegrees = 0.0;  // smallest angle between a triangle on a rim and the cap plane, on the annulus side
};

/**
 * @brief Trims the offset surface at the cap planes: for each plane, the
 * points the surface hands to that plane's rim (Surface::pointRim) and lying
 * past the plane come off, the triangles the plane passes through are cut
 * along it, and the cut lands on the plane. A point of the rim's own within
 * snapFraction of its mean edge of the plane is moved onto the plane instead
 * of being cut a hair from, so that no edge shorter than a fraction of the
 * local edge and no triangle lying flat in the plane is left on the rim
 * (measured 2026-09-21: a VTK clip left edges of 1e-8 and flat triangles in
 * the plane there, which TetGen merged into degenerate facets and then
 * refused as facets folded onto the annulus). Points that belong to no
 * plane, or to another, are kept whatever side of the plane they lie on.
 * Points are compacted afterwards; pointRim is carried (-1 for a cut point).
 * The rims of the result are the boundary loops of its triangles.
 * @return 0 on success, 1 with error set otherwise.
 */
int TrimSurfaceAtCaps(Surface &surface, const std::vector<CapPlane> &planes, double snapFraction,
    TrimReport &report, std::string &error);

/**
 * @brief Counts the edges of a triangle surface by how many triangles use
 * them, for the report and the tests.
 */
void CountEdges(const std::vector<long long> &triangles, long long &numBoundary,
    long long &numNonManifold, long long &numMiswound);

/**
 * @brief The angle below which TetGen refuses two facets on one edge as
 * "nearly self-intersecting" (its -p/# tolerance, 0.1 degree by default).
 */
const double foldDegrees = 0.1;

/**
 * @brief Two triangles on one edge that lie on each other.
 */
struct FoldedPair
{
  long long a = -1, b = -1;              // the two triangles
  long long edge[2] = {-1, -1};          // the points of the edge they share
  double degrees = 0.0;                  // the angle between them: 0 is folded flat
};

/**
 * @brief Lists the edges on exactly two triangles, wound against each other,
 * whose triangles meet at an angle below maxDegrees: folded onto each other.
 * A crossing count does not see them (the two share an edge), and the volume
 * mesher refuses them below its tolerance (foldDegrees) as two facets
 * intersecting - measured 2026-09-23 on the user's model, two outer wall
 * triangles 0.03 degree apart, after every crossing count was zero.
 * Degenerate triangles (no normal) and edges on one or more than two
 * triangles are left to the other counts.
 * @param maxPairs How many pairs to keep, the most folded first; the count
 * returned is of all of them.
 * @param smallestDegrees Set to the smallest angle over every edge on two
 * triangles, folded or not (180 when there is none).
 * @return The number of folded edges.
 */
long long ListFoldedEdges(const std::vector<double> &points, const std::vector<long long> &triangles,
    double maxDegrees, size_t maxPairs, std::vector<FoldedPair> &pairs, double &smallestDegrees);

}  // namespace svoffset

#endif  // __SV_TETGENMESH_OFFSET_H__
