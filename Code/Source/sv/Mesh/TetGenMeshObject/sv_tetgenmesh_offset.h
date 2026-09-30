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
#include <unordered_map>
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
  long long numUnfolded = 0;             // edges whose two triangles lay on each other within a degree, flipped, collapsed or taken out with the spike point they stand on by the decimation's last pass
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

  /// The triangle of the field's surface x's offset stands on: the one whose
  /// term, its distance less the wall there at the scale, is the least at x,
  /// the piece of the interface the field takes its value from (a collar
  /// triangle counts; its id is past the interface's). -1 if none.
  long long Owner(const double x[3], double thicknessScale = 1.0) const;


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
  /// The search behind Local and Owner: the nearest triangle by distance and
  /// the one of the least term (distance less the scaled wall) at x.
  void Search(const double x[3], double thicknessScale, long long &nearest, double &distance,
      long long &owner, double &term) const;
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
 * @brief What decides where the wall can be prism layers extruded along the
 * interface's normals (the structured zone) and where it must be filled
 * with tetrahedra between the interface and the offset surfaces (the
 * junction zone). See docs/wall-mesh-purpose.md section 8.
 */
struct ZoneOptions
{
  // A layer point p + f t n stands on its level when the field at it, at
  // the scale f, is within this fraction of f t of zero. It cannot be above
  // (its own triangle is within f t of it); it is below where another part
  // of the wall covers it - the collar a thick parent's offset makes around
  // a thin branch's root, or a fold where the normals of a concave stretch
  // meet within the wall.
  double levelTolerance = 0.1;
  // The top of a prism must face the way its base does and keep this
  // fraction of the base's area, at every layer.
  double minTopArea = 0.1;
  // The junction zone grows by this many rings of triangles (sharing a
  // point) around what the checks reject, so that the tetrahedra have room
  // and the zipper between the two zones' surfaces stays off the trouble.
  // Three rather than two since the 178k model (2026-09-29): the piece of
  // a layer surface over the annular junction of a thin branch on a thick
  // parent was two triangles wide with two, too narrow to keep its shape
  // through the erosion.
  int marginRings = 3;
  // A structured island of fewer triangles than this joins the junction zone.
  long long minIsland = 50;
};

struct ZoneReport
{
  long long numTriangles = 0;
  long long numStructured = 0;             // prism layers
  long long numPointsCovered = 0;          // a layer point of theirs inside another part of the wall
  long long numPointsOffLevel = 0;         // a layer point off its level the other way (should not happen)
  long long numTrianglesInverted = 0;      // a layer triangle facing away from the base, or too small
  long long numTrianglesCrossing = 0;      // tops of otherwise valid prisms passing through each other
  long long numTrianglesMargin = 0;        // added to the junction zone by the margin rings
  long long numIslands = 0;                // structured islands too small to keep
  long long numIslandTriangles = 0;
  long long numPinches = 0;                // points the zone boundary passed through twice, opened up
  long long numTrianglesPinched = 0;       // given to the junction zone at them
  long long numJunctionRegions = 0;        // connected pieces of the junction zone
  long long numEvaluations = 0;
};

/**
 * @brief Decides, triangle by triangle of the interface, whether the wall
 * over it can be prism layers along the normals (structured, 1) or has to
 * be filled with tetrahedra (0): the layer points p + (k/N) t n of its
 * corners must all stand on their levels of the field, the layer triangles
 * must face the way the base does, and the tops of the prisms kept must not
 * pass through one another; the junction zone is then widened by the
 * margin rings and small structured islands are given up.
 * @param numLayers N, the layers through the thickness.
 * @param structured One per interface triangle.
 * @return 0 on success, 1 with error set otherwise.
 */
int ClassifyPrismZone(const Interface &input, const OffsetField &field, int numLayers,
    const ZoneOptions &options, std::vector<unsigned char> &structured, ZoneReport &report, std::string &error);

/**
 * @brief The prism layers over the structured zone, as tetrahedra (section 8
 * step 2). The points are the interface's, in its order, followed by the
 * layer points p + (k/N) t n that a structured triangle touches, allocated
 * as met (layerPoint says which); every other layer point does not exist.
 * Each prism is three tetrahedra whose side diagonals run from the corner
 * of the smallest interface id (so neighbouring prisms share their side
 * triangles and no prism is left with a cyclic set of diagonals), all wound
 * with positive volume. The boundary faces are wound to face out of the
 * tetrahedron they belong to.
 */
struct PrismMesh
{
  int numLayers = 0;
  std::vector<double> points;                 // three per point
  std::vector<long long> layerPoint;          // (k-1)*numInterfacePoints + i -> point id of p_i + (k/N) t n, -1 if not made
  std::vector<long long> tetrahedra;          // four point ids each
  std::vector<long long> tetrahedronLayer;    // the layer (1..N) each tetrahedron lies in
  std::vector<long long> topTriangles;        // the outer surface of the zone: three ids, facing out of the wall
  std::vector<std::vector<long long> > layerTriangles;   // [k-1]: the zone's layer k surface (k = N is the tops), facing out
  std::vector<long long> sideTriangles;       // the zone's walls toward the junction zone: three ids, facing the junction zone
  std::vector<long long> sideTriangleLayer;   // the layer (1..N) of each side triangle
  std::vector<long long> rimTriangles;        // the zone's walls at the cap rims (the interface's own boundary): three ids, facing out
  // The zone boundary on the interface, as directed edges (a, b) of
  // interface point ids in the structured triangle's winding, one entry per
  // boundary edge. The zone's layer surfaces traverse these edges the same
  // way (they are wound like the interface, whose normal points into the
  // wall, which is out of the wall's outer side), so a strip zipped to them
  // must take them the other way round: ZipChains wants the chain reversed.
  std::vector<long long> zoneBoundaryEdges;
};

/**
 * @brief Builds the prism layers over the structured triangles (section 8
 * step 2): N stacked prisms per triangle, each three tetrahedra.
 * @param structured One per interface triangle, as ClassifyPrismZone leaves it.
 * @return 0 on success, 1 with error set otherwise.
 */
int BuildPrismLayers(const Interface &input, const std::vector<unsigned char> &structured, int numLayers,
    PrismMesh &out, std::string &error);

/**
 * @brief What is left of a layer surface over the junction zone (section 8
 * step 3), and how much of it went.
 */
struct ZonePieceReport
{
  long long numTriangles = 0;          // of the level surface
  long long numStructuredOwned = 0;    // dropped: their centre's offset stands on a structured triangle (or a collar of a structured rim)
  long long numOverPrisms = 0;         // dropped: a corner or the centre of theirs, lowered a little, lies inside a prism tetrahedron of this layer (the triangle hangs over the prism zone)
  long long numEroded = 0;             // dropped: within the erosion rings of a dropped triangle
  long long numCrossing = 0;           // dropped: passing through the prism zone's layer surface
  long long numEars = 0;               // dropped: on two boundary edges, so that the boundary the zipper follows is not jagged
  long long numHoles = 0;              // dropped patches enclosed by kept triangles: a hole in the piece leaves the shell open, so the prisms that own or pass through it have to go to the junction zone
  long long numHoleTriangles = 0;      // in them
  long long numScraps = 0;             // dropped: kept pieces of fewer than a dozen triangles, which no zone loop could be zipped to
  long long numScrapTriangles = 0;
  long long numPinches = 0;            // points where the piece's boundary passed twice (two fans of kept triangles around them), opened by dropping all but the largest fan
  long long numPinchTriangles = 0;
  long long numKept = 0;
  long long numChains = 0;             // boundary chains of the piece
};

/**
 * @brief Keeps the triangles of a layer surface whose centre's offset stands
 * on a junction-zone triangle of the interface (by OffsetField::Owner at the
 * level's fraction; a collar triangle counts as its rim's zone), takes off
 * the given rings around what was dropped so that the piece ends short of
 * the prism zone's surface rather than over it, and takes off whatever still
 * passes through that surface (guardPoints/guardTriangles: the prism zone's
 * layer surface at the same level), then the ears - triangles on two
 * boundary edges - for the given passes, so that the boundary the zipper
 * follows has no notches for a strip triangle to fold over (measured on the
 * synthetic junction: one notch, four crossings). The piece's points are
 * compacted.
 * @param rimStructured One per cap rim: whether the interface triangles at
 * that rim are all structured, for the collar-owned centres.
 * @return 0 on success, 1 with error set otherwise.
  * @param loopSegments The lifted zone boundary edges at this level, six
 * numbers each (a, b): given, the erosion is by distance instead of rows -
 * the kept triangles whose centre lies within erosionDistance times the
 * nearest segment's length of it go (peeled one by one as the rows are,
 * so that the piece keeps its shape). Rows are the level's triangles,
 * which over the collar of a thick parent around a thin branch are few
 * where the interface's are many (the level surface is compressed in the
 * crotch), and two rows ate the whole band there (measured 2026-09-29 on
 * the 178k model, a branch 0.38 thick on a parent 0.85 thick); the
 * distance follows the loop's own edge instead.
 * @param overPrismTetrahedra Set, if given, per triangle dropped for
 * hanging over the prism zone to the index (into layerTetrahedra) of the
 * tetrahedron that held its point: a patch of them inside a junction
 * region's piece is a hole where a structured triangle's prisms stand in
 * the way (a branch root passing through its parent's prisms beyond the
 * margin), and those prisms' triangles have to go to the junction zone.
 * @param layerTetrahedra The prism tetrahedra of this layer (four ids into
 * guardPoints each): a kept triangle whose corners or centre, lowered by
 * a twentieth of the thickness, lie inside one hangs over the prism zone
 * and goes, whatever owns it (its strip would pass through the prism
 * wall; measured 2026-09-29 on the 178k model, 108 strip-wall crossings.
 * The wall's plane was tried for this and cut the pieces off around every
 * curved vessel).
*/
int TrimSurfaceToZone(const Surface &level, double fraction, const OffsetField &field, long long numInterfaceTriangles,
    const std::vector<unsigned char> &structured, const std::vector<unsigned char> &rimStructured, int erosionRings,
    const std::vector<double> &guardPoints, const std::vector<long long> &guardTriangles, int earPasses,
    Surface &piece, ZonePieceReport &report, std::string &error,
    std::vector<long long> *crossingGuards = nullptr, std::vector<long long> *holeOwners = nullptr,
    const std::vector<long long> *layerTetrahedra = nullptr, const std::vector<double> *loopSegments = nullptr,
    double erosionDistance = 1.0, std::vector<long long> *overPrismTetrahedra = nullptr);

/**
 * @brief The boundary of a triangle set as chains of point ids, each in the
 * direction the triangles traverse its edges: a closed chain repeats no
 * point (closed[c] is 1), an open one runs from a point of one boundary
 * end to the other. A point on more than two boundary edges ends the chains
 * there.
 */
void BoundaryChains(const std::vector<long long> &triangles, std::vector<std::vector<long long> > &chains,
    std::vector<unsigned char> &closed);

/**
 * @brief Zips two chains of point ids over one point array into a strip of
 * triangles: both run the same way, and the strip takes each of its edges
 * the way chain A gives them (a_i -> a_i+1) and the edges of chain B
 * backwards (b_j+1 -> b_j), so that a surface traversing A's edges backwards
 * and one traversing B's edges forwards are wound consistently with it. At
 * each step the shorter of the two possible diagonals is taken. Closed
 * chains are started at the nearest pair of points and go round once; open
 * chains go from their first points to their last.
 * @return 0 on success, 1 with error set otherwise.
 */
int ZipChains(const std::vector<double> &points, const std::vector<long long> &chainA,
    const std::vector<long long> &chainB, bool closed, std::vector<long long> &triangles, std::string &error,
    double *agreement = nullptr);

/**
 * @brief The closed surface around the junction zone's volume (section 8
 * step 4), for the volume mesher: the zone's interface triangles reversed,
 * the prism zone's walls toward it, and at every layer the layer surface's
 * piece over the zone zipped to the prism zone's layer ring. The points are
 * the prism mesh's (the interface and its layer points, in its order)
 * followed by the pieces' at each level.
 */
// One zipper strip of the junction shell: a zone boundary loop at a layer
// and the piece chain it was zipped to, for the log.
struct ZipReport
{
  int level = 0;
  double centre[3] = {0.0, 0.0, 0.0};   // of the lifted loop
  long long loopPoints = 0;
  long long chainPoints = 0;
  double meanDistance = 0.0;            // from the loop's points to the chain
  double agreement = 0.0;               // of the two chains' directions (ZipChains), positive when they run the same way
  long long stripTriangles = 0;
  long long crossingTriangles = 0;      // of the strip's, passing through another shell triangle
  long long foldedEdges = 0;            // edges of the strip's triangles whose two triangles lie on each other
  std::vector<long long> loop;          // the zone boundary loop's interface points, for WidenJunctionZone when the strip is at fault
};

struct JunctionShell
{
  std::vector<double> points;            // only the points the shell's triangles use
  std::vector<long long> prismPoint;     // per shell point: its id in the prism mesh, -1 for a point of a layer surface's piece
  std::vector<long long> triangles;
  std::vector<int> markers;              // per triangle: 1 the interface, 2 the outer surface, 100+k the layer k surface, 300+k the prism zone's wall at layer k
  std::vector<ZonePieceReport> pieces;   // [k-1]
  std::vector<ZipReport> zips;           // one per loop and level
  // The bands (the pieces with their zipper strips) relaxed before the
  // check: edge flips that raise the smallest angle and a smoothing of the
  // pieces' points projected back onto the level, the prism rings fixed.
  long long numBandFlips = 0;
  long long numBandPointsMoved = 0;
  long long numBandTriangles = 0;
  long long numBandUnder10Before = 0;    // band triangles with an angle under 10 degrees, before and after
  long long numBandUnder10After = 0;
  double bandSmallestAngleBefore = 180.0;
  double bandSmallestAngleAfter = 180.0;
  bool relaxationReverted = false;       // the relaxation made the shell cross itself or fold and was undone whole (after the local undo below failed)
  long long numRelaxationUndone = 0;     // triangles whose relaxation was undone locally (their points put back, their flips undone) for crossing or folding
  long long relaxationCrossingsAfter = 0, relaxationFoldsAfter = 0;   // what the relaxed shell had before the undo, against the shell's own counts before
  long long numBandCollapses = 0;        // short edges of the bands collapsed (two triangles gone each)
  long long numBandVerticesRemoved = 0;  // free points on three triangles taken out (two triangles gone each)
  long long numBandTrianglesRemoved = 0;
  bool collapseReverted = false;         // the collapses and removals made the shell cross or fold and were undone whole
  long long numCappedChains = 0;         // short boundary chains of the pieces that no zone loop takes (a tunnel's mouth where two walls nearly touch, or a hole left by the trim), closed with a fan of triangles each
  long long numCapTriangles = 0;
  long long numCrossingTriangles = 0;    // shell triangles passing through another (svenvelope::CountCrossingTriangles over the whole shell)
  long long numFoldedEdges = 0;          // shell edges whose two triangles lie on each other within foldDegrees (ListFoldedEdges over the whole shell)
  long long numZipperTriangles = 0;
  long long numJunctionTriangles = 0;    // interface triangles of the zone
  long long numJunctionOnRims = 0;       // of them, on a cap rim: not handled yet, the build refuses them
  std::vector<long long> failedLoopPoints;   // when the build fails for a zone boundary loop without a piece chain: the loop's interface points, for WidenJunctionZone
  std::vector<long long> crossingStructuredTriangles;   // structured interface triangles whose layer surface passed through a piece at some level (the piece's triangles there were dropped, which can leave a hole in it)
  std::vector<long long> holeOwnerTriangles;   // structured interface triangles owning a triangle of a hole in a piece (a dropped patch enclosed by kept ones) at some level
  std::vector<long long> overPrismStructuredTriangles;   // structured interface triangles whose prisms held a point of a piece triangle dropped for hanging over the prism zone, at some level
  std::vector<long long> failedTriangles;    // when the build fails for a piece chain no loop takes (a hole): the structured triangles owning or passing through the piece near it, for WidenJunctionZone
};

/**
 * @brief Assembles the junction shell (section 8 step 4). The level surfaces
 * come in the order of the layers (the last is the outer surface), already
 * trimmed at the caps; each is trimmed to the zone (TrimSurfaceToZone) and
 * zipped to the prism zone's ring at its layer. A junction zone touching a
 * cap rim is refused for now. The shell's points are compacted to those its
 * triangles use, with prismPoint saying which are the prism mesh's (the
 * zone's interface points and the walls' layer points), so that the mesher
 * is not handed the whole prism zone's points.
 * @return 0 on success, 1 with error set otherwise.
 */
/// How BuildJunctionShell relaxes the bands (the layer surfaces' pieces
/// and the zipper strips); see RelaxJunctionBands. Zero rounds leaves the
/// bands as zipped; a zero fraction turns that operation off.
struct BandRelaxation
{
  int rounds = 3;                 // rounds of edge flips and smoothing per level
  double jitter = 0.15;           // the smoothing target moves in the surface by this fraction of the mean edge, hashed direction
  double collapseFraction = 0.3;  // an edge under this fraction of its two triangles' other edges is collapsed (two free points)
  bool removeDegree3 = true;      // a free point on three triangles, one of them under 10 degrees, is taken out
  bool freeThinCaps = true;       // a triangle under 5 degrees keeps neither its normal nor its creases when its point moves
};

/**
 * @brief Widens the junction zone around the given interface points: the
 * junction region they touch (junction triangles connected through shared
 * points) grows by the given rings of triangles. For a zone boundary loop
 * whose layer surface piece came out without a boundary chain of its own
 * (a piece too narrow to keep its shape through the erosion, measured
 * 2026-09-29 on the 178k model at a junction of two vessels 0.1 thick: one
 * row wide at the upper levels with a margin of three rings), so that the
 * next build has a wider piece there.
 * @param seedPoints Interface points whose junction regions grow.
 * @param seedTriangles Interface triangles given to the junction zone first
 * (structured ones whose layer surface passed through a piece), whose
 * regions then grow as well.
 * @param structured Per interface triangle; the widened ones are cleared.
 * @return How many triangles were given to the junction zone.
 */
long long WidenJunctionZone(const Interface &input, const std::vector<long long> &seedPoints,
    const std::vector<long long> &seedTriangles, int rings, std::vector<unsigned char> &structured);

int BuildJunctionShell(const Interface &input, const OffsetField &field, const std::vector<unsigned char> &structured,
    const PrismMesh &prisms, const std::vector<Surface> &levels, const std::vector<double> &fractions,
    int erosionRings, int earPasses, JunctionShell &out, std::string &error, bool innerLevels = true,
    const BandRelaxation &relaxation = BandRelaxation());

/**
 * The relaxation of a junction shell's bands (its layer surfaces' pieces
 * and zipper strips), as BuildJunctionShell runs it: rounds of edge flips
 * and a smoothing of the points onto their levels, with the points flagged
 * fixed (the prism mesh's) kept. On a compacted shell flag the points whose
 * prismPoint is not -1. Fills the shell's numBand* and relaxation* fields.
 * @return 0, or 1 when the flags do not match the points.
 */
int RelaxJunctionBands(JunctionShell &out, const OffsetField &field, const std::vector<double> &fractions,
    const std::vector<unsigned char> &fixed, const BandRelaxation &options);

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
 * @brief Lists the pairs of triangles on one edge that meet at an angle
 * below maxDegrees: folded onto each other. A crossing count does not see
 * them (the two share an edge), and the volume mesher refuses them below its
 * tolerance (foldDegrees) as two facets intersecting - measured 2026-09-23 on
 * the user's model, two outer wall triangles 0.03 degree apart, after every
 * crossing count was zero. The angle is between the two triangles'
 * half-planes on the edge (180 flat out, 0 folded) and does not depend on
 * their winding, and every pair on an edge is measured, so an edge on three
 * triangles - the rim of a layer surface in the wall's shell, shared with
 * the annulus on either side - is judged like any other (found in review
 * 2026-09-23: an edge-count of two let a fold there reach the mesher).
 * Degenerate triangles (the third corner on the edge's line) are left out.
 * @param maxPairs How many pairs to keep, the most folded first; the count
 * returned is of all of them.
 * @param smallestDegrees Set to the smallest angle over every pair on an
 * edge, folded or not (180 when there is none).
 * @return The number of folded pairs.
 */
long long ListFoldedEdges(const std::vector<double> &points, const std::vector<long long> &triangles,
    double maxDegrees, size_t maxPairs, std::vector<FoldedPair> &pairs, double &smallestDegrees);

/**
 * @brief The decimation's crossing guard: the live triangles of a surface in
 * a uniform grid, so that the triangles an operation would make can be
 * tested against those near them (by svenvelope::TrianglesCross) before it
 * is made. Points never move during the decimation, but triangles are
 * re-indexed, so a triangle is removed from its cells before it changes and
 * added after. The triangles an operation replaces are marked as leaving
 * and not tested against; which candidates one test has already seen is
 * kept apart from that, per test, so that several new triangles of one
 * operation are each tested against every triangle near them (with one mark
 * for both, the second test skipped whatever the first had looked at -
 * found in review 2026-09-23).
 */
class CrossingGuard
{
public:
  /// The surface the guard watches; the three are read, never written, and
  /// must outlive the guard. A triangle t is live while dead[t] is zero.
  CrossingGuard(const std::vector<double> &points, const std::vector<long long> &triangles,
      const std::vector<unsigned char> &dead);

  /// Indexes every triangle in a grid of the given cell size over the
  /// points' box (the cell grows until the grid has under 4e8 cells).
  void Build(double cell);

  void Add(long long t);
  void Remove(long long t);

  /// Starts an operation: the triangles marked as leaving are forgotten.
  void BeginOperation();

  /// Marks t as replaced by the current operation: no test is made against it.
  void Leave(long long t);

  /// Whether the triangle T (point ids) passes through a live triangle near
  /// it that is not leaving.
  bool Crosses(const long long T[3]);

private:
  void Box(const long long *T, int lo[3], int hi[3]) const;
  long long Index(int i, int j, int k) const;

  const std::vector<double> &points_;
  const std::vector<long long> &triangles_;
  const std::vector<unsigned char> &dead_;
  double origin_[3] = {0.0, 0.0, 0.0};
  double cell_ = 1.0;
  int n_[3] = {1, 1, 1};
  std::unordered_map<long long, std::vector<long long> > cells_;
  std::vector<long long> leaving_, seen_;   // per triangle: the operation it leaves in, the test that last saw it
  long long operation_ = 0, test_ = 0;
};

}  // namespace svoffset

#endif  // __SV_TETGENMESH_OFFSET_H__
