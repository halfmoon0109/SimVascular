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


/** @file sv_tetgenmesh_envelope.h
 *  @brief The outer envelope of a closed, self-intersecting triangle surface.
 *  @details The wall's outer surface is the inner surface extruded along its
 *  normals, and at a concave junction the extrusion runs through itself: the
 *  sheets extruded from the two vessels cross each other along a crease, and
 *  everything past the crease lies inside the wall. The boundary of the solid
 *  the extrusion encloses is the envelope: the parts of the sheets that have
 *  nothing of the solid on their outer side.
 *
 *  This computes that envelope exactly from the surface's own geometry. Every
 *  pair of triangles that cross is intersected, the crossing segments are
 *  keyed so that the two triangles on an edge and the two on either side of
 *  a crease share the crease points by id, each crossed triangle is split
 *  into pieces along its crease segments, and each piece is kept or dropped
 *  by the winding number of the closed surface on its outer side. Nothing is
 *  moved, nothing is measured against a margin, and there is no hole left to
 *  close: the pieces on the two sides of a crease end on the same points.
 *
 *  It is written against plain arrays so that it can be compiled and run on
 *  its own, without VTK; the caller in sv_tetgenmesh_utils.cxx does the
 *  extrusion and the conversion.
 */

#ifndef __SV_TETGENMESH_ENVELOPE_H__
#define __SV_TETGENMESH_ENVELOPE_H__

#include <string>
#include <vector>

namespace svenvelope {

/**
 * @brief A closed, consistently wound triangle surface.
 * @note The first numSheetTriangles triangles are the sheet the envelope is
 * taken from; the triangles after them close it (the fans over the vessel
 * ends) so that the winding number is defined, and are only ever crossed by a
 * ray - they are never cut and never part of the result.
 *
 * Consistently wound means every edge is traversed once each way by the two
 * triangles on it; the build refuses a surface that is not. Whether the
 * winding puts the geometric normals out of the solid or into it is not
 * assumed: it is read off the surface's signed volume, and the envelope is
 * the same either way, wound the way the input was.
 */
struct Surface
{
  std::vector<double> points;        // three per point
  std::vector<long long> triangles;  // three point ids per triangle
  long long numSheetTriangles = 0;
  // Optional, three per sheet triangle: the direction the triangle faced
  // before it was extruded (the inner triangle's normal). With it, a sheet
  // triangle that has turned over in the extrusion is known exactly; without
  // it, one is guessed from disagreeing with its neighbours.
  std::vector<double> sheetNormal;
  // Optional, three per point: where each point stood before it was
  // extruded - its foot on the interface - for every point the sheet
  // triangles use (the points that close the surface need none). With it,
  // each piece is also put to the prisms the sheet triangles swept out
  // between their feet and themselves: the solid the extrusion stands for
  // is the union of those prisms, and the winding number cannot see that
  // where the extrusion folds, the turned-over part of a prism counts -1
  // and cancels the prism of another sheet that does cover the spot.
  std::vector<double> footPoints;
};

/**
 * @brief The envelope: the pieces of the sheet kept, on the input points
 * followed by the crease points.
 */
struct Envelope
{
  std::vector<double> points;          // the input points first, at their ids, then the crease points
  std::vector<int> pointKind;          // per point: 0 input, 1 on an input edge, 2 a triple point
  std::vector<long long> triangles;    // three point ids per kept piece
  std::vector<long long> source;       // the input sheet triangle each piece came from
  std::vector<unsigned char> whole;    // 1 if the piece is an input triangle kept uncut
  std::vector<unsigned char> creaseEdge; // three per piece: 1 if edge j (from corner j to j+1) lies on a crease, the segment along which two sheets cross
  // The pieces dropped, with the winding number on their outer side, for
  // looking at what was cut away; a piece no ray could classify has winding
  // number -999 here, and a piece dropped as the wall of a pocket -998.
  std::vector<long long> droppedTriangles;
  std::vector<long long> droppedSource;
  std::vector<int> droppedWinding;
};

/**
 * @brief What happened, for the log and for deciding whether the result is
 * one a volume mesher will accept.
 */
struct Report
{
  long long numPairsTested = 0;        // triangle pairs whose boxes overlap
  long long numPairsCrossing = 0;      // pairs that cross along a segment
  long long numCreasePoints = 0;       // points added on edges
  long long numTriplePoints = 0;       // points added where two creases cross inside a triangle
  long long numTrianglesSplit = 0;     // sheet triangles cut into pieces
  long long numPieces = 0;             // pieces of the split triangles
  long long numPiecesKept = 0;
  long long numWholeKept = 0;          // uncut triangles kept
  long long numWholeDropped = 0;       // uncut triangles inside the solid
  long long numInverted = 0;           // sheet triangles wound against their neighbours' side, for the log
  long long numRays = 0;
  long long numRayRetries = 0;         // rays that grazed an edge and were shot again
  long long numUndecided = 0;          // pieces no ray could classify; they are dropped
  long long numPockets = 0;            // components of the kept surface, off the rims, that the rest of the closed surface winds around: the slits of a fold, dropped
  long long numPocketPieces = 0;       // pieces dropped with them
  long long numShreds = 0;             // components of the kept surface, off the rims, that are not closed: the kept wall of a slit whose other wall is turned over, dropped
  long long numShredPieces = 0;        // pieces dropped with them
  long long numPrismFaces = 0;         // faces put to the prisms (only with Surface::footPoints)
  long long numCoveredDropped = 0;     // pieces with winding number zero on their outer side dropped because a prism sweeps over that side: the walls of a fold's slit
  long long numTurnedKept = 0;         // pieces kept the other way round: a prism covers their outer side and nothing covers their inner side
  long long numArrangementFaults = 0;  // a triangle whose pieces do not add up to it, a segment with one end, and the like
  long long numNonManifoldEdges = 0;   // edges of the kept surface on more than two pieces
  long long numMiswoundEdges = 0;      // edges of the kept surface traversed the same way twice
  long long numBoundaryEdges = 0;      // edges of the kept surface on one piece: the cap rims, if all is well
  double seconds = 0.0;
  int windingSense = 1;                // +1 when the input's geometric normals point out of the solid, -1 when into it
  double signedVolume = 0.0;           // of the closed input, by the divergence theorem; its sign is windingSense
  std::string firstFault;              // where and what the first fault was
  double firstFaultAt[3] = {0.0, 0.0, 0.0};
  std::vector<long long> windingHistogram; // pieces per winding number, index = winding + windingOffset
  int windingOffset = 0;
};

/**
 * @brief Builds the outer envelope of a closed surface.
 * @param surface The closed surface; see Surface.
 * @param envelope Set to the envelope.
 * @param report Set to what happened.
 * @param error Set to why, when the build fails outright.
 * @return 0 if the envelope was built (it may still carry faults, counted in
 * the report), 1 if it could not be.
 */
int BuildOuterEnvelope(const Surface &surface, Envelope &envelope, Report &report,
    std::string &error);

/**
 * @brief What the sliver cleanup did.
 */
struct CleanReport
{
  long long numSliversBefore = 0;   // pieces above the aspect limit going in
  long long numSliversAfter = 0;    // and coming out
  double worstBefore = 0.0;         // the largest aspect ratio going in
  double worstAfter = 0.0;          // and coming out
  long long numCollapsed = 0;       // edges collapsed
  long long numFlipped = 0;         // edges flipped
  long long numSnapped = 0;         // apexes put on the long edge, splitting the triangle across it
  long long numRemoved = 0;         // points taken out, their ring triangulated afresh
  long long numNecksCut = 0;        // handles cut: tiny 3-cycles of edges that were no triangle, the two sides closed with a cone each; or a bubble on such a cycle taken off
  long long numUnpinched = 0;       // points made where the surface touched itself at a point, one for each extra fan
  long long numNoMoveAllowed = 0;   // slivers left because no move passed the checks
  long long numLeftCreased = 0;     // of those, the ones touching a crease
  double worstAfterAt[3] = {0.0, 0.0, 0.0};  // the centre of the worst piece coming out
  int numPasses = 0;
  double seconds = 0.0;
};

/**
 * @brief Takes the sliver pieces out of an envelope without moving a crease.
 * @details Cutting a triangle along a crease that passes close to one of its
 * corners, or nearly along one of its edges, leaves a piece with almost no
 * altitude, and a volume mesher has to stand a tetrahedron as flat as that
 * piece on it. The crease is exact and stays where it is; the extruded
 * points around it are free, because the outer surface has no
 * point-for-point relation to the interface. So a sliver is removed by
 * collapsing its shortest edge onto the end that must not move - a point on
 * a boundary loop, a point the caller fixes, or a crease point - or by
 * putting its apex on its long edge and splitting the triangle across that
 * edge at the foot, or, when neither is allowed, by flipping its longest
 * edge, which moves nothing; and when none of those is allowed - a
 * collapse needs the two points' links to meet only at the two apexes,
 * which in the crumple of a fold they often do not - by taking a point of
 * the sliver out altogether and triangulating the ring of its neighbours
 * afresh, on both sides of the crease chord when the point was a crease
 * point with exactly two crease edges.
 *
 * A collapse, a snap or a removal is allowed only when it keeps the surface a
 * manifold, leaves no triangle it touches worse than the worst of those it
 * touched or than the limit (so the local worst never grows), makes no new
 * fold across an edge that is not a crease and deepens none that was there,
 * turns no sound triangle large against the move by more than sixty
 * degrees, and moves the surface by no more than the allowance: for a piece
 * a crease cut, the piece's own altitude (or its short edge, when that is
 * short too), since such a piece and the crumple of a fold around it are
 * artefacts of the cutting and the extrusion; for a sliver of the sheet
 * itself, a twentieth of the mean edge around it, and the point may not
 * slide further than that either, so that no ridge is chamfered. Two
 * triangles no bigger than the move are not judged for folding or turning;
 * the crossing check settles them. A crease point may only be collapsed
 * onto another crease point or put on a crease edge, or off its crease by
 * a tenth of the longest crease edge at it, so a crease is only ever
 * shortened, joined or moved within its own resolution, never bent off its
 * line. A point put on a crease edge becomes a crease point, and the
 * envelope's crease-edge marks follow every move. A flip is allowed only
 * when the two triangles on the edge are nearly coplanar, the edge is not a
 * crease, both new triangles are better than the worse of the old, and the
 * surface moves by no more than a twentieth of an edge. Whether the result
 * crosses itself is not checked here; the caller counts that and decides.
 * @param envelope The envelope, cleaned in place. Pieces it reshapes get
 * whole set to 2; points it drops stay in the point list unreferenced.
 * @param fixedPoint One flag per envelope point, or empty, for points that
 * must not move besides the crease and boundary points.
 * @param aspectLimit The aspect ratio (1 for an equilateral triangle) above
 * which a piece is a sliver.
 * @param report Set to what was done.
 * @return 0 if the envelope was cleaned, 1 if it was malformed (an edge on
 * more than two pieces, ids out of range) and was left as it was.
 */
int CleanEnvelopeSlivers(Envelope &envelope, const std::vector<unsigned char> &fixedPoint,
    double aspectLimit, CleanReport &report);

/**
 * @brief Counts the triangles of a surface that pass through another triangle
 * of it, the way a volume mesher would refuse them.
 * @note Two triangles that share an edge are not tested; two that share a
 * corner cross when the edge opposite the corner of one passes through the
 * other; two that share nothing cross when they meet along a segment.
 * @param points Three per point.
 * @param triangles Three point ids per triangle.
 * @param crossing Set to one per triangle that crosses another.
 * @param firstAt Set to the centre of the first crossing triangle found.
 * @return The number of crossing triangles.
 */
long long CountCrossingTriangles(const std::vector<double> &points,
    const std::vector<long long> &triangles, std::vector<unsigned char> &crossing,
    double firstAt[3]);

/**
 * @brief Two triangles that cross, and the segment along which they do.
 */
struct CrossingPair
{
  long long a = -1;
  long long b = -1;
  double from[3] = {0.0, 0.0, 0.0};
  double to[3] = {0.0, 0.0, 0.0};
};

/**
 * @brief Lists the pairs of triangles that cross, by the judgement of
 * CountCrossingTriangles (a shared corner or an edge touching the other
 * triangle is not a crossing).
 * @param maxPairs How many pairs to keep at most; the count returned is of
 * all of them.
 * @return The number of crossing pairs.
 */
long long ListCrossingPairs(const std::vector<double> &points,
    const std::vector<long long> &triangles, size_t maxPairs,
    std::vector<CrossingPair> &pairs);

/**
 * @brief Whether two triangles cross, by the judgement of
 * CountCrossingTriangles: they meet along a segment of some length that is
 * not a shared edge. A shared corner alone, or an edge of one touching the
 * other at a point, is not a crossing; a degenerate triangle never crosses.
 * @param points Three per point.
 * @param ta Three point ids of the first triangle.
 * @param tb Three point ids of the second.
 */
bool TrianglesCross(const std::vector<double> &points, const long long ta[3], const long long tb[3]);

}

#endif
