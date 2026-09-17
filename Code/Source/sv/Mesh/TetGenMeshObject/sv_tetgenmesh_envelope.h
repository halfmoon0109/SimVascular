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
 */
struct Surface
{
  std::vector<double> points;        // three per point
  std::vector<long long> triangles;  // three point ids per triangle
  long long numSheetTriangles = 0;
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
  // The pieces dropped, with the winding number on their outer side, for
  // looking at what was cut away; a piece no ray could classify has winding
  // number -999 here.
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
  long long numArrangementFaults = 0;  // a triangle whose pieces do not add up to it, a segment with one end, and the like
  long long numNonManifoldEdges = 0;   // edges of the kept surface on more than two pieces
  long long numMiswoundEdges = 0;      // edges of the kept surface traversed the same way twice
  long long numBoundaryEdges = 0;      // edges of the kept surface on one piece: the cap rims, if all is well
  double seconds = 0.0;
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

}

#endif
