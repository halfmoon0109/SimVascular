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
 * @file sv_tetgenmesh_refine.h
 * @brief Edge bisection of a tetrahedral mesh with labelled surface
 * triangles, to refine the wall faces of the fluid mesh to the wall thickness
 * without refining the fluid.
 *
 * Standard library only: no VTK.
 *
 * Why. The solid wall is extruded from the fluid mesh's wall faces, N prism
 * layers of t/N each, and a prism's tetrahedra have a smallest dihedral near
 * atan((t/N) / h) for a face edge h. So the fluid's wall faces had to be
 * meshed at h ~ t, and a thin vessel (t 0.1 to 0.2) got a fluid mesh far
 * finer than the flow needs. The fluid and the solid must share their
 * interface nodes (the FSI solver ties coincident nodes and has no
 * non-matching interface), so the wall faces cannot simply be remeshed for
 * the solid: instead their triangles are split in four, k times, and the
 * split is carried into the fluid only as far as it must go. Bisecting an
 * edge at its midpoint and cutting every cell on that edge in two is a
 * conforming operation (no hanging node anywhere), and in a boundary layer
 * mesh extruded from the wall faces the only cells on a wall edge are the
 * first sublayer's; the interior fluid mesh is untouched.
 *
 * Within one level every edge of the triangles of a face to refine is
 * bisected, the longest first, so each triangle ends in the four-triangle
 * longest-edge partition (the quadrilateral in the middle is cut from the
 * midpoint of the first edge bisected to the vertex opposite it); its
 * angles do not degrade over the levels the way an arbitrary order's can.
 * An edge between a face refined deeper and one refined less is bisected
 * to the deeper face's level, and the shallower face's triangle on it is cut
 * in two along the way (a conforming "green" split).
 */

#ifndef __SV_TETGENMESH_REFINE_H__
#define __SV_TETGENMESH_REFINE_H__

#include <map>
#include <string>
#include <vector>

namespace svrefine
{

/**
 * @brief A tetrahedral mesh with surface (two-dimensional) cells, each
 * carrying integer attributes that children inherit from their parent.
 */
struct Mesh
{
  std::vector<double> points;                 // three per point
  std::vector<long long> tetrahedra;          // four point ids per tetrahedron, positively oriented
  std::vector<std::vector<int> > tetrahedronInts;   // integer attributes, one vector per attribute, one value per tetrahedron
  std::vector<long long> triangles;           // three point ids per surface triangle
  std::vector<std::vector<int> > triangleInts;      // integer attributes, one vector per attribute, one value per triangle
  // Every point added by a bisection, as the two points it is the midpoint
  // of: parents[2*(id - numOriginalPoints)], parents[... + 1]. The caller
  // interpolates any per-point data of its own from these.
  std::vector<long long> parents;
  long long numOriginalPoints = 0;
};

/**
 * @brief An edge, as its two point ids in either order.
 */
struct Edge
{
  long long a, b;
};

/**
 * @brief What happened, for the log.
 */
struct Report
{
  long long numEdgesBisected = 0;
  long long numPointsAdded = 0;
  long long numTetrahedraBefore = 0, numTetrahedraAfter = 0;
  long long numTrianglesBefore = 0, numTrianglesAfter = 0;
  long long numTetrahedraSplit = 0;           // tetrahedra cut in two (a tetrahedron cut twice counts twice)
  std::vector<long long> edgesPerLevel;       // edges bisected at each level, for RefineTriangles
};

/**
 * @brief The depth to which the triangles with each label are refined, from
 * their size against the thickness wanted there.
 *
 * A label's triangles are refined ceil(log2(h / (ratio * t))) times, no
 * less than zero and no more than maxDepth, h being the mean edge of the
 * label's triangles and t the thickness for the label (thicknessOfLabel, or
 * defaultThickness when the label has no entry). Labels in skipLabels (the
 * caps, the layer surfaces, the sidewalls) get depth zero.
 * @param mesh The mesh; its triangles' labels are triangleInts[labelIndex].
 * @param table Set per label to {label, triangles, mean edge, thickness, depth}, for the log.
 * @return The depth per triangle.
 */
struct DepthRow
{
  int label = 0;
  long long numTriangles = 0;
  double meanEdge = 0.0;
  double thickness = 0.0;
  int depth = 0;
};
std::vector<int> DepthOfTriangles(const Mesh &mesh, size_t labelIndex, const std::map<int, double> &thicknessOfLabel,
    double defaultThickness, double ratio, int maxDepth, const std::vector<int> &skipLabels, std::vector<DepthRow> &table);

/**
 * @brief Bisect the given edges at their midpoints, the longest first.
 *
 * Every tetrahedron and every triangle on an edge is cut in two by the
 * midpoint; children keep their parent's orientation and attributes. An
 * edge that is not an edge of any cell is skipped. The edges are bisected
 * in order of decreasing length, as measured before any of them is cut.
 * @return The number of edges bisected.
 */
long long BisectEdges(Mesh &mesh, const std::vector<Edge> &edges, Report &report);

/**
 * @brief Refine the surface triangles to the given depths: at each level
 * every edge of every triangle whose depth is at least that level is
 * bisected (BisectEdges), so a triangle of depth k ends in 4^k and the
 * tetrahedra on its edges in the pieces the bisections leave.
 * @param depth One per triangle; maintained through the refinement
 * (children inherit) and returned at the final triangle count.
 * @return The number of edges bisected in all.
 */
long long RefineTriangles(Mesh &mesh, std::vector<int> &depth, Report &report);

/**
 * @brief Conformity, for tests and for the caller's check: every face of a
 * tetrahedron is on one or two tetrahedra, every surface triangle is a face
 * of exactly one tetrahedron (when the mesh has tetrahedra), no tetrahedron
 * has a non-positive volume.
 * @param error Set to what is wrong.
 * @return true when the mesh is conforming.
 */
bool CheckConforming(const Mesh &mesh, std::string &error);

}  // namespace svrefine

#endif  // __SV_TETGENMESH_REFINE_H__
