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
 * @file sv_tetgenmesh_refine.cxx
 * @brief Edge bisection of a tetrahedral mesh with labelled surface
 * triangles. See the header.
 */

#include "sv_tetgenmesh_refine.h"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <unordered_map>

namespace svrefine
{

namespace
{

typedef long long ll;

double Distance(const double *p, const double *q)
{
  return std::sqrt((p[0] - q[0])*(p[0] - q[0]) + (p[1] - q[1])*(p[1] - q[1]) + (p[2] - q[2])*(p[2] - q[2]));
}

double SignedVolume(const std::vector<double> &pts, const ll *t)
{
  const double *a = &pts[(size_t)3*t[0]], *b = &pts[(size_t)3*t[1]], *c = &pts[(size_t)3*t[2]], *d = &pts[(size_t)3*t[3]];
  double u[3], v[3], w[3];
  for (int k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; w[k] = d[k] - a[k]; }
  return (u[0]*(v[1]*w[2] - v[2]*w[1]) - u[1]*(v[0]*w[2] - v[2]*w[0]) + u[2]*(v[0]*w[1] - v[1]*w[0]))/6.0;
}

// Cells incident to each point, kept up to date as cells are cut: a child
// that takes its parent's slot is already listed at the points it keeps,
// and is added at the midpoint; the other child is a new cell, listed at
// all its points.
struct Incidence
{
  std::vector<std::vector<ll> > ofPoint;
  void Build(ll numPoints, const std::vector<ll> &cells, int size)
  {
    ofPoint.assign((size_t)numPoints, std::vector<ll>());
    for (size_t c = 0; c + (size_t)size - 1 < cells.size(); c += (size_t)size)
    {
      for (int j = 0; j < size; j++) ofPoint[(size_t)cells[c + (size_t)j]].push_back((ll)(c/(size_t)size));
    }
  }
  // the cells on both a and b
  void OnEdge(ll a, ll b, std::vector<ll> &out) const
  {
    out.clear();
    const std::vector<ll> &la = ofPoint[(size_t)a], &lb = ofPoint[(size_t)b];
    for (size_t i = 0; i < la.size(); i++)
    {
      if (std::find(lb.begin(), lb.end(), la[i]) != lb.end()) out.push_back(la[i]);
    }
  }
  void Replace(ll point, ll oldCell, ll newCell)
  {
    std::vector<ll> &l = ofPoint[(size_t)point];
    std::vector<ll>::iterator it = std::find(l.begin(), l.end(), oldCell);
    if (it != l.end()) *it = newCell;
  }
};

}  // namespace

//------------------
// DepthOfTriangles
//------------------
std::vector<int> DepthOfTriangles(const Mesh &mesh, size_t labelIndex, const std::map<int, double> &thicknessOfLabel,
    double defaultThickness, double ratio, int maxDepth, const std::vector<int> &skipLabels, std::vector<DepthRow> &table)
{
  table.clear();
  const ll nt = (ll)(mesh.triangles.size()/3);
  std::vector<int> depth((size_t)nt, 0);
  if (labelIndex >= mesh.triangleInts.size() || (ll)mesh.triangleInts[labelIndex].size() != nt || !(ratio > 0.0)) return depth;
  const std::vector<int> &labels = mesh.triangleInts[labelIndex];
  std::map<int, std::pair<ll, double> > edgeSum;   // label -> (triangles, sum of mean edges)
  for (ll t = 0; t < nt; t++)
  {
    const ll *tri = &mesh.triangles[(size_t)3*t];
    double mean = 0.0;
    for (int j = 0; j < 3; j++) mean += Distance(&mesh.points[(size_t)3*tri[j]], &mesh.points[(size_t)3*tri[(j+1)%3]])/3.0;
    std::pair<ll, double> &s = edgeSum[labels[(size_t)t]];
    s.first++;
    s.second += mean;
  }
  std::map<int, int> depthOfLabel;
  for (std::map<int, std::pair<ll, double> >::const_iterator it = edgeSum.begin(); it != edgeSum.end(); ++it)
  {
    DepthRow row;
    row.label = it->first;
    row.numTriangles = it->second.first;
    row.meanEdge = it->second.second/(double)it->second.first;
    std::map<int, double>::const_iterator th = thicknessOfLabel.find(row.label);
    row.thickness = (th != thicknessOfLabel.end()) ? th->second : defaultThickness;
    row.depth = 0;
    const bool skip = std::find(skipLabels.begin(), skipLabels.end(), row.label) != skipLabels.end();
    if (!skip && row.thickness > 0.0 && row.meanEdge > ratio*row.thickness)
    {
      row.depth = (int)std::ceil(std::log2(row.meanEdge/(ratio*row.thickness)) - 1.0e-9);
      row.depth = std::max(0, std::min(maxDepth, row.depth));
    }
    depthOfLabel[row.label] = row.depth;
    table.push_back(row);
  }
  for (ll t = 0; t < nt; t++) depth[(size_t)t] = depthOfLabel[labels[(size_t)t]];
  return depth;
}

//-------------
// BisectEdges
//-------------
long long BisectEdges(Mesh &mesh, const std::vector<Edge> &edges, Report &report)
{
  const ll numPointsBefore = (ll)(mesh.points.size()/3);
  if (mesh.numOriginalPoints == 0 && mesh.parents.empty()) mesh.numOriginalPoints = numPointsBefore;
  report.numTetrahedraBefore = (ll)(mesh.tetrahedra.size()/4);
  report.numTrianglesBefore = (ll)(mesh.triangles.size()/3);
  // the longest first, as measured now
  std::vector<std::pair<double, size_t> > order;
  order.reserve(edges.size());
  for (size_t e = 0; e < edges.size(); e++)
  {
    const Edge &ed = edges[e];
    if (ed.a < 0 || ed.b < 0 || ed.a >= numPointsBefore || ed.b >= numPointsBefore || ed.a == ed.b) continue;
    order.push_back(std::make_pair(-Distance(&mesh.points[(size_t)3*ed.a], &mesh.points[(size_t)3*ed.b]), e));
  }
  std::sort(order.begin(), order.end());
  Incidence tets, tris;
  tets.Build(numPointsBefore, mesh.tetrahedra, 4);
  tris.Build(numPointsBefore, mesh.triangles, 3);
  std::set<std::pair<ll, ll> > done;
  std::vector<ll> onEdge;
  ll numBisected = 0;
  for (size_t o = 0; o < order.size(); o++)
  {
    const Edge &ed = edges[order[o].second];
    const ll a = std::min(ed.a, ed.b), b = std::max(ed.a, ed.b);
    if (!done.insert(std::make_pair(a, b)).second) continue;
    std::vector<ll> tetsOn, trisOn;
    tets.OnEdge(a, b, tetsOn);
    tris.OnEdge(a, b, trisOn);
    if (tetsOn.empty() && trisOn.empty()) continue;
    // the midpoint
    const ll m = (ll)(mesh.points.size()/3);
    for (int k = 0; k < 3; k++) mesh.points.push_back(0.5*(mesh.points[(size_t)3*a + k] + mesh.points[(size_t)3*b + k]));
    mesh.parents.push_back(a);
    mesh.parents.push_back(b);
    tets.ofPoint.push_back(std::vector<ll>());
    tris.ofPoint.push_back(std::vector<ll>());
    // every tetrahedron on the edge: (a, b, c, d) -> (m, b, c, d) in its
    // slot and (a, m, c, d) new; both keep the orientation, m lying on ab
    for (size_t i = 0; i < tetsOn.size(); i++)
    {
      const ll t = tetsOn[i];
      // by index, not by pointer: the push_backs below may move the array
      ll parent[4], child[4];
      for (int j = 0; j < 4; j++) parent[j] = mesh.tetrahedra[(size_t)4*t + (size_t)j];
      for (int j = 0; j < 4; j++) child[j] = (parent[j] == a) ? m : parent[j];
      const ll c = (ll)(mesh.tetrahedra.size()/4);
      for (int j = 0; j < 4; j++) mesh.tetrahedra.push_back(child[j]);
      for (size_t q = 0; q < mesh.tetrahedronInts.size(); q++) mesh.tetrahedronInts[q].push_back(mesh.tetrahedronInts[q][(size_t)t]);
      for (int j = 0; j < 4; j++) mesh.tetrahedra[(size_t)4*t + (size_t)j] = (parent[j] == b) ? m : parent[j];
      // incidence: the slot t leaves b and joins m; the new cell c joins all its points
      tets.Replace(b, t, -1);
      tets.ofPoint[(size_t)m].push_back(t);
      for (int j = 0; j < 4; j++) tets.ofPoint[(size_t)child[j]].push_back(c);
      report.numTetrahedraSplit++;
    }
    for (size_t i = 0; i < trisOn.size(); i++)
    {
      const ll t = trisOn[i];
      ll parent[3], child[3];
      for (int j = 0; j < 3; j++) parent[j] = mesh.triangles[(size_t)3*t + (size_t)j];
      for (int j = 0; j < 3; j++) child[j] = (parent[j] == a) ? m : parent[j];
      const ll c = (ll)(mesh.triangles.size()/3);
      for (int j = 0; j < 3; j++) mesh.triangles.push_back(child[j]);
      for (size_t q = 0; q < mesh.triangleInts.size(); q++) mesh.triangleInts[q].push_back(mesh.triangleInts[q][(size_t)t]);
      for (int j = 0; j < 3; j++) mesh.triangles[(size_t)3*t + (size_t)j] = (parent[j] == b) ? m : parent[j];
      tris.Replace(b, t, -1);
      tris.ofPoint[(size_t)m].push_back(t);
      for (int j = 0; j < 3; j++) tris.ofPoint[(size_t)child[j]].push_back(c);
    }
    // the -1 slots left at b
    for (int which = 0; which < 2; which++)
    {
      std::vector<ll> &l = (which ? tris : tets).ofPoint[(size_t)b];
      l.erase(std::remove(l.begin(), l.end(), (ll)-1), l.end());
    }
    numBisected++;
  }
  report.numEdgesBisected += numBisected;
  report.numPointsAdded += (ll)(mesh.points.size()/3) - numPointsBefore;
  report.numTetrahedraAfter = (ll)(mesh.tetrahedra.size()/4);
  report.numTrianglesAfter = (ll)(mesh.triangles.size()/3);
  return numBisected;
}

//-----------------
// RefineTriangles
//-----------------
long long RefineTriangles(Mesh &mesh, std::vector<int> &depth, Report &report)
{
  const ll before = (ll)(mesh.tetrahedra.size()/4), trianglesBefore = (ll)(mesh.triangles.size()/3);
  if (depth.size() != (size_t)trianglesBefore) return 0;
  int maxDepth = 0;
  for (size_t t = 0; t < depth.size(); t++) maxDepth = std::max(maxDepth, depth[t]);
  ll total = 0;
  // depth rides along as a triangle attribute so that children inherit it
  mesh.triangleInts.push_back(depth);
  const size_t depthIndex = mesh.triangleInts.size() - 1;
  for (int level = 1; level <= maxDepth; level++)
  {
    std::set<std::pair<ll, ll> > unique;
    std::vector<Edge> edges;
    const std::vector<int> &d = mesh.triangleInts[depthIndex];
    for (size_t t = 0; t + 2 < mesh.triangles.size(); t += 3)
    {
      if (d[t/3] < level) continue;
      for (int j = 0; j < 3; j++)
      {
        const ll a = mesh.triangles[t + (size_t)j], b = mesh.triangles[t + (size_t)((j+1)%3)];
        if (unique.insert(std::make_pair(std::min(a, b), std::max(a, b))).second)
        {
          Edge e; e.a = a; e.b = b;
          edges.push_back(e);
        }
      }
    }
    Report levelReport;
    const ll n = BisectEdges(mesh, edges, levelReport);
    report.edgesPerLevel.push_back(n);
    report.numEdgesBisected += n;
    report.numPointsAdded += levelReport.numPointsAdded;
    report.numTetrahedraSplit += levelReport.numTetrahedraSplit;
    total += n;
  }
  depth = mesh.triangleInts[depthIndex];
  mesh.triangleInts.pop_back();
  report.numTetrahedraBefore = before;
  report.numTrianglesBefore = trianglesBefore;
  report.numTetrahedraAfter = (ll)(mesh.tetrahedra.size()/4);
  report.numTrianglesAfter = (ll)(mesh.triangles.size()/3);
  return total;
}

//-----------------
// CheckConforming
//-----------------
bool CheckConforming(const Mesh &mesh, std::string &error)
{
  const ll np = (ll)(mesh.points.size()/3);
  struct Key
  {
    ll a, b, c;
    bool operator<(const Key &o) const { return a != o.a ? a < o.a : (b != o.b ? b < o.b : c < o.c); }
  };
  auto keyOf = [](ll x, ll y, ll z)
  {
    ll v[3] = {x, y, z};
    std::sort(v, v + 3);
    Key k; k.a = v[0]; k.b = v[1]; k.c = v[2];
    return k;
  };
  std::map<Key, int> faces;
  static const int faceOf[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  for (size_t t = 0; t + 3 < mesh.tetrahedra.size(); t += 4)
  {
    const ll *T = &mesh.tetrahedra[t];
    for (int j = 0; j < 4; j++)
    {
      if (T[j] < 0 || T[j] >= np) { error = "a tetrahedron names a point that is not there"; return false; }
    }
    if (!(SignedVolume(mesh.points, T) > 0.0)) { error = "a tetrahedron of non-positive volume at index " + std::to_string(t/4); return false; }
    for (int f = 0; f < 4; f++) faces[keyOf(T[faceOf[f][0]], T[faceOf[f][1]], T[faceOf[f][2]])]++;
  }
  for (std::map<Key, int>::const_iterator it = faces.begin(); it != faces.end(); ++it)
  {
    if (it->second > 2) { error = "a face on more than two tetrahedra"; return false; }
  }
  for (size_t t = 0; t + 2 < mesh.triangles.size(); t += 3)
  {
    const ll *T = &mesh.triangles[t];
    for (int j = 0; j < 3; j++)
    {
      if (T[j] < 0 || T[j] >= np) { error = "a triangle names a point that is not there"; return false; }
    }
    if (T[0] == T[1] || T[1] == T[2] || T[0] == T[2]) { error = "a triangle naming a point twice"; return false; }
    if (!mesh.tetrahedra.empty())
    {
      std::map<Key, int>::const_iterator it = faces.find(keyOf(T[0], T[1], T[2]));
      if (it == faces.end()) { error = "a surface triangle that is not a face of a tetrahedron, at index " + std::to_string(t/3); return false; }
      if (it->second != 1) { error = "a surface triangle on a face between two tetrahedra, at index " + std::to_string(t/3); return false; }
    }
  }
  for (size_t q = 0; q < mesh.tetrahedronInts.size(); q++)
  {
    if (mesh.tetrahedronInts[q].size() != mesh.tetrahedra.size()/4) { error = "a tetrahedron attribute of the wrong length"; return false; }
  }
  for (size_t q = 0; q < mesh.triangleInts.size(); q++)
  {
    if (mesh.triangleInts[q].size() != mesh.triangles.size()/3) { error = "a triangle attribute of the wrong length"; return false; }
  }
  return true;
}

}  // namespace svrefine
