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

// A standalone test of the edge bisection core (sv_tetgenmesh_refine): no
// VTK, no TetGen. Build and run from the TetGenMeshObject directory with
//
//   g++ -O2 -std=c++17 -I. Testing/test_refine.cxx sv_tetgenmesh_refine.cxx sv_tetgenmesh_offset.cxx
//       sv_tetgenmesh_envelope.cxx -o test_refine && ./test_refine
//
// (the offset core is linked for test 6, the prism layers on a refined tube)
//
// The mesh is a boundary layer as vmtk extrudes it: a wall surface of right
// triangles on a grid (two model faces), two sublayers of prisms below it,
// each prism three tetrahedra, sidewall triangles on the open sides and
// the inner surface at the bottom. The wall faces are refined and the
// checks are that the mesh stays conforming, that the second sublayer is
// untouched, that the pieces inherit their labels, and that the triangles'
// angles hold through the levels.

#include "sv_tetgenmesh_refine.h"
#include "sv_tetgenmesh_offset.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <string>
#include <vector>

typedef long long ll;
static int numFailed = 0;
static void Check(bool ok, const char *what)
{
  printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) numFailed++;
}

static double Volume(const std::vector<double> &p, const ll *t)
{
  const double *a = &p[3*t[0]], *b = &p[3*t[1]], *c = &p[3*t[2]], *d = &p[3*t[3]];
  double u[3], v[3], w[3];
  for (int k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; w[k] = d[k] - a[k]; }
  return (u[0]*(v[1]*w[2] - v[2]*w[1]) - u[1]*(v[0]*w[2] - v[2]*w[0]) + u[2]*(v[0]*w[1] - v[1]*w[0]))/6.0;
}

static double TotalVolume(const svrefine::Mesh &m)
{
  double v = 0.0;
  for (size_t t = 0; t + 3 < m.tetrahedra.size(); t += 4) v += Volume(m.points, &m.tetrahedra[t]);
  return v;
}

// the smallest angle over the triangles with a label, in degrees
static double SmallestAngle(const svrefine::Mesh &m, int label, size_t labelIndex)
{
  double smallest = 180.0;
  for (size_t t = 0; t + 2 < m.triangles.size(); t += 3)
  {
    if (m.triangleInts[labelIndex][t/3] != label) continue;
    for (int j = 0; j < 3; j++)
    {
      const double *a = &m.points[3*m.triangles[t + j]], *b = &m.points[3*m.triangles[t + (j+1)%3]], *c = &m.points[3*m.triangles[t + (j+2)%3]];
      double u[3], v[3], lu = 0.0, lv = 0.0, dot = 0.0;
      for (int k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; lu += u[k]*u[k]; lv += v[k]*v[k]; dot += u[k]*v[k]; }
      smallest = std::min(smallest, std::acos(std::max(-1.0, std::min(1.0, dot/std::sqrt(lu*lv))))*180.0/3.14159265358979323846);
    }
  }
  return smallest;
}

static double LongestEdge(const svrefine::Mesh &m, int label, size_t labelIndex)
{
  double longest = 0.0;
  for (size_t t = 0; t + 2 < m.triangles.size(); t += 3)
  {
    if (m.triangleInts[labelIndex][t/3] != label) continue;
    for (int j = 0; j < 3; j++)
    {
      const double *a = &m.points[3*m.triangles[t + j]], *b = &m.points[3*m.triangles[t + (j+1)%3]];
      longest = std::max(longest, std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2])));
    }
  }
  return longest;
}

static ll CountLabel(const svrefine::Mesh &m, int label, size_t labelIndex)
{
  ll n = 0;
  for (size_t t = 0; t < m.triangleInts[labelIndex].size(); t++) if (m.triangleInts[labelIndex][t] == label) n++;
  return n;
}

// A prism on the triangle (a, b, c) of the layer below (ids + offset) as
// three tetrahedra, the quads' diagonals from the smaller bottom id, so
// that neighbours agree; each tetrahedron turned positive.
static void AddPrism(svrefine::Mesh &m, ll a, ll b, ll c, ll off, int layer)
{
  // rotate so that a is the smallest
  while (!(a < b && a < c)) { ll t = a; a = b; b = c; c = t; }
  ll tets[3][4];
  if (b < c)
  {
    ll T[3][4] = {{a, b, c, c + off}, {a, b, c + off, b + off}, {a, b + off, c + off, a + off}};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 4; j++) tets[i][j] = T[i][j];
  }
  else
  {
    ll T[3][4] = {{a, b, c, b + off}, {a, c, b + off, c + off}, {a, b + off, c + off, a + off}};
    for (int i = 0; i < 3; i++) for (int j = 0; j < 4; j++) tets[i][j] = T[i][j];
  }
  for (int i = 0; i < 3; i++)
  {
    if (Volume(m.points, tets[i]) < 0.0) std::swap(tets[i][1], tets[i][2]);
    for (int j = 0; j < 4; j++) m.tetrahedra.push_back(tets[i][j]);
    m.tetrahedronInts[0].push_back(layer);
  }
}

// The sidewall of a boundary edge (a, b) of the layer below: the quad
// (a, b, b+off, a+off) as the two triangles the prism's diagonal makes.
static void AddSidewall(svrefine::Mesh &m, ll a, ll b, ll off, int label)
{
  const ll lo = std::min(a, b), hi = std::max(a, b);
  // diagonal lo - hi+off
  ll t1[3] = {lo, hi, hi + off}, t2[3] = {lo, hi + off, lo + off};
  for (int j = 0; j < 3; j++) m.triangles.push_back(t1[j]);
  m.triangleInts[0].push_back(label);
  for (int j = 0; j < 3; j++) m.triangles.push_back(t2[j]);
  m.triangleInts[0].push_back(label);
}

// The boundary layer under a grid of n by n cells of spacing 1 (the wall at
// z = 0, faces 10 on the left half and 11 on the right), sublayers at
// depths 0.1 and 0.25.
static svrefine::Mesh BuildLayer(int n)
{
  svrefine::Mesh m;
  m.tetrahedronInts.resize(1);
  m.triangleInts.resize(1);
  const ll perLayer = (ll)(n + 1)*(n + 1);
  const double z[3] = {0.0, -0.1, -0.25};
  for (int layer = 0; layer < 3; layer++)
  {
    for (int j = 0; j <= n; j++) for (int i = 0; i <= n; i++) { m.points.push_back(i); m.points.push_back(j); m.points.push_back(z[layer]); }
  }
  auto id = [&](int i, int j) { return (ll)j*(n + 1) + i; };
  std::vector<std::vector<ll> > wall;   // the wall triangles, with face ids
  std::vector<int> wallFace;
  for (int j = 0; j < n; j++)
  {
    for (int i = 0; i < n; i++)
    {
      const int face = (i < n/2) ? 10 : 11;
      ll t1[3] = {id(i, j), id(i + 1, j), id(i + 1, j + 1)}, t2[3] = {id(i, j), id(i + 1, j + 1), id(i, j + 1)};
      wall.push_back(std::vector<ll>(t1, t1 + 3)); wallFace.push_back(face);
      wall.push_back(std::vector<ll>(t2, t2 + 3)); wallFace.push_back(face);
    }
  }
  // the wall triangles face up (+z): with the layer below them the prisms
  // go down; AddPrism turns each tetrahedron positive
  for (size_t w = 0; w < wall.size(); w++)
  {
    for (int j = 0; j < 3; j++) m.triangles.push_back(wall[w][j]);
    m.triangleInts[0].push_back(wallFace[w]);
    AddPrism(m, wall[w][0], wall[w][1], wall[w][2], perLayer, 1);
    AddPrism(m, wall[w][0] + perLayer, wall[w][1] + perLayer, wall[w][2] + perLayer, perLayer, 2);
    // the inner surface at the bottom, facing down
    for (int j = 2; j >= 0; j--) m.triangles.push_back(wall[w][j] + 2*perLayer);
    m.triangleInts[0].push_back(1);
  }
  // sidewalls on the four open sides, both sublayers
  for (int i = 0; i < n; i++)
  {
    for (int layer = 0; layer < 2; layer++)
    {
      const ll o = layer*perLayer;
      AddSidewall(m, id(i, 0) + o, id(i + 1, 0) + o, perLayer, 9999);
      AddSidewall(m, id(i, n) + o, id(i + 1, n) + o, perLayer, 9999);
      AddSidewall(m, id(0, i) + o, id(0, i + 1) + o, perLayer, 9999);
      AddSidewall(m, id(n, i) + o, id(n, i + 1) + o, perLayer, 9999);
    }
  }
  m.numOriginalPoints = (ll)(m.points.size()/3);
  return m;
}

// the tetrahedra of the second sublayer, as sorted point ids, to see that
// they are untouched
static std::set<std::vector<ll> > SecondLayer(const svrefine::Mesh &m)
{
  std::set<std::vector<ll> > out;
  for (size_t t = 0; t + 3 < m.tetrahedra.size(); t += 4)
  {
    if (m.tetrahedronInts[0][t/4] != 2) continue;
    std::vector<ll> v(&m.tetrahedra[t], &m.tetrahedra[t] + 4);
    std::sort(v.begin(), v.end());
    out.insert(v);
  }
  return out;
}

// The smallest dihedral angle of a tetrahedron, in degrees.
static double MinDihedral(const std::vector<double> &pts, const ll *q)
{
  static const int faces[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  double n[4][3];
  for (int f = 0; f < 4; f++)
  {
    const double *a = &pts[3*q[faces[f][0]]], *b = &pts[3*q[faces[f][1]]], *c = &pts[3*q[faces[f][2]]];
    double u[3], v[3];
    for (int k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; }
    n[f][0] = u[1]*v[2] - u[2]*v[1]; n[f][1] = u[2]*v[0] - u[0]*v[2]; n[f][2] = u[0]*v[1] - u[1]*v[0];
    double l = std::sqrt(n[f][0]*n[f][0] + n[f][1]*n[f][1] + n[f][2]*n[f][2]);
    if (l > 0.0) for (int k = 0; k < 3; k++) n[f][k] /= l;
  }
  double smallest = 180.0;
  for (int f = 0; f < 4; f++)
  {
    for (int g = f + 1; g < 4; g++)
    {
      // outward normals: the dihedral is 180 minus the angle between them
      double d = -(n[f][0]*n[g][0] + n[f][1]*n[g][1] + n[f][2]*n[g][2]);
      smallest = std::min(smallest, std::acos(std::max(-1.0, std::min(1.0, d)))*180.0/3.14159265358979323846);
    }
  }
  return smallest;
}

// A tube of radius r, length L, wall t, meshed with 'around' points around
// and a spacing of about h along; the normals point out of the lumen into
// the wall and the triangles are wound to match.
static svoffset::Interface MakeTube(double r, double L, double t, int around, double h)
{
  svoffset::Interface f;
  const int along = std::max(1, (int)std::ceil(L/h));
  for (int j = 0; j <= along; j++)
  {
    for (int i = 0; i < around; i++)
    {
      const double a = 2.0*3.14159265358979323846*i/around;
      f.points.push_back(r*std::cos(a)); f.points.push_back(r*std::sin(a)); f.points.push_back(L*j/along);
      f.normals.push_back(std::cos(a)); f.normals.push_back(std::sin(a)); f.normals.push_back(0.0);
      f.thickness.push_back(t);
    }
  }
  auto id = [&](int i, int j) { return (ll)j*around + (i % around); };
  for (int j = 0; j < along; j++)
  {
    for (int i = 0; i < around; i++)
    {
      ll t1[3] = {id(i, j), id(i + 1, j), id(i + 1, j + 1)}, t2[3] = {id(i, j), id(i + 1, j + 1), id(i, j + 1)};
      for (int w = 0; w < 2; w++)
      {
        ll *tri = w ? t2 : t1;
        // wound so that the right-hand normal points outward (into the wall)
        const double *a = &f.points[3*tri[0]], *b = &f.points[3*tri[1]], *c = &f.points[3*tri[2]];
        double u[3], v[3], n[3];
        for (int k = 0; k < 3; k++) { u[k] = b[k] - a[k]; v[k] = c[k] - a[k]; }
        n[0] = u[1]*v[2] - u[2]*v[1]; n[1] = u[2]*v[0] - u[0]*v[2]; n[2] = u[0]*v[1] - u[1]*v[0];
        if (n[0]*a[0] + n[1]*a[1] < 0.0) std::swap(tri[1], tri[2]);
        for (int k = 0; k < 3; k++) f.triangles.push_back(tri[k]);
      }
    }
  }
  return f;
}

// The interface refined with the core: the triangles bisected to the
// depth, the midpoints' normals and thickness from their parents, and the
// midpoints put back on the cylinder when asked (as the glue projects them
// onto the model surface).
static svoffset::Interface RefineTube(const svoffset::Interface &in, double thickness, double ratio, bool project, double r, int &levels)
{
  svrefine::Mesh m;
  m.points = in.points;
  m.triangles = in.triangles;
  m.triangleInts.resize(1);
  m.triangleInts[0].assign(in.triangles.size()/3, 1);
  m.numOriginalPoints = (ll)(in.points.size()/3);
  std::map<int, double> th;
  th[1] = thickness;
  std::vector<int> skip;
  std::vector<svrefine::DepthRow> table;
  std::vector<int> depth = svrefine::DepthOfTriangles(m, 0, th, thickness, ratio, 3, skip, table);
  levels = table.empty() ? 0 : table[0].depthMax;
  svrefine::Report rep;
  svrefine::RefineTriangles(m, depth, rep);
  svoffset::Interface out;
  out.points = m.points;
  out.normals = in.normals;
  out.thickness = in.thickness;
  for (size_t p = (size_t)m.numOriginalPoints; p < m.points.size()/3; p++)
  {
    const ll a = m.parents[2*(p - (size_t)m.numOriginalPoints)], b = m.parents[2*(p - (size_t)m.numOriginalPoints) + 1];
    double n[3], l = 0.0;
    for (int k = 0; k < 3; k++) { n[k] = 0.5*(out.normals[3*a + k] + out.normals[3*b + k]); l += n[k]*n[k]; }
    l = std::sqrt(l);
    for (int k = 0; k < 3; k++) out.normals.push_back(l > 0.0 ? n[k]/l : 0.0);
    out.thickness.push_back(0.5*(out.thickness[a] + out.thickness[b]));
    if (project)
    {
      double *x = &out.points[3*p];
      const double rr = std::sqrt(x[0]*x[0] + x[1]*x[1]);
      if (rr > 0.0) { x[0] *= r/rr; x[1] *= r/rr; }
    }
  }
  out.triangles = m.triangles;
  return out;
}

static double PrismSmallestDihedral(const svoffset::Interface &f, int numLayers, ll &numTets, ll &under10)
{
  std::vector<unsigned char> structured(f.triangles.size()/3, 1);
  svoffset::PrismMesh pm;
  std::string err;
  if (svoffset::BuildPrismLayers(f, structured, numLayers, pm, err) != 0) { printf("  prisms failed: %s\n", err.c_str()); numTets = 0; under10 = 0; return 0.0; }
  numTets = (ll)(pm.tetrahedra.size()/4);
  under10 = 0;
  double smallest = 180.0;
  for (ll t = 0; t < numTets; t++)
  {
    const double d = MinDihedral(pm.points, &pm.tetrahedra[4*t]);
    smallest = std::min(smallest, d);
    if (d < 10.0) under10++;
  }
  return smallest;
}

int main()
{
  std::string error;
  {
    printf("test 1: the layer itself, and one edge bisected\n");
    svrefine::Mesh m = BuildLayer(4);
    const ll tetsBefore = (ll)(m.tetrahedra.size()/4), trisBefore = (ll)(m.triangles.size()/3);
    printf("  %lld points, %lld tetrahedra, %lld triangles\n", (ll)(m.points.size()/3), tetsBefore, trisBefore);
    const bool conforming1 = svrefine::CheckConforming(m, error);
    Check(conforming1, ("the layer as built is conforming: " + error).c_str());
    const double volume = TotalVolume(m);
    Check(std::abs(volume - 16.0*0.25) < 1e-12, "its volume is the grid's area times the layer depth");
    std::set<std::vector<ll> > second = SecondLayer(m);
    // the wall edge from (1,1) to (2,1): on two wall triangles and the first sublayer's tetrahedra only
    std::vector<svrefine::Edge> edges(1);
    edges[0].a = 1*5 + 1; edges[0].b = 1*5 + 2;
    svrefine::Report rep;
    ll n = svrefine::BisectEdges(m, edges, rep);
    printf("  bisected %lld: %lld tetrahedra split, %lld points added; %lld tetrahedra, %lld triangles\n", n, rep.numTetrahedraSplit, rep.numPointsAdded, rep.numTetrahedraAfter, rep.numTrianglesAfter);
    Check(n == 1 && rep.numPointsAdded == 1, "one edge, one midpoint");
    const bool conforming2 = svrefine::CheckConforming(m, error);
    Check(conforming2, ("conforming after the bisection: " + error).c_str());
    Check(std::abs(TotalVolume(m) - volume) < 1e-12, "the volume is kept");
    Check(rep.numTrianglesAfter == trisBefore + 2, "the two wall triangles on the edge are cut in two");
    Check(SecondLayer(m) == second, "the second sublayer is untouched");
    Check(m.parents.size() == 2 && std::min(m.parents[0], m.parents[1]) == 6 && std::max(m.parents[0], m.parents[1]) == 7, "the midpoint knows its parents");
    const double *mp = &m.points[3*(m.points.size()/3 - 1)];
    Check(std::abs(mp[0] - 1.5) < 1e-12 && std::abs(mp[1] - 1.0) < 1e-12 && std::abs(mp[2]) < 1e-12, "and lies at the middle of the edge");
    // the same edge again is not bisected twice
    n = svrefine::BisectEdges(m, edges, rep);
    Check(n == 0, "an edge that is no longer an edge (its ends no longer share a cell) is skipped");
  }
  {
    printf("test 2: face 10 refined one level, face 11 not\n");
    svrefine::Mesh m = BuildLayer(4);
    const double volume = TotalVolume(m);
    std::set<std::vector<ll> > second = SecondLayer(m);
    const ll n10 = CountLabel(m, 10, 0), n11 = CountLabel(m, 11, 0), n9999 = CountLabel(m, 9999, 0), n1 = CountLabel(m, 1, 0);
    std::vector<int> depth(m.triangles.size()/3, 0);
    for (size_t t = 0; t < depth.size(); t++) depth[t] = (m.triangleInts[0][t] == 10) ? 1 : 0;
    svrefine::Report rep;
    ll n = svrefine::RefineTriangles(m, depth, rep);
    printf("  %lld edges bisected, %lld tetrahedra -> %lld, %lld triangles -> %lld; face 10: %lld -> %lld, face 11: %lld -> %lld, sidewalls %lld -> %lld, inner %lld -> %lld\n",
        n, rep.numTetrahedraBefore, rep.numTetrahedraAfter, rep.numTrianglesBefore, rep.numTrianglesAfter,
        n10, CountLabel(m, 10, 0), n11, CountLabel(m, 11, 0), n9999, CountLabel(m, 9999, 0), n1, CountLabel(m, 1, 0));
    const bool conforming3 = svrefine::CheckConforming(m, error);
    Check(conforming3, ("conforming: " + error).c_str());
    Check(std::abs(TotalVolume(m) - volume) < 1e-12, "the volume is kept");
    Check(CountLabel(m, 10, 0) == 4*n10, "face 10's triangles are four each");
    // face 11 shares the column of edges at x = 2 with face 10: the four triangles of face 11 on those edges are cut in two
    Check(CountLabel(m, 11, 0) == n11 + 4, "face 11's triangles on the shared edges are cut in two, the others kept");
    Check(CountLabel(m, 1, 0) == n1, "the inner surface is untouched");
    Check(SecondLayer(m) == second, "the second sublayer is untouched");
    // sidewalls: face 10 covers x in [0, 2]: on y = 0 and y = 4 two wall edges each, on x = 0 four; the first sublayer's sidewall triangles on those edges (one per edge, the one with the wall edge) are cut
    Check(CountLabel(m, 9999, 0) == n9999 + 8, "the sidewall triangles on the refined rim edges are cut in two");
    Check(depth.size() == m.triangles.size()/3, "the depths are returned at the new triangle count");
    Check(std::abs(LongestEdge(m, 10, 0) - std::sqrt(2.0)/2.0) < 1e-12, "face 10's longest edge is halved");
    Check(std::abs(SmallestAngle(m, 10, 0) - 45.0) < 1e-9, "right isosceles triangles split into right isosceles triangles (the longest edge first)");
  }
  {
    printf("test 3: both faces two levels\n");
    svrefine::Mesh m = BuildLayer(4);
    const double volume = TotalVolume(m);
    std::set<std::vector<ll> > second = SecondLayer(m);
    const ll n10 = CountLabel(m, 10, 0), n11 = CountLabel(m, 11, 0);
    std::vector<int> depth(m.triangles.size()/3, 0);
    for (size_t t = 0; t < depth.size(); t++) depth[t] = (m.triangleInts[0][t] == 10 || m.triangleInts[0][t] == 11) ? 2 : 0;
    svrefine::Report rep;
    ll n = svrefine::RefineTriangles(m, depth, rep);
    printf("  %lld edges bisected (%lld, %lld per level), %lld tetrahedra -> %lld, %lld triangles -> %lld\n", n, rep.edgesPerLevel[0], rep.edgesPerLevel[1],
        rep.numTetrahedraBefore, rep.numTetrahedraAfter, rep.numTrianglesBefore, rep.numTrianglesAfter);
    const bool conforming4 = svrefine::CheckConforming(m, error);
    Check(conforming4, ("conforming: " + error).c_str());
    Check(std::abs(TotalVolume(m) - volume) < 1e-12, "the volume is kept");
    Check(CountLabel(m, 10, 0) == 16*n10 && CountLabel(m, 11, 0) == 16*n11, "sixteen triangles each");
    Check(SecondLayer(m) == second, "the second sublayer is untouched");
    Check(std::abs(SmallestAngle(m, 10, 0) - 45.0) < 1e-9 && std::abs(SmallestAngle(m, 11, 0) - 45.0) < 1e-9, "the angles hold over two levels");
    Check(std::abs(LongestEdge(m, 11, 0) - std::sqrt(2.0)/4.0) < 1e-12, "the longest edge is quartered");
    // the first sublayer's tetrahedra: 3 per wall triangle before; the wall triangles are 16 each now and every one is a face of a tetrahedron (checked above)
    ll first = 0;
    double smallestFirst = 180.0, smallestSecond = 180.0;
    for (size_t t = 0; t < m.tetrahedronInts[0].size(); t++)
    {
      const double d = MinDihedral(m.points, &m.tetrahedra[4*t]);
      if (m.tetrahedronInts[0][t] == 1) { first++; smallestFirst = std::min(smallestFirst, d); }
      else smallestSecond = std::min(smallestSecond, d);
    }
    printf("  first sublayer: %lld tetrahedra for %lld wall triangles (%.2f each), smallest dihedral %.2f degrees (the untouched second sublayer: %.2f)\n",
        first, CountLabel(m, 10, 0) + CountLabel(m, 11, 0), (double)first/(double)(CountLabel(m, 10, 0) + CountLabel(m, 11, 0)), smallestFirst, smallestSecond);
    // the transition's cost, measured: the wedge's tetrahedron with one wall
    // edge is fanned 2^k ways on its lower edge; at a sublayer a tenth of
    // the edge thick its smallest dihedral is a few degrees to begin with
    Check(smallestFirst > 0.2*smallestSecond, "the first sublayer's smallest dihedral after two levels is above a fifth of the untouched sublayer's (measured 1.9 against 8.5 degrees)");
  }
  {
    printf("test 3b: the depth of the first sublayer's slivers over the levels (the reason for the depth cap)\n");
    for (int k = 1; k <= 3; k++)
    {
      svrefine::Mesh m = BuildLayer(4);
      std::vector<int> depth(m.triangles.size()/3, 0);
      for (size_t t = 0; t < depth.size(); t++) depth[t] = (m.triangleInts[0][t] == 10 || m.triangleInts[0][t] == 11) ? k : 0;
      svrefine::Report rep;
      svrefine::RefineTriangles(m, depth, rep);
      double smallest = 180.0, before = 180.0;
      ll under1 = 0, first = 0;
      for (size_t t = 0; t < m.tetrahedronInts[0].size(); t++)
      {
        const double d = MinDihedral(m.points, &m.tetrahedra[4*t]);
        if (m.tetrahedronInts[0][t] == 1) { first++; smallest = std::min(smallest, d); if (d < 1.0) under1++; }
        else before = std::min(before, d);
      }
      printf("  %d level(s): first sublayer %lld tetrahedra, smallest dihedral %.2f degrees (untouched %.2f), under 1 degree: %lld\n", k, first, smallest, before, under1);
    }
  }
  {
    printf("test 4: a thin triangle's angles over the levels\n");
    // a thin wedge of one prism: the triangle (0,0) (1,0) (0.5,0.1)
    svrefine::Mesh m;
    m.tetrahedronInts.resize(1); m.triangleInts.resize(1);
    double P[6][3] = {{0,0,0}, {1,0,0}, {0.5,0.1,0}, {0,0,-0.1}, {1,0,-0.1}, {0.5,0.1,-0.1}};
    for (int i = 0; i < 6; i++) for (int k = 0; k < 3; k++) m.points.push_back(P[i][k]);
    ll tri[3] = {0, 1, 2};
    for (int j = 0; j < 3; j++) m.triangles.push_back(tri[j]);
    m.triangleInts[0].push_back(10);
    AddPrism(m, 0, 1, 2, 3, 1);
    m.numOriginalPoints = 6;
    const double before = SmallestAngle(m, 10, 0);
    std::vector<int> depth(1, 3);
    svrefine::Report rep;
    svrefine::RefineTriangles(m, depth, rep);
    const double after = SmallestAngle(m, 10, 0);
    printf("  smallest angle %.3f -> %.3f degrees over three levels, %lld triangles\n", before, after, CountLabel(m, 10, 0));
    const bool conforming5 = svrefine::CheckConforming(m, error);
    Check(conforming5, ("conforming: " + error).c_str());
    Check(CountLabel(m, 10, 0) == 64, "sixty-four triangles");
    Check(after >= 0.5*before - 1e-9, "the smallest angle is no less than half the original (the longest-edge partition)");
  }
  {
    printf("test 5: the depth from the thickness\n");
    svrefine::Mesh m = BuildLayer(4);
    std::map<int, double> thickness;
    thickness[10] = 0.3;
    thickness[11] = 1.0;
    std::vector<int> skip;
    skip.push_back(9999);
    skip.push_back(1);
    std::vector<svrefine::DepthRow> table;
    std::vector<int> depth = svrefine::DepthOfTriangles(m, 0, thickness, 0.5, 1.0, 3, skip, table);
    std::map<int, int> byLabel;
    for (size_t r = 0; r < table.size(); r++)
    {
      printf("  label %d: %lld triangles, mean edge %.4f, thickness %.3g, depth %d..%d\n", table[r].label, table[r].numTriangles, table[r].meanEdge, table[r].thickness, table[r].depthMin, table[r].depthMax);
      byLabel[table[r].label] = table[r].depthMax;
      Check(table[r].depthMin == table[r].depthMax, "on a uniform grid every triangle of a label gets the same depth");
    }
    // the mean edge of a right isosceles grid triangle is (1 + 1 + sqrt 2)/3 = 1.138
    Check(byLabel[10] == 2 && byLabel[11] == 1, "1.138 against 0.3 wants two levels, against 1.0 one");
    Check(byLabel[9999] == 0 && byLabel[1] == 0, "the sidewalls and the inner surface are skipped");
    Check(depth.size() == m.triangles.size()/3, "one depth per triangle");
    std::vector<int> none;
    depth = svrefine::DepthOfTriangles(m, 0, thickness, 0.5, 3.0, 3, none, table);
    int most = 0;
    for (size_t t = 0; t < depth.size(); t++) most = std::max(most, depth[t]);
    Check(most == 1, "a ratio of 3 asks one level of face 10 only (1.138 > 0.9) and none of the rest");
  }
  {
    printf("test 7: neighbouring faces three levels apart are balanced to one level\n");
    svrefine::Mesh m = BuildLayer(4);
    std::vector<int> depth(m.triangles.size()/3, 0);
    for (size_t t = 0; t < depth.size(); t++) depth[t] = (m.triangleInts[0][t] == 10) ? 3 : 0;
    std::vector<int> skip;
    skip.push_back(9999);
    skip.push_back(1);
    // without the balancing: face 11's triangles on the shared edges are fanned eight ways
    svrefine::Mesh raw = m;
    std::vector<int> rawDepth = depth;
    svrefine::Report rawRep;
    svrefine::RefineTriangles(raw, rawDepth, rawRep);
    const double rawAngle = SmallestAngle(raw, 11, 0);
    svrefine::Report rep;
    svrefine::RefineTriangles(m, depth, rep, 0, skip);
    const double angle = SmallestAngle(m, 11, 0);
    printf("  face 11's smallest angle: %.2f degrees without the balancing, %.2f with it (%lld triangles raised); face 10 %.2f\n", rawAngle, angle, rep.numTrianglesRaised, SmallestAngle(m, 10, 0));
    const bool conforming7 = svrefine::CheckConforming(m, error);
    Check(conforming7, ("conforming: " + error).c_str());
    Check(rawAngle < 5.0, "unbalanced, the eight-way fan leaves angles under 5 degrees");
    Check(angle > 15.0, "balanced, face 11 keeps its angles above 15 degrees (one green cut per child at most: a right isosceles child cut from a leg's midpoint has 18.4)");
    Check(rep.numTrianglesRaised > 0, "triangles were raised");
    // the raised triangles: within one level of every neighbour
    bool within = true;
    std::map<std::pair<ll, ll>, int> edgeDepth;
    for (size_t t = 0; t + 2 < m.triangles.size(); t += 3)
    {
      if (m.triangleInts[0][t/3] != 10 && m.triangleInts[0][t/3] != 11) continue;
      for (int j = 0; j < 3; j++)
      {
        const ll a = m.triangles[t + j], b = m.triangles[t + (j+1)%3];
        std::pair<ll, ll> key(std::min(a, b), std::max(a, b));
        std::map<std::pair<ll, ll>, int>::iterator it = edgeDepth.find(key);
        if (it == edgeDepth.end()) edgeDepth[key] = depth[t/3];
        else if (std::abs(it->second - depth[t/3]) > 1) within = false;
      }
    }
    Check(within, "after the refinement no two wall triangles on one edge differ by more than one level");
    Check(CountLabel(m, 9999, 0) > 64, "the sidewalls on the rim edges were cut but not balanced (they keep depth 0)");
  }
  {
    printf("test 8: the conformity check rejects what it must\n");
    svrefine::Mesh m = BuildLayer(2);
    const bool ok0 = svrefine::CheckConforming(m, error);
    Check(ok0, "the layer passes");
    // a hanging node: cut one tetrahedron's wall edge without its neighbours
    svrefine::Mesh h = m;
    {
      ll t = 0;
      for (size_t q = 0; q < h.tetrahedronInts[0].size(); q++) if (h.tetrahedronInts[0][q] == 1) { t = (ll)q; break; }
      const ll a = h.tetrahedra[4*t], b = h.tetrahedra[4*t + 1];
      const ll mid = (ll)(h.points.size()/3);
      for (int k = 0; k < 3; k++) h.points.push_back(0.5*(h.points[3*a + k] + h.points[3*b + k]));
      ll child[4] = {mid, h.tetrahedra[4*t + 1], h.tetrahedra[4*t + 2], h.tetrahedra[4*t + 3]};
      h.tetrahedra[4*t + 1] = mid;
      for (int j = 0; j < 4; j++) h.tetrahedra.push_back(child[j]);
      h.tetrahedronInts[0].push_back(1);
    }
    const bool okH = svrefine::CheckConforming(h, error);
    printf("  hanging node: %s\n", error.c_str());
    Check(!okH, "a tetrahedron cut without its neighbours is rejected");
    // a surface triangle that is no face: swap one corner of a wall triangle
    svrefine::Mesh s = m;
    s.triangles[2] = s.triangles[5];
    const bool okS = svrefine::CheckConforming(s, error);
    printf("  stray triangle: %s\n", error.c_str());
    Check(!okS, "a surface triangle that is not a face of a tetrahedron is rejected");
    // an inverted tetrahedron
    svrefine::Mesh v = m;
    std::swap(v.tetrahedra[0], v.tetrahedra[1]);
    const bool okV = svrefine::CheckConforming(v, error);
    Check(!okV, "a tetrahedron of negative volume is rejected");
  }
  {
    printf("test 6: the prism layers of a thin tube meshed coarser than its wall, before and after the refinement\n");
    // r 1, wall 0.1, three layers: the fluid surface at a spacing of 0.3
    // (three times the wall) gives prisms with a smallest dihedral near
    // atan((0.1/3)/0.3) = 6.3 degrees; split to the wall thickness they
    // should come back to the 18 degrees a fluid mesh at the thickness gives.
    const double r = 1.0, t = 0.1;
    svoffset::Interface coarse = MakeTube(r, 3.0, t, 20, 0.3);
    svoffset::Interface fine = MakeTube(r, 3.0, t, 63, 0.1);
    ll nCoarse, u10Coarse, nFine, u10Fine, nRef, u10Ref, nProj, u10Proj;
    const double dCoarse = PrismSmallestDihedral(coarse, 3, nCoarse, u10Coarse);
    const double dFine = PrismSmallestDihedral(fine, 3, nFine, u10Fine);
    int levels = 0;
    svoffset::Interface refined = RefineTube(coarse, t, 1.0, false, r, levels);
    const double dRef = PrismSmallestDihedral(refined, 3, nRef, u10Ref);
    int levels2 = 0;
    svoffset::Interface projected = RefineTube(coarse, t, 1.0, true, r, levels2);
    const double dProj = PrismSmallestDihedral(projected, 3, nProj, u10Proj);
    printf("  coarse (spacing 0.3, %zu triangles): %lld prism tetrahedra, smallest dihedral %.2f degrees, under 10: %lld\n", coarse.triangles.size()/3, nCoarse, dCoarse, u10Coarse);
    printf("  fine (spacing 0.1, %zu triangles): %lld, smallest %.2f, under 10: %lld\n", fine.triangles.size()/3, nFine, dFine, u10Fine);
    printf("  coarse refined %d level(s) at ratio 1 (%zu triangles): %lld, smallest %.2f, under 10: %lld\n", levels, refined.triangles.size()/3, nRef, dRef, u10Ref);
    printf("  the same with the midpoints put back on the cylinder: %lld, smallest %.2f, under 10: %lld\n", nProj, dProj, u10Proj);
    Check(dCoarse < 8.0, "the coarse fluid surface gives flat prisms (under 8 degrees)");
    Check(levels == 2, "a spacing of three thicknesses asks two levels at ratio 1");
    Check(nRef == 16*nCoarse, "the prisms are sixteen per original triangle after two levels");
    Check(dRef > 2.5*dCoarse, "the refined surface's prisms are at least two and a half times better in the smallest dihedral");
    Check(dRef > 0.8*dFine, "and within 20 percent of a fluid surface meshed at the thickness");
    Check(dProj > 0.8*dFine, "with the midpoints on the cylinder as well");
    Check(u10Ref == 0 && u10Proj == 0, "no prism tetrahedron under 10 degrees after the refinement");
  }
  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
