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


// Standalone checks of sv_tetgenmesh_offset on synthetic interfaces. Build
// and run with a plain C++17 compiler and TetGen from the source tree, no
// VTK or SimVascular needed (from this directory):
//
//   T=../../../../../ThirdParty/tetgen/simvascular_tetgen
//   c++ -std=c++17 -O2 -DTETLIBRARY -I.. -I$T test_offset.cxx ../sv_tetgenmesh_offset.cxx
//       ../sv_tetgenmesh_envelope.cxx $T/tetgen.cxx $T/predicates.cxx -o test_offset && ./test_offset
//
// (one command line; TetGen takes a while to compile)
//
// Every case builds an interface open at its vessel ends, takes its offset
// surface, trims it at the cap planes as the glue does, and checks that the
// result is a manifold whose only boundary is the cap rims, that no two of
// its triangles cross, that the wall it makes over the interface is the
// thickness asked for, and that it puts nothing where the wall is solid.

#include "sv_tetgenmesh_offset.h"
#include "sv_tetgenmesh_envelope.h"

#ifndef TETLIBRARY
#define TETLIBRARY   // the build line defines it too, for tetgen.cxx
#endif
#include "tetgen.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <string>
#include <vector>

using svoffset::Interface;
using svoffset::Options;
using svoffset::Surface;
using svoffset::Report;
typedef long long ll;

static int numFailed = 0;

static void Check(bool ok, const char *what)
{
  if (!ok)
  {
    numFailed++;
  }
  printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
}

static double Dot(const double *a, const double *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static void Sub(const double *a, const double *b, double *r) { r[0] = a[0]-b[0]; r[1] = a[1]-b[1]; r[2] = a[2]-b[2]; }
static void Cross(const double *a, const double *b, double *r) { r[0] = a[1]*b[2]-a[2]*b[1]; r[1] = a[2]*b[0]-a[0]*b[2]; r[2] = a[0]*b[1]-a[1]*b[0]; }
static double Norm(const double *a) { return std::sqrt(Dot(a, a)); }
static double Dist(const double *a, const double *b) { double d[3]; Sub(a, b, d); return Norm(d); }

// The Delaunay tetrahedralization by TetGen, points only.
static bool DelaunayWithTetGen(const std::vector<double> &points, std::vector<ll> &tets, void *, std::string &error)
{
  tetgenio in, out;
  in.firstnumber = 0;
  in.numberofpoints = (int)(points.size()/3);
  in.pointlist = new REAL[points.size()];
  for (size_t i = 0; i < points.size(); i++) in.pointlist[i] = points[i];
  char switches[] = "Q";
  try { tetrahedralize(switches, &in, &out); } catch (...) { error = "TetGen threw"; return false; }
  std::map<std::array<ll, 3>, ll> byPos;
  auto keyOf = [](const double *p) { std::array<ll, 3> k = {(ll)std::llround(p[0]*1e6), (ll)std::llround(p[1]*1e6), (ll)std::llround(p[2]*1e6)}; return k; };
  for (size_t i = 0; i < points.size()/3; i++) byPos[keyOf(&points[3*i])] = (ll)i;
  std::vector<ll> toInput(out.numberofpoints, -1);
  for (int i = 0; i < out.numberofpoints; i++)
  {
    std::map<std::array<ll, 3>, ll>::iterator it = byPos.find(keyOf(&out.pointlist[3*i]));
    if (it == byPos.end() || Dist(&out.pointlist[3*i], &points[3*it->second]) > 1e-9) { error = "an output point is not an input point"; return false; }
    toInput[i] = it->second;
  }
  tets.resize((size_t)4*out.numberoftetrahedra);
  for (int t = 0; t < out.numberoftetrahedra; t++) for (int m = 0; m < 4; m++) tets[(size_t)4*t+m] = toInput[out.tetrahedronlist[4*t+m]];
  return true;
}

// A tube along z, open at both ends, normals outward, wound so that the
// right-hand normal points outward.
static void AddTube(Interface &iface, double cx, double cy, double radius, double z0, double z1,
    int numAround, int numAlong, double thickness)
{
  ll base = (ll)(iface.points.size()/3);
  for (int j = 0; j <= numAlong; j++)
  {
    double z = z0 + (z1 - z0)*j/(double)numAlong;
    for (int i = 0; i < numAround; i++)
    {
      double theta = 2.0*M_PI*i/(double)numAround;
      double n[3] = {std::cos(theta), std::sin(theta), 0.0};
      iface.points.push_back(cx + radius*n[0]);
      iface.points.push_back(cy + radius*n[1]);
      iface.points.push_back(z);
      for (int k = 0; k < 3; k++) iface.normals.push_back(n[k]);
      iface.thickness.push_back(thickness);
    }
  }
  auto id = [&](int i, int j) { return base + (ll)j*numAround + (ll)((i + numAround) % numAround); };
  for (int j = 0; j < numAlong; j++)
  {
    for (int i = 0; i < numAround; i++)
    {
      iface.triangles.push_back(id(i, j)); iface.triangles.push_back(id(i+1, j)); iface.triangles.push_back(id(i+1, j+1));
      iface.triangles.push_back(id(i, j)); iface.triangles.push_back(id(i+1, j+1)); iface.triangles.push_back(id(i, j+1));
    }
  }
}

// Trims the surface at the cap planes of the tubes (z = z0 and z = z1), as
// the glue trims at the caps: every triangle whose centre lies beyond a plane
// comes off. The tubes here all share the same two planes.
static void TrimAtPlanes(Surface &surface, double z0, double z1)
{
  std::vector<ll> kept;
  for (size_t i = 0; i + 2 < surface.triangles.size(); i += 3)
  {
    double z = 0.0;
    for (int j = 0; j < 3; j++) z += surface.points[(size_t)3*surface.triangles[i+j] + 2]/3.0;
    if (z < z0 || z > z1) continue;
    kept.insert(kept.end(), &surface.triangles[i], &surface.triangles[i] + 3);
  }
  surface.triangles = kept;
}

struct Measure
{
  ll numBoundary, numNonManifold, numMiswound, numCrossing;
  double minRatio, maxRatio;   // the wall over the interface points away from the ends, against the thickness
  ll numBelow90;
  ll numAbove100;              // triangles of aspect ratio above 100
  ll numTriangles;
};

static double PointTri(const double *p, const double *a, const double *b, const double *c)
{
  double ab[3], ac[3], ap[3]; Sub(b,a,ab); Sub(c,a,ac); Sub(p,a,ap);
  double d1 = Dot(ab,ap), d2 = Dot(ac,ap);
  if (d1 <= 0 && d2 <= 0) return Dist(p,a);
  double bp[3]; Sub(p,b,bp); double d3 = Dot(ab,bp), d4 = Dot(ac,bp);
  if (d3 >= 0 && d4 <= d3) return Dist(p,b);
  double vc = d1*d4 - d3*d2;
  if (vc <= 0 && d1 >= 0 && d3 <= 0) { double v = (d1-d3 != 0) ? d1/(d1-d3) : 0; double q[3] = {a[0]+v*ab[0], a[1]+v*ab[1], a[2]+v*ab[2]}; return Dist(p,q); }
  double cp[3]; Sub(p,c,cp); double d5 = Dot(ab,cp), d6 = Dot(ac,cp);
  if (d6 >= 0 && d5 <= d6) return Dist(p,c);
  double vb = d5*d2 - d1*d6;
  if (vb <= 0 && d2 >= 0 && d6 <= 0) { double w = (d2-d6 != 0) ? d2/(d2-d6) : 0; double q[3] = {a[0]+w*ac[0], a[1]+w*ac[1], a[2]+w*ac[2]}; return Dist(p,q); }
  double va = d3*d6 - d5*d4;
  if (va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0) { double w = (d4-d3)/((d4-d3)+(d5-d6)); double q[3] = {b[0]+w*(c[0]-b[0]), b[1]+w*(c[1]-b[1]), b[2]+w*(c[2]-b[2])}; return Dist(p,q); }
  double den = va+vb+vc; if (den == 0) return Dist(p,a);
  double v = vb/den, w = vc/den; double q[3] = {a[0]+ab[0]*v+ac[0]*w, a[1]+ab[1]*v+ac[1]*w, a[2]+ab[2]*v+ac[2]*w}; return Dist(p,q);
}

// The distance from p to the surface, by brute force (the tests are small).
static double DistanceToSurface(const Surface &s, const double *p)
{
  double best = 1e300;
  for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
  {
    best = std::min(best, PointTri(p, &s.points[(size_t)3*s.triangles[i]], &s.points[(size_t)3*s.triangles[i+1]], &s.points[(size_t)3*s.triangles[i+2]]));
  }
  return best;
}

static Measure MeasureSurface(const Interface &iface, const Surface &s, double zMin, double zMax, double endMargin)
{
  Measure m;
  svoffset::CountEdges(s.triangles, m.numBoundary, m.numNonManifold, m.numMiswound);
  std::vector<unsigned char> crossing;
  double at[3];
  m.numCrossing = svenvelope::CountCrossingTriangles(s.points, s.triangles, crossing, at);
  m.minRatio = 1e300; m.maxRatio = 0.0; m.numBelow90 = 0;
  for (size_t i = 0; i < iface.points.size()/3; i++)
  {
    const double *p = &iface.points[3*i];
    if (p[2] < zMin + endMargin || p[2] > zMax - endMargin) continue;
    double ratio = DistanceToSurface(s, p)/iface.thickness[i];
    m.minRatio = std::min(m.minRatio, ratio);
    m.maxRatio = std::max(m.maxRatio, ratio);
    if (ratio < 0.9) m.numBelow90++;
  }
  m.numAbove100 = 0;
  m.numTriangles = (ll)(s.triangles.size()/3);
  for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
  {
    const double *a = &s.points[(size_t)3*s.triangles[i]], *b = &s.points[(size_t)3*s.triangles[i+1]], *c = &s.points[(size_t)3*s.triangles[i+2]];
    double la = Dist(b,c), lb = Dist(c,a), lc = Dist(a,b), per = la+lb+lc; double e1[3], e2[3], n[3]; Sub(b,a,e1); Sub(c,a,e2); Cross(e1,e2,n);
    double inr = (per > 0) ? Norm(n)/per : 0; double longest = std::max(la, std::max(lb, lc));
    if (!(inr > 0) || longest/(2*inr) > 100.0) m.numAbove100++;
  }
  return m;
}

static bool Build(const Interface &iface, const Options &options, Surface &s, Report &r)
{
  std::string error;
  int rc = svoffset::BuildOffsetSurface(iface, options, DelaunayWithTetGen, nullptr, s, r, error);
  if (rc != 0)
  {
    printf("  FAIL build: %s\n", error.c_str());
    numFailed++;
    return false;
  }
  printf("  built: %lld cloud points, %lld tetrahedra (%lld cut), contour %lld -> %lld triangles after %lld collapses; boundary %lld, non-manifold %lld, miswound %lld; %.1f s\n",
      r.numCloudPoints, r.numTetrahedra, r.numTetrahedraCut, r.numContourTriangles, r.numTriangles, r.numCollapsed,
      r.numBoundaryEdges, r.numNonManifoldEdges, r.numMiswoundEdges, r.secondsField + r.secondsDelaunay + r.secondsContour + r.secondsDecimate);
  return true;
}

// 1. A tube: the offset is a tube a wall wider.
static void TestTube()
{
  printf("test 1: a tube of radius 1 with a wall of 0.3\n");
  Interface iface;
  AddTube(iface, 0.0, 0.0, 1.0, 0.0, 10.0, 32, 40, 0.3);
  Options options;
  Surface s; Report r;
  if (!Build(iface, options, s, r)) return;
  Check(r.numRims == 2, "two cap rims found");
  Check(r.numNonManifoldEdges == 0 && r.numMiswoundEdges == 0, "the offset surface is a manifold, consistently wound");
  // every point past a cap plane belongs to that cap, and none well inside the tube belongs to any
  {
    ll numUnowned = 0, numWrong = 0, numOwnedInside = 0;
    for (size_t i = 0; i < s.points.size()/3; i++)
    {
      double z = s.points[3*i + 2];
      ll rim = s.pointRim[i];
      if (z > 10.0 || z < 0.0)
      {
        if (rim < 0) numUnowned++;
        else
        {
          double rimZ = iface.points[3*(size_t)s.rims[(size_t)rim][0] + 2];
          if ((z > 10.0 && rimZ < 5.0) || (z < 0.0 && rimZ > 5.0)) numWrong++;
        }
      }
      else if (z > 2.0 && z < 8.0 && rim >= 0) numOwnedInside++;
    }
    Check(numUnowned == 0 && numWrong == 0, "every point past a cap plane belongs to that cap");
    Check(numOwnedInside == 0, "no point well inside the tube belongs to a cap");
  }
  TrimAtPlanes(s, 0.0, 10.0);
  Measure m = MeasureSurface(iface, s, 0.0, 10.0, 0.6);
  printf("  after the trim: boundary %lld, non-manifold %lld, miswound %lld, crossing %lld; wall over the interface %.3f..%.3f of the thickness, %lld points below 0.9; %lld of %lld triangles above aspect 100\n",
      m.numBoundary, m.numNonManifold, m.numMiswound, m.numCrossing, m.minRatio, m.maxRatio, m.numBelow90, m.numAbove100, m.numTriangles);
  Check(m.numBoundary > 0 && m.numNonManifold == 0 && m.numMiswound == 0, "trimmed: manifold, open only at the rims");
  Check(m.numCrossing == 0, "no two triangles cross");
  Check(m.numBelow90 == 0, "the wall over the interface is at least 0.9 of the thickness away from the ends");
  Check(m.maxRatio < 1.15, "and no more than 1.15 of it");
  // every surface point away from the ends stands between r + 0.9t and r + 1.1t off the axis
  ll numOff = 0;
  for (size_t i = 0; i < s.points.size()/3; i++)
  {
    const double *p = &s.points[3*i];
    if (p[2] < 0.6 || p[2] > 9.4) continue;
    double rho = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    if (rho < 1.0 + 0.9*0.3 || rho > 1.0 + 1.1*0.3) numOff++;
  }
  Check(numOff == 0, "every surface point lies a wall off the tube, none in the lumen");
  Check(m.numAbove100*100 < m.numTriangles, "fewer than 1% of the triangles have an aspect ratio above 100");
  Check(r.numTriangles < 6*(ll)(iface.triangles.size()/3), "the decimated surface has fewer than six triangles per interface triangle");
}

// 2. A thin coarse tube: the wall is thinner than the interface's edge.
static void TestThinTube()
{
  printf("test 2: a tube of radius 0.3 with a wall of 0.1, edges of 0.15\n");
  Interface iface;
  AddTube(iface, 0.0, 0.0, 0.3, 0.0, 6.0, 12, 40, 0.1);
  Options options;
  Surface s; Report r;
  if (!Build(iface, options, s, r)) return;
  TrimAtPlanes(s, 0.0, 6.0);
  Measure m = MeasureSurface(iface, s, 0.0, 6.0, 0.3);
  printf("  after the trim: boundary %lld, non-manifold %lld, miswound %lld, crossing %lld; wall %.3f..%.3f, %lld below 0.9; %lld of %lld above aspect 100\n",
      m.numBoundary, m.numNonManifold, m.numMiswound, m.numCrossing, m.minRatio, m.maxRatio, m.numBelow90, m.numAbove100, m.numTriangles);
  Check(m.numBoundary > 0 && m.numNonManifold == 0 && m.numMiswound == 0, "trimmed: manifold, open only at the rims");
  Check(m.numCrossing == 0, "no two triangles cross");
  Check(m.numBelow90 == 0, "the wall over the interface is at least 0.9 of the thickness");
  Check(m.maxRatio < 1.2, "and no more than 1.2 of it");
  ll numOff = 0;
  for (size_t i = 0; i < s.points.size()/3; i++)
  {
    const double *p = &s.points[3*i];
    if (p[2] < 0.3 || p[2] > 5.7) continue;
    double rho = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    if (rho < 0.3 + 0.85*0.1 || rho > 0.3 + 1.15*0.1) numOff++;
  }
  Check(numOff == 0, "every surface point lies a wall off the tube, none in the lumen");
}

// 3. Two tubes side by side, closer than their walls: the septum between them
// is solid and the surface creases where the two offsets meet.
static void TestSeptum()
{
  printf("test 3: two tubes whose walls overlap between them (radius 1 wall 0.3, radius 0.5 wall 0.1, interfaces 0.3 apart)\n");
  Interface iface;
  AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3);
  AddTube(iface, 1.8, 0.0, 0.5, 0.0, 8.0, 20, 32, 0.1);
  Options options;
  Surface s; Report r;
  if (!Build(iface, options, s, r)) return;
  Check(r.numRims == 4, "four cap rims found");
  TrimAtPlanes(s, 0.0, 8.0);
  Measure m = MeasureSurface(iface, s, 0.0, 8.0, 0.6);
  printf("  after the trim: boundary %lld, non-manifold %lld, miswound %lld, crossing %lld; wall %.3f..%.3f, %lld below 0.9; %lld of %lld above aspect 100\n",
      m.numBoundary, m.numNonManifold, m.numMiswound, m.numCrossing, m.minRatio, m.maxRatio, m.numBelow90, m.numAbove100, m.numTriangles);
  Check(m.numBoundary > 0 && m.numNonManifold == 0 && m.numMiswound == 0, "trimmed: manifold, open only at the rims");
  Check(m.numCrossing == 0, "no two triangles cross");
  Check(m.numBelow90 == 0, "the wall over every interface point is at least 0.9 of its thickness");
  // the septum: points midway between the two interfaces are deep inside
  double worst = 1e300;
  for (int j = 2; j < 30; j++)
  {
    double p[3] = {1.15, 0.0, 8.0*j/32.0};
    worst = std::min(worst, DistanceToSurface(s, p));
  }
  printf("  the surface stays at least %.3f from the middle of the septum (the crease is 0.28 off it)\n", worst);
  Check(worst > 0.2, "nothing of the surface runs through the septum between the tubes");
  // and the surface reaches around both tubes: points on the far sides
  ll numOff = 0;
  for (size_t i = 0; i < s.points.size()/3; i++)
  {
    const double *p = &s.points[3*i];
    if (p[2] < 0.6 || p[2] > 7.4) continue;
    double rho1 = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    double rho2 = std::sqrt((p[0]-1.8)*(p[0]-1.8) + p[1]*p[1]);
    bool onOffset1 = std::abs(rho1 - 1.3) < 0.05, onOffset2 = std::abs(rho2 - 0.6) < 0.03;
    if (!onOffset1 && !onOffset2) numOff++;
  }
  Check(numOff == 0, "every surface point lies on one of the two offsets");
}

// The number of boundary loops of a surface.

// Trims as the glue does with the ownership: every triangle past a cap's
// plane whose corners belong to that cap comes off; triangles of a
// neighbouring vessel that run past the plane stay. The cap planes here are
// z = z0 and z = z1 of each tube, told apart by the rim's own z.
// The cap planes of an interface as the glue computes them: each rim's
// centroid and Newell normal, the outward direction being the side the
// interface is not on.
static std::vector<svoffset::CapPlane> PlanesFromRims(const Interface &iface, const Surface &s)
{
  std::vector<svoffset::CapPlane> planes;
  for (size_t r = 0; r < s.rims.size(); r++)
  {
    const std::vector<ll> &loop = s.rims[r];
    svoffset::CapPlane plane;
    plane.rim = (ll)r;
    for (int k = 0; k < 3; k++) plane.origin[k] = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
      for (int k = 0; k < 3; k++) plane.origin[k] += iface.points[3*(size_t)loop[m] + k]/(double)loop.size();
    double normal[3] = {0.0, 0.0, 0.0}, radius = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
    {
      const double *a = &iface.points[3*(size_t)loop[m]], *b = &iface.points[3*(size_t)loop[(m+1)%loop.size()]];
      normal[0] += (a[1]-b[1])*(a[2]+b[2]);
      normal[1] += (a[2]-b[2])*(a[0]+b[0]);
      normal[2] += (a[0]-b[0])*(a[1]+b[1]);
      radius = std::max(radius, Dist(a, plane.origin));
    }
    double L = Norm(normal);
    for (int k = 0; k < 3; k++) plane.outward[k] = -normal[k]/L;
    double plus = 0.0, minus = 0.0;
    for (size_t i = 0; i < iface.points.size()/3; i++)
    {
      double off[3]; Sub(&iface.points[3*i], plane.origin, off);
      if (Norm(off) > 3.0*radius) continue;
      double h = Dot(off, plane.outward);
      if (h > 0.05*radius) plus += 1.0; else if (h < -0.05*radius) minus += 1.0;
    }
    if (plus > minus) for (int k = 0; k < 3; k++) plane.outward[k] = -plane.outward[k];
    planes.push_back(plane);
  }
  return planes;
}

// The boundary loops of a surface, and for each the plane (index) all of its
// points lie on within tol, -1 for none.
static std::vector<int> LoopPlanes(const Surface &s, const std::vector<svoffset::CapPlane> &planes, double tol, int &numLoops)
{
  std::map<std::pair<ll,ll>, int> count;
  for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
    for (int j = 0; j < 3; j++) { ll a = s.triangles[i+j], b = s.triangles[i+(j+1)%3]; count[std::make_pair(std::min(a,b), std::max(a,b))]++; }
  std::map<ll,ll> next;
  for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
    for (int j = 0; j < 3; j++) { ll a = s.triangles[i+j], b = s.triangles[i+(j+1)%3]; if (count[std::make_pair(std::min(a,b), std::max(a,b))] == 1) next[a] = b; }
  std::set<ll> seen; std::vector<int> onPlane; numLoops = 0;
  for (std::map<ll,ll>::iterator it = next.begin(); it != next.end(); ++it)
  {
    if (seen.count(it->first)) continue;
    std::vector<ll> loop; ll cur = it->first;
    while (!seen.count(cur)) { seen.insert(cur); loop.push_back(cur); std::map<ll,ll>::iterator nx = next.find(cur); if (nx == next.end()) break; cur = nx->second; }
    numLoops++;
    // the plane it lies on; two caps can share a plane (both tubes' low ends
    // here), so of those the one whose origin is nearest the loop's centre
    double centre[3] = {0.0, 0.0, 0.0};
    for (size_t m = 0; m < loop.size(); m++) for (int k = 0; k < 3; k++) centre[k] += s.points[3*(size_t)loop[m] + k]/(double)loop.size();
    int which = -1; double nearest = 1e300;
    for (size_t p = 0; p < planes.size(); p++)
    {
      double worst = 0.0;
      for (size_t m = 0; m < loop.size(); m++) { double off[3]; Sub(&s.points[3*(size_t)loop[m]], planes[p].origin, off); worst = std::max(worst, std::abs(Dot(off, planes[p].outward))); }
      if (worst > tol) continue;
      double d = Dist(centre, planes[p].origin);
      if (d < nearest) { nearest = d; which = (int)p; }
    }
    onPlane.push_back(which);
  }
  return onPlane;
}

// 4. Two tubes side by side ending at different heights, close enough that a
// window around the shorter one's cap would reach the longer one's wall: the
// trim must take only what belongs to each cap.
static void TestNeighbouringEnds()
{
  // The axes are 2.2 apart: the walls stay 0.3 clear of each other, yet
  // the long tube's interface runs within a collar's length of the short
  // tube's rim, which is where a cap's reach must be decided along the
  // interface and not by a cylinder about the cap (the user's model has two
  // such vessels at (0.4, -54.4, 104.6), and the cylinder handed the
  // neighbour's wall to the cap and cut a hole in it).
  const double dx = 2.2;
  printf("test 4: a tube ending at z = 8 beside one running on to z = 12, axes %.1f apart (radius 1 wall 0.3, radius 0.5 wall 0.1)\n", dx);
  Interface iface;
  AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3);
  AddTube(iface, dx, 0.0, 0.5, 0.0, 12.0, 20, 48, 0.1);
  Options options;
  Surface s; Report r;
  if (!Build(iface, options, s, r)) return;
  Check(r.numRims == 4 && s.rims.size() == 4, "four cap rims found and reported");
  // ownership past the short tube's end
  ll numUnownedPast = 0;
  for (size_t i = 0; i < s.points.size()/3; i++)
  {
    const double *p = &s.points[3*i];
    double rho1 = std::sqrt(p[0]*p[0] + p[1]*p[1]);
    if (p[2] > 8.0 && rho1 < 1.6 && s.pointRim[i] < 0) numUnownedPast++;
  }
  Check(numUnownedPast == 0, "every point of the short tube's dome past its cap belongs to a cap");
  // and nothing of the long tube's wall is handed to the short tube's cap
  ll topRim = -1;
  for (size_t r = 0; r < s.rims.size(); r++)
  {
    if (std::abs(iface.points[3*(size_t)s.rims[r][0] + 2] - 8.0) < 1e-6) topRim = (ll)r;
  }
  Check(topRim >= 0, "the short tube's cap rim at z = 8 is among the rims");
  ll numStolen = 0;
  for (size_t i = 0; i < s.points.size()/3; i++)
  {
    const double *p = &s.points[3*i];
    double rho2 = std::sqrt((p[0]-dx)*(p[0]-dx) + p[1]*p[1]);
    if (rho2 < 0.65 && p[2] > 8.0 + 1e-6 && s.pointRim[i] == topRim) numStolen++;
  }
  printf("  points of the long tube's wall past z = 8 owned by the short tube's cap: %lld\n", numStolen);
  Check(numStolen == 0, "the neighbouring tube's wall belongs to no cap of the short tube");
  // The trim itself, as the glue calls it. Without the snap band the cut
  // lands a hair from the contour points that lie exactly in the cap planes
  // (the points offset along the rim's own normals do), leaving rim edges of
  // nothing and triangles flat in the plane; with it, nothing shorter than a
  // fraction of an edge and nothing leaning over the annulus.
  std::vector<svoffset::CapPlane> planes = PlanesFromRims(iface, s);
  {
    Surface bare = s;
    svoffset::TrimReport bareReport; std::string bareError;
    if (svoffset::TrimSurfaceAtCaps(bare, planes, 0.0, bareReport, bareError) == 0)
      printf("  trimmed without the snap band: shortest rim edge %.2e of the local edge, smallest rim angle %.4f degrees, %lld triangles split, %lld dropped\n",
          bareReport.shortestRimEdgeRatio, bareReport.smallestRimAngleDegrees, bareReport.numSplit, bareReport.numDropped);
  }
  svoffset::TrimReport trim; std::string trimError;
  int trimRc = svoffset::TrimSurfaceAtCaps(s, planes, 0.1, trim, trimError);
  Check(trimRc == 0, trimRc == 0 ? "the surface was trimmed at its four cap planes" : trimError.c_str());
  if (trimRc != 0) return;
  ll snapped = 0, cut = 0;
  for (size_t p = 0; p < planes.size(); p++) { snapped += trim.numSnapped[p]; cut += trim.numCut[p]; }
  printf("  trimmed: %lld points cut, %lld snapped onto a plane, %lld triangles removed, %lld split, %lld dropped (an edge or a point on the kept side, or flat in the plane); %lld -> %lld points, %lld -> %lld triangles; %lld rim edges, shortest %.3f of the local edge, smallest rim angle %.2f degrees\n",
      cut, snapped, trim.numRemoved, trim.numSplit, trim.numDropped, trim.numPointsBefore, trim.numPointsAfter, trim.numTrianglesBefore, trim.numTrianglesAfter,
      trim.numRimEdges, trim.shortestRimEdgeRatio, trim.smallestRimAngleDegrees);
  int loops = 0;
  std::vector<int> loopPlane = LoopPlanes(s, planes, 1e-6, loops);
  std::set<int> planesHit; int stray = 0;
  for (size_t l = 0; l < loopPlane.size(); l++) { if (loopPlane[l] < 0) stray++; else planesHit.insert(loopPlane[l]); }
  ll nb, nn, nm;
  svoffset::CountEdges(s.triangles, nb, nn, nm);
  std::vector<unsigned char> crossing;
  double at[3];
  ll numCrossing = svenvelope::CountCrossingTriangles(s.points, s.triangles, crossing, at);
  printf("  after the trim: %d boundary loops (%zu planes each with one, %d on no plane), non-manifold %lld, miswound %lld, crossing %lld\n", loops, planesHit.size(), stray, nn, nm, numCrossing);
  Check(loops == 4 && planesHit.size() == 4 && stray == 0, "the trim leaves exactly the four cap rims, one on each plane");
  Check(nn == 0 && nm == 0 && numCrossing == 0, "and a manifold without crossings");
  Check(trim.shortestRimEdgeRatio >= 0.05, "no rim edge is shorter than a twentieth of the local edge");
  Check(trim.smallestRimAngleDegrees >= 1.0, "no rim triangle leans over the annulus within a degree");
  // nothing TetGen would merge: no two points of the surface within 1e-4 of an edge
  {
    double shortest = 1e300;
    for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
      for (int j = 0; j < 3; j++) shortest = std::min(shortest, Dist(&s.points[3*(size_t)s.triangles[i+j]], &s.points[3*(size_t)s.triangles[i+(j+1)%3]]));
    printf("  shortest edge anywhere after the trim: %.4g\n", shortest);
    Check(shortest > 1e-3, "no edge of the trimmed surface is shorter than 1e-3");
  }
  // the long tube's wall between z = 8.5 and 11.5 is untouched
  double worst = 1e300;
  for (size_t i = 0; i < iface.points.size()/3; i++)
  {
    const double *p = &iface.points[3*i];
    if (std::abs(std::sqrt((p[0]-dx)*(p[0]-dx) + p[1]*p[1]) - 0.5) > 1e-6) continue;
    if (p[2] < 8.5 || p[2] > 11.5) continue;
    worst = std::min(worst, DistanceToSurface(s, p)/iface.thickness[i]);
  }
  printf("  the long tube's wall past the short tube's end: at least %.3f of the thickness\n", worst);
  Check(worst > 0.9, "the neighbouring tube's wall past the short tube's cap is intact");
}

// 5. Layers: offsets at 1/3, 2/3 and 1 of the thickness of a tube, each
// trimmed at the cap planes, must nest without crossing one another or the
// interface, each a layer thick from the interface, at one size (the chord
// tolerance scaled by the reciprocal of the fraction).
static void TestLayers()
{
  printf("test 5: three layer surfaces of a tube of radius 1 with a wall of 0.3\n");
  Interface iface;
  AddTube(iface, 0.0, 0.0, 1.0, 0.0, 10.0, 32, 40, 0.3);
  const int numLayers = 3;
  std::vector<Surface> levels;
  std::vector<double> allPts;
  std::vector<ll> allTris;
  // the interface itself, for the crossing check
  allPts = iface.points;
  allTris = iface.triangles;
  Surface interfaceAsSurface;
  interfaceAsSurface.points = iface.points;
  interfaceAsSurface.triangles = iface.triangles;
  bool ok = true;
  for (int k = 1; k <= numLayers && ok; k++)
  {
    double fraction = (double)k/numLayers;
    Interface scaled = iface;
    for (size_t i = 0; i < scaled.thickness.size(); i++) scaled.thickness[i] *= fraction;
    Options options;
    options.chordTolerance = 0.05/fraction;
    Surface s; Report r;
    if (!Build(scaled, options, s, r)) { ok = false; break; }
    std::vector<svoffset::CapPlane> planes = PlanesFromRims(scaled, s);
    svoffset::TrimReport trim; std::string err;
    if (svoffset::TrimSurfaceAtCaps(s, planes, 0.1, trim, err) != 0) { printf("  FAIL trim of layer %d: %s\n", k, err.c_str()); numFailed++; ok = false; break; }
    // the layer is k/N of the wall from the interface, away from the ends
    double worst = 1e300, best = 0.0;
    for (size_t i = 0; i < s.points.size()/3; i++)
    {
      const double *p = &s.points[3*i];
      if (p[2] < 0.6 || p[2] > 9.4) continue;
      double d = DistanceToSurface(interfaceAsSurface, p)/(0.3*fraction);
      worst = std::min(worst, d); best = std::max(best, d);
    }
    // its size: the mean edge against the interface's
    double meanEdge = 0.0; ll numEdges = 0;
    for (size_t i = 0; i + 2 < s.triangles.size(); i += 3) for (int j = 0; j < 3; j++) { meanEdge += Dist(&s.points[3*(size_t)s.triangles[i+j]], &s.points[3*(size_t)s.triangles[i+(j+1)%3]]); numEdges++; }
    meanEdge /= (double)numEdges;
    int loops = 0;
    std::vector<int> loopPlane = LoopPlanes(s, planes, 1e-6, loops);
    printf("  layer %d (%.2f of the wall): %lld points, %lld triangles, mean edge %.3f, %d rims; distance from the interface %.3f..%.3f of its own offset\n",
        k, fraction, (ll)(s.points.size()/3), (ll)(s.triangles.size()/3), meanEdge, loops, worst, best);
    Check(loops == 2, "the layer is trimmed to two rims");
    Check(worst > 0.9 && best < 1.1, "the layer stands its own fraction of the wall off the interface");
    Check(meanEdge > 0.12 && meanEdge < 0.3, "the layer is decimated to about the interface's size, not its own thinner offset");
    ll base = (ll)(allPts.size()/3);
    allPts.insert(allPts.end(), s.points.begin(), s.points.end());
    for (size_t i = 0; i < s.triangles.size(); i++) allTris.push_back(s.triangles[i] + base);
    levels.push_back(s);
  }
  if (!ok) return;
  std::vector<unsigned char> crossing;
  double at[3];
  ll numCrossing = svenvelope::CountCrossingTriangles(allPts, allTris, crossing, at);
  printf("  all layers with the interface: %lld crossing triangles\n", numCrossing);
  Check(numCrossing == 0, "no layer surface crosses another, the outer surface or the interface");
}

int main()
{
  TestTube();
  TestThinTube();
  TestSeptum();
  TestNeighbouringEnds();
  TestLayers();
  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
