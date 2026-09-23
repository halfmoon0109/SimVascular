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

// An end-to-end test of the wall fill without VTK: the interface of a
// synthetic junction (a thin branch on a thick parent) goes through the
// offset core (outer surface and layer surfaces from one field), the trim at
// the cap planes, the wall shell (interface reversed, layer facets, outer,
// and an annulus by angle between consecutive rims at each cap, as
// TGenUtils_BuildWallShellSurface stitches them) and TetGen with the switches
// FillWallMeshWithTetGen uses (-p, -Y, nojettison, -q, -nn), and what the fill
// flow then relies on is checked: the shell closed and free of crossings, the
// mesher accepting it, every interface triangle back as a boundary face on its
// own three points, every layer facet back inside the wall, and the wall's
// volume about the interface's area times the thickness. Written 2026-09-23
// after the layer surfaces, built each from its own scaled interface, crossed
// the outer surface on the user's model and TetGen refused the shell - which
// no test of the surfaces one at a time could have shown.
//
// Build and run (one command line; TetGen takes a while to compile):
//   g++ -O2 -std=c++17 -DTETLIBRARY -I. -I../../../../ThirdParty/tetgen/simvascular_tetgen
//     Testing/test_fill.cxx sv_tetgenmesh_offset.cxx sv_tetgenmesh_envelope.cxx
//     ../../../../ThirdParty/tetgen/simvascular_tetgen/tetgen.cxx
//     ../../../../ThirdParty/tetgen/simvascular_tetgen/predicates.cxx -o test_fill && ./test_fill
#include "sv_tetgenmesh_offset.h"
#include "sv_tetgenmesh_envelope.h"
#ifndef TETLIBRARY
#define TETLIBRARY
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
typedef long long ll;
using namespace svoffset;

static int numFailed = 0;
static void Check(bool ok, const char *what)
{
  printf("  %s   %s\n", ok ? "ok" : "FAIL", what);
  if (!ok) numFailed++;
}
static double Dot(const double *a, const double *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static void Sub(const double *a, const double *b, double *r) { r[0] = a[0]-b[0]; r[1] = a[1]-b[1]; r[2] = a[2]-b[2]; }
static void Cross(const double *a, const double *b, double *r) { r[0] = a[1]*b[2]-a[2]*b[1]; r[1] = a[2]*b[0]-a[0]*b[2]; r[2] = a[0]*b[1]-a[1]*b[0]; }
static double Norm(const double *a) { return std::sqrt(Dot(a, a)); }
static double Dist(const double *a, const double *b) { double d[3]; Sub(a, b, d); return Norm(d); }

static bool Delaunay(const std::vector<double> &points, std::vector<ll> &tets, void *, std::string &error)
{
  tetgenio in, out; in.firstnumber = 0; in.numberofpoints = (int)(points.size()/3); in.pointlist = new REAL[points.size()];
  for (size_t i = 0; i < points.size(); i++) in.pointlist[i] = points[i];
  char sw[] = "Q"; try { tetrahedralize(sw, &in, &out); } catch (...) { error = "tetgen"; return false; }
  std::map<std::array<ll,3>, ll> by; auto key = [](const double *p) { std::array<ll,3> k = {(ll)std::llround(p[0]*1e6), (ll)std::llround(p[1]*1e6), (ll)std::llround(p[2]*1e6)}; return k; };
  for (size_t i = 0; i < points.size()/3; i++) by[key(&points[3*i])] = (ll)i;
  std::vector<ll> to(out.numberofpoints, -1);
  for (int i = 0; i < out.numberofpoints; i++) { auto it = by.find(key(&out.pointlist[3*i])); if (it == by.end()) { error = "point"; return false; } to[i] = it->second; }
  tets.resize(4*(size_t)out.numberoftetrahedra); for (int t = 0; t < out.numberoftetrahedra; t++) for (int m = 0; m < 4; m++) tets[4*t+m] = to[out.tetrahedronlist[4*t+m]];
  return true;
}


static std::vector<CapPlane> PlanesFromRims(const Interface &iface, const Surface &s)
{
  std::vector<CapPlane> planes;
  for (size_t r = 0; r < s.rims.size(); r++)
  {
    const std::vector<ll> &loop = s.rims[r];
    CapPlane plane; plane.rim = (ll)r;
    for (int k = 0; k < 3; k++) plane.origin[k] = 0.0;
    for (size_t m = 0; m < loop.size(); m++) for (int k = 0; k < 3; k++) plane.origin[k] += iface.points[3*(size_t)loop[m] + k]/(double)loop.size();
    double normal[3] = {0.0, 0.0, 0.0}, radius = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
    {
      const double *a = &iface.points[3*(size_t)loop[m]], *b = &iface.points[3*(size_t)loop[(m+1)%loop.size()]];
      normal[0] += (a[1]-b[1])*(a[2]+b[2]); normal[1] += (a[2]-b[2])*(a[0]+b[0]); normal[2] += (a[0]-b[0])*(a[1]+b[1]);
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

// directed boundary loops (edges used by one triangle, traversed the way that triangle does)

static std::vector<std::vector<ll> > BoundaryLoops(const std::vector<ll> &tris)
{
  std::map<std::pair<ll,ll>, int> count;
  for (size_t i = 0; i + 2 < tris.size(); i += 3) for (int j = 0; j < 3; j++) { ll a = tris[i+j], b = tris[i+(j+1)%3]; count[std::make_pair(std::min(a,b), std::max(a,b))]++; }
  std::map<ll,ll> next;
  for (size_t i = 0; i + 2 < tris.size(); i += 3) for (int j = 0; j < 3; j++) { ll a = tris[i+j], b = tris[i+(j+1)%3]; if (count[std::make_pair(std::min(a,b), std::max(a,b))] == 1) next[a] = b; }
  std::set<ll> seen; std::vector<std::vector<ll> > loops;
  for (auto it = next.begin(); it != next.end(); ++it)
  {
    if (seen.count(it->first)) continue;
    std::vector<ll> loop; ll cur = it->first;
    while (!seen.count(cur)) { seen.insert(cur); loop.push_back(cur); auto nx = next.find(cur); if (nx == next.end()) break; cur = nx->second; }
    loops.push_back(loop);
  }
  return loops;
}

// one loop per plane: the loop whose centre is nearest the plane origin, among loops lying on the plane within tol

static bool MatchLoops(const std::vector<double> &pts, const std::vector<std::vector<ll> > &loops, const std::vector<CapPlane> &planes, double tol, std::vector<std::vector<ll> > &perPlane, const char *what)
{
  perPlane.assign(planes.size(), std::vector<ll>());
  std::vector<int> assigned(loops.size(), -1);
  for (size_t p = 0; p < planes.size(); p++)
  {
    int best = -1; double nearest = 1e300;
    for (size_t l = 0; l < loops.size(); l++)
    {
      if (assigned[l] >= 0) continue;
      double centre[3] = {0,0,0}, worst = 0.0;
      for (size_t m = 0; m < loops[l].size(); m++) for (int k = 0; k < 3; k++) centre[k] += pts[3*loops[l][m] + k]/(double)loops[l].size();
      for (size_t m = 0; m < loops[l].size(); m++) { double off[3]; Sub(&pts[3*loops[l][m]], planes[p].origin, off); worst = std::max(worst, std::fabs(Dot(off, planes[p].outward))); }
      if (worst > tol) continue;
      double d = Dist(centre, planes[p].origin);
      if (d < nearest) { nearest = d; best = (int)l; }
    }
    if (best < 0) { printf("  %s: no rim loop on cap plane %zu at (%.3f, %.3f, %.3f) within %.3g\n", what, p, planes[p].origin[0], planes[p].origin[1], planes[p].origin[2], tol); return false; }
    assigned[best] = (int)p; perPlane[p] = loops[best];
  }
  int unassigned = 0; for (size_t l = 0; l < loops.size(); l++) if (assigned[l] < 0) unassigned++;
  if (unassigned) printf("  %s: %d rim loops matched no cap plane\n", what, unassigned);
  return unassigned == 0;
}

// port of TGenUtils_StitchCapAnnulus

static bool StitchAnnulus(const std::vector<double> &pts, std::vector<ll> inner, std::vector<ll> outerRim, const double outward[3], std::vector<ll> &cells, int &numDegenerate)
{
  if (inner.size() < 3 || outerRim.size() < 3) { printf("  annulus: rims of %zu and %zu points\n", inner.size(), outerRim.size()); return false; }
  double center[3] = {0,0,0};
  for (size_t m = 0; m < inner.size(); m++) for (int k = 0; k < 3; k++) center[k] += pts[3*inner[m] + k]/(double)inner.size();
  double axis[3] = {0,0,0}; int smallest = 0; for (int k = 1; k < 3; k++) if (std::fabs(outward[k]) < std::fabs(outward[smallest])) smallest = k; axis[smallest] = 1.0;
  double u[3], v[3]; Cross(axis, outward, u); double lu = Norm(u); if (lu <= 0) return false; for (int k = 0; k < 3; k++) u[k] /= lu;
  Cross(outward, u, v); double lv = Norm(v); if (lv <= 0) return false; for (int k = 0; k < 3; k++) v[k] /= lv;
  auto angleOf = [&](ll id) { double off[3]; Sub(&pts[3*id], center, off); return std::atan2(Dot(off, v), Dot(off, u)); };
  auto wrap = [](double d) { while (d > M_PI) d -= 2*M_PI; while (d <= -M_PI) d += 2*M_PI; return d; };
  for (int side = 0; side < 2; side++)
  {
    std::vector<ll> &loop = side == 0 ? inner : outerRim; double turning = 0.0;
    for (size_t m = 0; m < loop.size(); m++) turning += wrap(angleOf(loop[(m+1)%loop.size()]) - angleOf(loop[m]));
    if (std::fabs(std::fabs(turning) - 2*M_PI) > 0.5) { printf("  annulus: a rim of %zu points turns %.3f radians\n", loop.size(), turning); return false; }
    if (turning < 0) std::reverse(loop.begin(), loop.end());
  }
  size_t n = inner.size(), m = outerRim.size(), startOuter = 0; double startAngle = angleOf(inner[0]), bestGap = 0;
  for (size_t j = 0; j < m; j++) { double gap = std::fabs(wrap(angleOf(outerRim[j]) - startAngle)); if (j == 0 || gap < bestGap) { bestGap = gap; startOuter = j; } }
  std::vector<double> innerSweep(n+1, 0.0), outerSweep(m+1, 0.0);
  for (size_t k = 0; k < n; k++) innerSweep[k+1] = innerSweep[k] + wrap(angleOf(inner[(k+1)%n]) - angleOf(inner[k]));
  for (size_t k = 0; k < m; k++) outerSweep[k+1] = outerSweep[k] + wrap(angleOf(outerRim[(startOuter+k+1)%m]) - angleOf(outerRim[(startOuter+k)%m]));
  for (size_t k = 0; k <= n; k++) innerSweep[k] /= innerSweep[n];
  for (size_t k = 0; k <= m; k++) outerSweep[k] /= outerSweep[m];
  size_t i = 0, j = 0;
  while (i < n || j < m)
  {
    bool advanceInner = (i >= n) ? false : (j >= m) ? true : (innerSweep[i+1] <= outerSweep[j+1]);
    ll tri[3];
    if (advanceInner) { tri[0] = inner[i%n]; tri[1] = inner[(i+1)%n]; tri[2] = outerRim[(startOuter+j)%m]; i++; }
    else { tri[0] = outerRim[(startOuter+j)%m]; tri[1] = outerRim[(startOuter+j+1)%m]; tri[2] = inner[i%n]; j++; }
    double e1[3], e2[3], nrm[3]; Sub(&pts[3*tri[1]], &pts[3*tri[0]], e1); Sub(&pts[3*tri[2]], &pts[3*tri[0]], e2); Cross(e1, e2, nrm);
    if (Norm(nrm) <= 0.0) numDegenerate++; else if (Dot(nrm, outward) < 0.0) std::swap(tri[1], tri[2]);
    cells.insert(cells.end(), tri, tri + 3);
  }
  return true;
}


// A parent tube of radius 1 along z from 0 to 12 with a branch of radius
// rBranch along +x at z = 6, the hole in the parent stitched to the branch's
// base ring by angle; the parent's wall is tParent, the branch's tBranch.
static void MakeJunction(double tParent, double tBranch, double rBranch, int aroundP, int aroundB, Interface &f)
{
  int alongP = 48, alongB = 24;
  // parent tube along z
  for (int j = 0; j <= alongP; j++) for (int i = 0; i < aroundP; i++)
  {
    double th = 2*M_PI*i/aroundP, z = 12.0*j/alongP;
    f.points.push_back(cos(th)); f.points.push_back(sin(th)); f.points.push_back(z);
    f.normals.push_back(cos(th)); f.normals.push_back(sin(th)); f.normals.push_back(0);
    f.thickness.push_back(tParent);
  }
  std::vector<ll> parentTris;
  double hole = rBranch + 0.15;
  auto inHole = [&](ll v) { const double *p = &f.points[3*v]; return p[0] > 0 && p[1]*p[1] + (p[2]-6.0)*(p[2]-6.0) < hole*hole; };
  for (int j = 0; j < alongP; j++) for (int i = 0; i < aroundP; i++)
  {
    ll a = j*aroundP + i, b = j*aroundP + (i+1)%aroundP, c = (j+1)*aroundP + i, d = (j+1)*aroundP + (i+1)%aroundP;
    ll t1[3] = {a, b, d}, t2[3] = {a, d, c};
    if (!(inHole(a) || inHole(b) || inHole(d))) parentTris.insert(parentTris.end(), t1, t1+3);
    if (!(inHole(a) || inHole(d) || inHole(c))) parentTris.insert(parentTris.end(), t2, t2+3);
  }
  // the hole's boundary loop (edges used once among parent triangles, on the hole side)
  std::map<std::pair<ll,ll>, int> cnt; std::map<ll,ll> next;
  for (size_t i = 0; i + 2 < parentTris.size(); i += 3) for (int j = 0; j < 3; j++) { ll a = parentTris[i+j], b = parentTris[i+(j+1)%3]; cnt[std::make_pair(std::min(a,b), std::max(a,b))]++; }
  for (size_t i = 0; i + 2 < parentTris.size(); i += 3) for (int j = 0; j < 3; j++) { ll a = parentTris[i+j], b = parentTris[i+(j+1)%3]; if (cnt[std::make_pair(std::min(a,b), std::max(a,b))] == 1) { const double *p = &f.points[3*a]; if (p[2] > 0.01 && p[2] < 11.99) next[a] = b; } }
  std::vector<ll> holeLoop; { ll start = next.begin()->first, cur = start; do { holeLoop.push_back(cur); cur = next[cur]; } while (cur != start && holeLoop.size() < 10000); }
  // branch tube along +x from the parent surface
  ll base = f.points.size()/3;
  for (int j = 0; j <= alongB; j++) for (int i = 0; i < aroundB; i++)
  {
    double ph = 2*M_PI*i/aroundB, y = rBranch*sin(ph), z = 6.0 + rBranch*cos(ph);
    double x0 = std::sqrt(std::max(0.0, 1.0 - y*y)); double x = x0 + (4.0 - x0)*j/alongB;
    f.points.push_back(x); f.points.push_back(y); f.points.push_back(z);
    f.normals.push_back(0); f.normals.push_back(sin(ph)); f.normals.push_back(cos(ph));
    f.thickness.push_back(tBranch);
  }
  std::vector<ll> tris = parentTris;
  for (int j = 0; j < alongB; j++) for (int i = 0; i < aroundB; i++)
  {
    ll a = base + j*aroundB + i, b = base + j*aroundB + (i+1)%aroundB, c = base + (j+1)*aroundB + i, d = base + (j+1)*aroundB + (i+1)%aroundB;
    // outward about the x axis: normal = (0, sin, cos); winding (a, c, d), (a, d, b)? check: at ph=0 (top, z=6.3), i+1 has larger y; along +x is c. (a, b, d): edge a->b (+y), b->d (+x): cross(+y,+x) = -z -> inward. so use (a, d, b) and (a, c, d).
    tris.push_back(a); tris.push_back(d); tris.push_back(b);
    tris.push_back(a); tris.push_back(c); tris.push_back(d);
  }
  // stitch the hole loop to the branch's base ring by angle about the branch axis
  auto angleOf = [&](ll v) { const double *p = &f.points[3*v]; return atan2(p[2]-6.0, p[1]); };
  std::vector<ll> ring; for (int i = 0; i < aroundB; i++) ring.push_back(base + i);
  auto byAngle = [&](std::vector<ll> &loop) { std::sort(loop.begin(), loop.end(), [&](ll u, ll v) { return angleOf(u) < angleOf(v); }); };
  byAngle(holeLoop); byAngle(ring);
  std::vector<ll> band;
  {
    size_t n = holeLoop.size(), m = ring.size(); size_t i = 0, j = 0;
    auto wrap = [](double a) { while (a > M_PI) a -= 2*M_PI; while (a < -M_PI) a += 2*M_PI; return a; };
    while (i < n || j < m)
    {
      ll h0 = holeLoop[i%n], h1 = holeLoop[(i+1)%n], r0 = ring[j%m], r1 = ring[(j+1)%m];
      bool advanceHole;
      if (i >= n) advanceHole = false; else if (j >= m) advanceHole = true;
      else advanceHole = wrap(angleOf(h1) - angleOf(h0)) + angleOf(h0) <= wrap(angleOf(r1) - angleOf(r0)) + angleOf(r0);
      if (advanceHole) { band.push_back(h0); band.push_back(h1); band.push_back(r0); i++; }
      else { band.push_back(h0); band.push_back(r1); band.push_back(r0); j++; }
    }
  }
  // orient the band against the parent's boundary edges: the parent traverses a->b, the band must traverse b->a
  {
    std::set<std::pair<ll,ll> > parentDir; for (size_t i = 0; i + 2 < parentTris.size(); i += 3) for (int j = 0; j < 3; j++) parentDir.insert(std::make_pair(parentTris[i+j], parentTris[i+(j+1)%3]));
    int same = 0, opposite = 0;
    for (size_t i = 0; i + 2 < band.size(); i += 3) for (int j = 0; j < 3; j++) { ll a = band[i+j], b = band[i+(j+1)%3]; if (parentDir.count(std::make_pair(a,b))) same++; if (parentDir.count(std::make_pair(b,a))) opposite++; }
    if (same > opposite) for (size_t i = 0; i + 2 < band.size(); i += 3) std::swap(band[i+1], band[i+2]);
    printf("band: %zu triangles (hole loop %zu, ring %zu), edges same %d opposite %d%s\n", band.size()/3, holeLoop.size(), ring.size(), same, opposite, same > opposite ? " -> flipped" : "");
  }
  tris.insert(tris.end(), band.begin(), band.end());
  f.triangles = tris;

}

struct Built { Surface surface; std::vector<std::vector<ll> > rimPerCap; TrimReport trim; };

static bool FinishLevel(const Interface &iface, const std::vector<CapPlane> *planesIn, std::vector<CapPlane> &planesOut, Built &b, const char *label)
{
  std::string err;
  planesOut = planesIn ? *planesIn : PlanesFromRims(iface, b.surface);
  if (planesIn) for (size_t p = 0; p < planesOut.size(); p++) planesOut[p].rim = (ll)p;
  if (TrimSurfaceAtCaps(b.surface, planesOut, 0.1, b.trim, err) != 0) { printf("  %s: trim failed: %s\n", label, err.c_str()); return false; }
  ll nb, nn, nm; CountEdges(b.surface.triangles, nb, nn, nm);
  std::vector<unsigned char> cr; double at[3]; ll crossings = svenvelope::CountCrossingTriangles(b.surface.points, b.surface.triangles, cr, at);
  std::vector<std::vector<ll> > loops = BoundaryLoops(b.surface.triangles);
  printf("  %s: %lld points, %lld triangles after the trim, %zu rim loops, %lld non-manifold, %lld miswound, %lld crossing; %lld rim edges merged, shortest rim edge %.3f of the mean, rim angle %.1f degrees\n",
      label, (ll)(b.surface.points.size()/3), (ll)(b.surface.triangles.size()/3), loops.size(), nn, nm, crossings, b.trim.numRimEdgesMerged, b.trim.shortestRimEdgeRatio, b.trim.smallestRimAngleDegrees);
  Check(nn == 0 && nm == 0 && crossings == 0, "the trimmed surface is a manifold without crossings");
  Check(loops.size() == planesOut.size(), "the trimmed surface has one rim per cap");
  if (nn || nm || crossings || loops.size() != planesOut.size()) return false;
  double tol = 0.0; for (size_t i = 0; i < iface.thickness.size(); i++) tol = std::max(tol, iface.thickness[i]); tol *= 0.02;
  return MatchLoops(b.surface.points, loops, planesOut, tol, b.rimPerCap, label);
}

// The whole fill of one junction with numLayers layers; true if everything held.
static void FillJunction(double tParent, double tBranch, double rBranch, int numLayers)
{
  printf("junction: parent wall %.2f, branch wall %.2f and radius %.2f, %d layer(s)\n", tParent, tBranch, rBranch, numLayers);
  Interface iface;
  MakeJunction(tParent, tBranch, rBranch, 48, 16, iface);
  ll numPts = (ll)(iface.points.size()/3);
  std::vector<double> fractions; for (int k = 1; k <= numLayers; k++) fractions.push_back((double)k/numLayers);
  std::vector<Surface> surfs; std::vector<Report> reps; std::string err; Options options;
  if (BuildOffsetSurfaces(iface, options, fractions, Delaunay, nullptr, surfs, reps, err) != 0) { Check(false, "the offset surfaces are built"); printf("  %s\n", err.c_str()); return; }
  Check(true, "the offset surfaces are built");
  Built outer; std::vector<CapPlane> planes;
  outer.surface = surfs.back();
  if (!FinishLevel(iface, nullptr, planes, outer, "outer")) return;
  std::vector<Built> levels(numLayers - 1);
  for (int k = 1; k < numLayers; k++)
  {
    char label[64]; snprintf(label, sizeof(label), "layer %d/%d", k, numLayers);
    levels[k-1].surface = surfs[k-1]; std::vector<CapPlane> lp;
    if (!FinishLevel(iface, &planes, lp, levels[k-1], label)) return;
  }
  // the shell: interface reversed | layers | outer, annuli between consecutive rims at each cap
  std::vector<double> pts = iface.points; std::vector<ll> tris; std::vector<int> role;
  std::vector<ll> levelBase(levels.size(), 0);
  for (size_t k = 0; k < levels.size(); k++) { levelBase[k] = (ll)(pts.size()/3); pts.insert(pts.end(), levels[k].surface.points.begin(), levels[k].surface.points.end()); }
  ll outerBase = (ll)(pts.size()/3); pts.insert(pts.end(), outer.surface.points.begin(), outer.surface.points.end());
  for (size_t i = 0; i + 2 < iface.triangles.size(); i += 3) { tris.push_back(iface.triangles[i+2]); tris.push_back(iface.triangles[i+1]); tris.push_back(iface.triangles[i]); role.push_back(1); }
  for (size_t k = 0; k < levels.size(); k++) for (size_t i = 0; i + 2 < levels[k].surface.triangles.size(); i += 3) { for (int j = 0; j < 3; j++) tris.push_back(levels[k].surface.triangles[i+j] + levelBase[k]); role.push_back(100 + (int)k + 1); }
  {
    // the outer faces away from the interface, or is reversed, as the shell builder decides by a vote
    int agree = 0, disagree = 0; size_t nOut = outer.surface.triangles.size()/3;
    for (size_t t = 0; t < nOut; t += 4)
    {
      const ll *tt = &outer.surface.triangles[3*t]; const double *p0 = &outer.surface.points[3*tt[0]], *p1 = &outer.surface.points[3*tt[1]], *p2 = &outer.surface.points[3*tt[2]];
      double e1[3], e2[3], nrm[3], c[3]; Sub(p1, p0, e1); Sub(p2, p0, e2); Cross(e1, e2, nrm); if (Norm(nrm) <= 0) continue;
      for (int k = 0; k < 3; k++) c[k] = (p0[k] + p1[k] + p2[k])/3.0;
      double best = 1e300; ll bi = -1; for (ll i = 0; i < numPts; i++) { double d = Dist(c, &iface.points[3*i]); if (d < best) { best = d; bi = i; } }
      double away[3]; Sub(c, &iface.points[3*bi], away); if (Dot(nrm, away) >= 0) agree++; else disagree++;
    }
    bool reverse = disagree > agree;
    for (size_t i = 0; i + 2 < outer.surface.triangles.size(); i += 3) { const ll *tt = &outer.surface.triangles[i]; if (reverse) { tris.push_back(tt[2] + outerBase); tris.push_back(tt[1] + outerBase); tris.push_back(tt[0] + outerBase); } else { tris.push_back(tt[0] + outerBase); tris.push_back(tt[1] + outerBase); tris.push_back(tt[2] + outerBase); } role.push_back(2); }
  }
  int numDegenerate = 0;
  for (size_t p = 0; p < planes.size(); p++)
  {
    std::vector<std::vector<ll> > rims; rims.push_back(outer.surface.rims[p]);
    for (size_t k = 0; k < levels.size(); k++) { std::vector<ll> shifted = levels[k].rimPerCap[p]; for (size_t m = 0; m < shifted.size(); m++) shifted[m] += levelBase[k]; rims.push_back(shifted); }
    { std::vector<ll> shifted = outer.rimPerCap[p]; for (size_t m = 0; m < shifted.size(); m++) shifted[m] += outerBase; rims.push_back(shifted); }
    for (size_t r = 0; r + 1 < rims.size(); r++)
    {
      size_t before = tris.size();
      if (!StitchAnnulus(pts, rims[r], rims[r+1], planes[p].outward, tris, numDegenerate)) { Check(false, "the annulus between consecutive rims is stitched"); return; }
      for (size_t a = before; a < tris.size(); a += 3) role.push_back(9999);
    }
  }
  ll nb, nn, nm; CountEdges(tris, nb, nn, nm);
  std::vector<unsigned char> cr; double at[3]; ll crossings = svenvelope::CountCrossingTriangles(pts, tris, cr, at);
  ll layerRimEdges = 0; for (size_t k = 0; k < levels.size(); k++) { ll b2, n2, m2; CountEdges(levels[k].surface.triangles, b2, n2, m2); layerRimEdges += b2; }
  printf("  shell: %lld points, %zu triangles; %lld boundary edges, %lld on three facets (%lld layer rim edges), %lld miswound, %lld crossing, %d annulus triangles of no area\n",
      (ll)(pts.size()/3), tris.size()/3, nb, nn, layerRimEdges, nm, crossings, numDegenerate);
  Check(nb == 0 && nm == 0 && nn == layerRimEdges && crossings == 0 && numDegenerate == 0, "the shell is closed, wound consistently and free of crossings");
  if (nb || nm || crossings) return;
  // TetGen as the fill flow calls it
  tetgenio in, out; in.firstnumber = 0; in.numberofpoints = (int)(pts.size()/3); in.pointlist = new REAL[pts.size()]; for (size_t i = 0; i < pts.size(); i++) in.pointlist[i] = pts[i];
  in.numberoffacets = (int)(tris.size()/3); in.facetlist = new tetgenio::facet[in.numberoffacets]; in.facetmarkerlist = new int[in.numberoffacets]();
  for (int i = 0; i < in.numberoffacets; i++) { tetgenio::facet *f = &in.facetlist[i]; f->numberofpolygons = 1; f->polygonlist = new tetgenio::polygon[1]; f->numberofholes = 0; f->holelist = nullptr; tetgenio::polygon *pg = &f->polygonlist[0]; pg->numberofvertices = 3; pg->vertexlist = new int[3]; for (int j = 0; j < 3; j++) pg->vertexlist[j] = (int)tris[3*i+j]; in.facetmarkerlist[i] = role[i]; }
  tetgenbehavior b; b.plc = 1; b.nobisect = 1; b.nojettison = 1; b.quality = 1; b.minratio = numLayers > 1 ? 2.0 : 1.414; b.mindihedral = numLayers > 1 ? 5.0 : 10.0; b.neighout = 2; b.quiet = 1;
  bool accepted = true;
  try { tetrahedralize(&b, &in, &out); } catch (int r) { printf("  TetGen error %d\n", r); accepted = false; }
  Check(accepted, "TetGen accepts the shell");
  if (!accepted) return;
  int nIn = in.numberofpoints;
  std::set<std::array<ll,3> > interfaceSet; for (size_t i = 0; i + 2 < iface.triangles.size(); i += 3) { std::array<ll,3> k = {iface.triangles[i], iface.triangles[i+1], iface.triangles[i+2]}; std::sort(k.begin(), k.end()); interfaceSet.insert(k); }
  std::map<int, ll> facesByMarker; ll onSteiner = 0; std::set<std::array<ll,3> > seenInterface; ll interfaceOff = 0;
  for (int i = 0; i < out.numberoftrifaces; i++)
  {
    int mk = out.trifacemarkerlist ? out.trifacemarkerlist[i] : -1; facesByMarker[mk]++;
    ll a = out.trifacelist[3*i], c2 = out.trifacelist[3*i+1], c = out.trifacelist[3*i+2];
    if (a >= nIn || c2 >= nIn || c >= nIn) onSteiner++;
    if (mk == 1) { std::array<ll,3> k = {a, c2, c}; std::sort(k.begin(), k.end()); if (interfaceSet.count(k)) seenInterface.insert(k); else interfaceOff++; }
  }
  double volume = 0.0, expected = 0.0, minDihedral = 180.0; ll below5 = 0;
  for (int t = 0; t < out.numberoftetrahedra; t++)
  {
    const int *q = &out.tetrahedronlist[4*t]; const double *P[4]; for (int m = 0; m < 4; m++) P[m] = &out.pointlist[3*q[m]];
    double e1[3], e2[3], e3[3], crs[3]; Sub(P[1], P[0], e1); Sub(P[2], P[0], e2); Sub(P[3], P[0], e3); Cross(e2, e3, crs); volume += std::fabs(Dot(e1, crs))/6.0;
    const int faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}}; double N[4][3];
    for (int f = 0; f < 4; f++) { double u[3], v[3]; Sub(P[faces[f][1]], P[faces[f][0]], u); Sub(P[faces[f][2]], P[faces[f][0]], v); Cross(u, v, N[f]); double l = Norm(N[f]); if (l > 0) for (int k = 0; k < 3; k++) N[f][k] /= l; }
    double tetMin = 180.0; for (int f = 0; f < 4; f++) for (int g = f+1; g < 4; g++) { double cosang = std::max(-1.0, std::min(1.0, -Dot(N[f], N[g]))); tetMin = std::min(tetMin, std::acos(cosang)*180.0/M_PI); }
    minDihedral = std::min(minDihedral, tetMin); if (tetMin < 5.0) below5++;
  }
  for (size_t i = 0; i + 2 < iface.triangles.size(); i += 3) { const ll *tt = &iface.triangles[i]; double e1[3], e2[3], crs[3]; Sub(&iface.points[3*tt[1]], &iface.points[3*tt[0]], e1); Sub(&iface.points[3*tt[2]], &iface.points[3*tt[0]], e2); Cross(e1, e2, crs); expected += 0.5*Norm(crs)*(iface.thickness[tt[0]] + iface.thickness[tt[1]] + iface.thickness[tt[2]])/3.0; }
  printf("  TetGen: %d points in, %d out, %d tetrahedra, %d boundary faces (%lld on Steiner points); faces by marker:", nIn, out.numberofpoints, out.numberoftetrahedra, out.numberoftrifaces, onSteiner);
  for (std::map<int, ll>::iterator it = facesByMarker.begin(); it != facesByMarker.end(); ++it) printf(" %d:%lld", it->first, it->second);
  printf("; volume %.3f against area times thickness %.3f (%.3f); smallest dihedral angle %.2f degrees, %lld tetrahedra under 5\n", volume, expected, volume/expected, minDihedral, below5);
  Check(seenInterface.size() == interfaceSet.size() && interfaceOff == 0, "every interface triangle comes back as a boundary face on its own points, and no other face carries its marker");
  Check(onSteiner == 0, "no facet was split by a Steiner point");
  for (size_t k = 0; k < levels.size(); k++)
  {
    char what[96]; snprintf(what, sizeof(what), "every triangle of layer surface %zu comes back as a face inside the wall", k + 1);
    Check(facesByMarker[100 + (int)k + 1] == (ll)(levels[k].surface.triangles.size()/3), what);
  }
  Check(facesByMarker[2] == (ll)(outer.surface.triangles.size()/3), "every outer triangle comes back as a boundary face");
  // a tube's wall is thicker than area times thickness by half the thickness over the radius: 1.25 for the parent here
  Check(volume/expected > 1.1 && volume/expected < 1.4, "the wall's volume is the interface's area times the thickness, allowing for the curvature");
}

int main()
{
  FillJunction(0.5, 0.1, 0.3, 1);
  FillJunction(0.5, 0.1, 0.3, 3);
  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
