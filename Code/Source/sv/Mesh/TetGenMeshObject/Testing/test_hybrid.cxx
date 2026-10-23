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

// A standalone test of the prism zone classification behind the hybrid wall
// mesh (docs/wall-mesh-purpose.md section 8): where the wall can be prism
// layers along the normals and where it has to be filled with tetrahedra.
// No VTK and no TetGen: the field alone decides.
//
// Build and run (TetGen is compiled in for the Delaunay of the level surfaces):
//   g++ -O2 -std=c++17 -DTETLIBRARY -I. -I../../../../ThirdParty/tetgen/simvascular_tetgen
//     Testing/test_hybrid.cxx sv_tetgenmesh_offset.cxx sv_tetgenmesh_envelope.cxx
//     ../../../../ThirdParty/tetgen/simvascular_tetgen/tetgen.cxx
//     ../../../../ThirdParty/tetgen/simvascular_tetgen/predicates.cxx -o test_hybrid && ./test_hybrid
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
  printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok) numFailed++;
}

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


static bool Classify(const Interface &iface, int numLayers, std::vector<unsigned char> &structured, ZoneReport &zr)
{
  OffsetField field; Report fr; std::string err;
  if (field.Build(iface, fr, err) != 0) { printf("  FAIL field: %s\n", err.c_str()); numFailed++; return false; }
  ZoneOptions zo;
  if (ClassifyPrismZone(iface, field, numLayers, zo, structured, zr, err) != 0) { printf("  FAIL classify: %s\n", err.c_str()); numFailed++; return false; }
  printf("  zone: %lld of %lld triangles structured; points covered %lld, off level %lld; triangles inverted %lld, tops crossing %lld, margin %lld; islands %lld (%lld triangles); junction regions %lld; %lld evaluations\n",
      zr.numStructured, zr.numTriangles, zr.numPointsCovered, zr.numPointsOffLevel, zr.numTrianglesInverted, zr.numTrianglesCrossing, zr.numTrianglesMargin,
      zr.numIslands, zr.numIslandTriangles, zr.numJunctionRegions, zr.numEvaluations);
  return true;
}


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

static double Dot(const double *a, const double *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static void Sub(const double *a, const double *b, double *r) { r[0] = a[0]-b[0]; r[1] = a[1]-b[1]; r[2] = a[2]-b[2]; }
static double Norm(const double *a) { return std::sqrt(Dot(a, a)); }
static double Dist(const double *a, const double *b) { double d[3]; Sub(a, b, d); return Norm(d); }

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

// The directed zone boundary edges of the prism mesh chained into loops of interface point ids.
static void ZoneLoops(const PrismMesh &pm, std::vector<std::vector<ll> > &loops)
{
  std::map<ll, ll> next; for (size_t e = 0; e + 1 < pm.zoneBoundaryEdges.size(); e += 2) next[pm.zoneBoundaryEdges[e]] = pm.zoneBoundaryEdges[e+1];
  std::set<ll> used;
  for (std::map<ll, ll>::iterator it = next.begin(); it != next.end(); ++it)
  {
    if (used.count(it->first)) continue;
    std::vector<ll> loop; ll cur = it->first;
    while (!used.count(cur)) { used.insert(cur); loop.push_back(cur); cur = next[cur]; }
    loops.push_back(loop);
  }
}

// The junction zone's surface at one level - the layer piece zipped to the prism zone's layer surface - must close with it.
static void CheckLevelPiece(const Interface &iface, const OffsetField &field, const Surface &level, double fraction, int k,
    const std::vector<unsigned char> &structured, const PrismMesh &pm)
{
  std::vector<unsigned char> rimStructured(field.Rims().size(), 1);
  {
    const std::vector<std::vector<ll> > &rims = field.Rims();
    std::map<ll, std::vector<ll> > incident; for (size_t t = 0; t + 2 < iface.triangles.size(); t += 3) for (int j = 0; j < 3; j++) incident[iface.triangles[t+j]].push_back((ll)(t/3));
    for (size_t r = 0; r < rims.size(); r++) for (size_t m = 0; m < rims[r].size(); m++) for (ll t : incident[rims[r][m]]) if (!structured[t]) rimStructured[r] = 0;
  }
  Surface piece; ZonePieceReport zp; std::string err;
  const std::vector<ll> &guard = pm.layerTriangles[k-1];
  // as BuildJunctionShell trims: the layer's prism tetrahedra for the triangles hanging over the zone, the lifted loops for the erosion by distance (one loop edge)
  const ll np = (ll)(iface.points.size()/3);
  std::vector<ll> layerTets; for (size_t t = 0; t + 3 < pm.tetrahedra.size(); t += 4) if (pm.tetrahedronLayer[t/4] == k) for (int m = 0; m < 4; m++) layerTets.push_back(pm.tetrahedra[t + m]);
  std::vector<double> loopSegments;
  for (size_t e = 0; e + 1 < pm.zoneBoundaryEdges.size(); e += 2) { ll la = pm.layerPoint[(size_t)(k-1)*np + pm.zoneBoundaryEdges[e]], lb = pm.layerPoint[(size_t)(k-1)*np + pm.zoneBoundaryEdges[e+1]]; if (la < 0 || lb < 0) continue; for (int j = 0; j < 3; j++) loopSegments.push_back(pm.points[3*la+j]); for (int j = 0; j < 3; j++) loopSegments.push_back(pm.points[3*lb+j]); }
  if (TrimSurfaceToZone(level, fraction, field, (ll)(iface.triangles.size()/3), structured, rimStructured, 1, pm.points, guard, 2, piece, zp, err, nullptr, nullptr, &layerTets, &loopSegments, 1.0) != 0) { printf("  FAIL trim: %s\n", err.c_str()); numFailed++; return; }
  printf("  level %d: %lld triangles, %lld owned by the prism zone, %lld eroded, %lld crossing, %lld ears, %lld kept in the piece, %lld boundary chains\n", k, zp.numTriangles, zp.numStructuredOwned, zp.numEroded, zp.numCrossing, zp.numEars, zp.numKept, zp.numChains);
  // one point array: the prism mesh's points, then the piece's
  std::vector<double> pts = pm.points; ll base = (ll)(pts.size()/3); pts.insert(pts.end(), piece.points.begin(), piece.points.end());
  std::vector<ll> tris = guard;
  for (size_t i = 0; i < piece.triangles.size(); i++) tris.push_back(piece.triangles[i] + base);
  // the chains: A from the zone boundary at this level, B from the piece
  std::vector<std::vector<ll> > loopsA; ZoneLoops(pm, loopsA);
  std::vector<std::vector<ll> > chainsB; std::vector<unsigned char> closedB; BoundaryChains(piece.triangles, chainsB, closedB);
  ll numZip = 0; bool zipped = true;
  for (size_t a = 0; a < loopsA.size(); a++)
  {
    // the zone loop lifted to this layer, reversed: the prism layer traverses it the interface's way, the strip must go the other
    std::vector<ll> A; for (size_t m = loopsA[a].size(); m > 0; m--) A.push_back(pm.layerPoint[(size_t)(k-1)*(iface.points.size()/3) + loopsA[a][m-1]]);
    // the nearest closed chain of the piece
    size_t bestB = chainsB.size(); double best = 1e300;
    for (size_t b = 0; b < chainsB.size(); b++)
    {
      if (!closedB[b]) continue;
      double d = 0; for (size_t m = 0; m < A.size(); m++) { double dm = 1e300; for (size_t q = 0; q < chainsB[b].size(); q++) dm = std::min(dm, Dist(&pts[3*A[m]], &piece.points[3*chainsB[b][q]])); d += dm; }
      if (d < best) { best = d; bestB = b; }
    }
    if (bestB == chainsB.size()) { zipped = false; continue; }
    std::vector<ll> B; for (size_t q = 0; q < chainsB[bestB].size(); q++) B.push_back(chainsB[bestB][q] + base);
    std::vector<ll> strip;
    if (ZipChains(pts, A, B, true, strip, err) != 0) { printf("  FAIL zip: %s\n", err.c_str()); zipped = false; continue; }
    numZip += (ll)(strip.size()/3);
    tris.insert(tris.end(), strip.begin(), strip.end());
  }
  ll nb, nn, nm; CountEdges(tris, nb, nn, nm);
  std::vector<unsigned char> cr; double at[3]; ll crossings = svenvelope::CountCrossingTriangles(pts, tris, cr, at);
  // the rim edges of the prism zone and the piece at this level are the only boundary allowed
  ll rimEdges = 0;
  for (size_t r = 0; r < field.Rims().size(); r++) rimEdges += (ll)field.Rims()[r].size();
  printf("  level %d closed with the prism zone: %lld zipper triangles, boundary %lld (cap rims %lld), non-manifold %lld, miswound %lld, crossing %lld%s\n", k, numZip, nb, rimEdges, nn, nm, crossings, crossings ? " first at" : "");
  if (crossings)
  {
    printf("    (%.3f, %.3f, %.3f)\n", at[0], at[1], at[2]);
    std::vector<svenvelope::CrossingPair> pairs; svenvelope::ListCrossingPairs(pts, tris, 12, pairs);
    const ll nGuard = (ll)(guard.size()/3), nPiece = (ll)(piece.triangles.size()/3);
    auto kind = [&](ll t) { return t < nGuard ? "prism" : (t < nGuard + nPiece ? "piece" : "strip"); };
    for (size_t q = 0; q < pairs.size(); q++)
      printf("    pair: %s %lld x %s %lld at (%.3f, %.3f, %.3f), length %.4f\n", kind(pairs[q].a), pairs[q].a, kind(pairs[q].b), pairs[q].b, 0.5*(pairs[q].from[0]+pairs[q].to[0]), 0.5*(pairs[q].from[1]+pairs[q].to[1]), 0.5*(pairs[q].from[2]+pairs[q].to[2]), Dist(pairs[q].from, pairs[q].to));
  }
  Check(zipped, "every zone loop found a piece chain to zip to");
  Check(nb == rimEdges && nn == 0 && nm == 0, "the layer surface over the junction, zipped to the prism zone's, closes up to the cap rims and is wound consistently");
  Check(crossings == 0, "and has no crossings");
}


// the farthest a junction-zone triangle's centre is from a point
static double FarthestJunction(const Interface &iface, const std::vector<unsigned char> &structured, const double c[3])
{
  double far = 0.0;
  for (size_t t = 0; t + 2 < iface.triangles.size()/3*3; t += 3)
  {
    if (structured[t/3]) continue;
    const ll *tt = &iface.triangles[t];
    double cen[3] = {0, 0, 0};
    for (int j = 0; j < 3; j++) for (int k = 0; k < 3; k++) cen[k] += iface.points[3*tt[j] + k]/3.0;
    double d = std::sqrt((cen[0]-c[0])*(cen[0]-c[0]) + (cen[1]-c[1])*(cen[1]-c[1]) + (cen[2]-c[2])*(cen[2]-c[2]));
    far = std::max(far, d);
  }
  return far;
}


static double TetVolume(const std::vector<double> &pts, const ll *q)
{
  const double *a = &pts[3*q[0]], *b = &pts[3*q[1]], *c = &pts[3*q[2]], *d = &pts[3*q[3]];
  double e1[3], e2[3], e3[3];
  for (int k = 0; k < 3; k++) { e1[k] = b[k]-a[k]; e2[k] = c[k]-a[k]; e3[k] = d[k]-a[k]; }
  double cr[3] = {e2[1]*e3[2]-e2[2]*e3[1], e2[2]*e3[0]-e2[0]*e3[2], e2[0]*e3[1]-e2[1]*e3[0]};
  return (e1[0]*cr[0] + e1[1]*cr[1] + e1[2]*cr[2])/6.0;
}

// The faces of a set of tetrahedra by how many tetrahedra share them; the
// boundary is the faces on one.
static double TetVolume6(const std::vector<double> &p, const ll *q)
{
  double a[3], b[3], c[3];
  for (int k = 0; k < 3; k++) { a[k] = p[3*q[1]+k]-p[3*q[0]+k]; b[k] = p[3*q[2]+k]-p[3*q[0]+k]; c[k] = p[3*q[3]+k]-p[3*q[0]+k]; }
  return (a[1]*b[2]-a[2]*b[1])*c[0] + (a[2]*b[0]-a[0]*b[2])*c[1] + (a[0]*b[1]-a[1]*b[0])*c[2];
}
static void FaceCensus(const std::vector<ll> &tets, std::map<std::array<ll,3>, int> &count)
{
  const int faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}};
  for (size_t t = 0; t + 3 < tets.size(); t += 4)
    for (int f = 0; f < 4; f++)
    {
      std::array<ll,3> k = {tets[t + faces[f][0]], tets[t + faces[f][1]], tets[t + faces[f][2]]};
      std::sort(k.begin(), k.end());
      count[k]++;
    }
}

static bool IsBoundaryFace(const std::map<std::array<ll,3>, int> &count, const ll *tri)
{
  std::array<ll,3> k = {tri[0], tri[1], tri[2]};
  std::sort(k.begin(), k.end());
  std::map<std::array<ll,3>, int>::const_iterator it = count.find(k);
  return it != count.end() && it->second == 1;
}

// Prisms over the structured zone, with the checks every zone must pass.
static void CheckPrisms(const Interface &iface, const std::vector<unsigned char> &structured, int numLayers, PrismMesh &pm)
{
  std::string err;
  if (BuildPrismLayers(iface, structured, numLayers, pm, err) != 0) { printf("  FAIL prisms: %s\n", err.c_str()); numFailed++; return; }
  ll numStructured = 0; for (size_t t = 0; t < structured.size(); t++) if (structured[t]) numStructured++;
  ll nTets = (ll)(pm.tetrahedra.size()/4);
  double volume = 0.0, minVolume = 1e300;
  for (ll t = 0; t < nTets; t++) { double v = TetVolume(pm.points, &pm.tetrahedra[4*t]); volume += v; minVolume = std::min(minVolume, v); }
  std::map<std::array<ll,3>, int> count; FaceCensus(pm.tetrahedra, count);
  ll numBoundary = 0, numOver = 0;
  for (std::map<std::array<ll,3>, int>::const_iterator it = count.begin(); it != count.end(); ++it) { if (it->second == 1) numBoundary++; if (it->second > 2) numOver++; }
  ll topsOnBoundary = 0, sidesOnBoundary = 0, rimsOnBoundary = 0;
  for (size_t i = 0; i + 2 < pm.topTriangles.size(); i += 3) if (IsBoundaryFace(count, &pm.topTriangles[i])) topsOnBoundary++;
  for (size_t i = 0; i + 2 < pm.sideTriangles.size(); i += 3) if (IsBoundaryFace(count, &pm.sideTriangles[i])) sidesOnBoundary++;
  for (size_t i = 0; i + 2 < pm.rimTriangles.size(); i += 3) if (IsBoundaryFace(count, &pm.rimTriangles[i])) rimsOnBoundary++;
  ll nTop = (ll)(pm.topTriangles.size()/3), nSide = (ll)(pm.sideTriangles.size()/3), nRim = (ll)(pm.rimTriangles.size()/3), nEdges = (ll)(pm.zoneBoundaryEdges.size()/2);
  // the zone boundary edges chain into loops: every point once as a source and once as a target
  std::map<ll,int> outDeg, inDeg; for (ll e = 0; e < nEdges; e++) { outDeg[pm.zoneBoundaryEdges[2*e]]++; inDeg[pm.zoneBoundaryEdges[2*e+1]]++; }
  bool loops = true; for (std::map<ll,int>::iterator it = outDeg.begin(); it != outDeg.end(); ++it) if (it->second != 1 || inDeg[it->first] != 1) loops = false;
  printf("  prisms: %lld tetrahedra over %lld triangles (%lld points), volume %.4f, smallest %.2e; boundary faces %lld = bottom %lld + tops %lld + sides %lld + rim %lld; faces on more than two %lld; zone boundary edges %lld\n",
      nTets, numStructured, (ll)(pm.points.size()/3), volume, minVolume, numBoundary, numStructured, nTop, nSide, nRim, numOver, nEdges);
  Check(nTets == 3*numLayers*numStructured, "three tetrahedra per prism, N prisms per triangle");
  Check(minVolume > 0.0, "every tetrahedron has positive volume");
  Check(numOver == 0, "no face is shared by more than two tetrahedra");
  Check(numBoundary == numStructured + nTop + nSide + nRim, "the boundary is the bottom, the tops, the sides and the rim sides");
  Check(topsOnBoundary == nTop && sidesOnBoundary == nSide && rimsOnBoundary == nRim, "every top, side and rim triangle is a boundary face");
  Check(nSide == 2*numLayers*nEdges, "two side triangles per layer per zone boundary edge");
  Check(nEdges == 0 || loops, "the zone boundary edges chain into loops");
}

// The smallest dihedral angle of a tetrahedron, in degrees.
static double MinDihedral(const std::vector<double> &pts, const ll *q)
{
  const double *P[4]; for (int m = 0; m < 4; m++) P[m] = &pts[3*q[m]];
  const int faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{0,2,1}}; double N[4][3];
  for (int f = 0; f < 4; f++) { double u[3], v[3]; Sub(P[faces[f][1]], P[faces[f][0]], u); Sub(P[faces[f][2]], P[faces[f][0]], v); N[f][0] = u[1]*v[2]-u[2]*v[1]; N[f][1] = u[2]*v[0]-u[0]*v[2]; N[f][2] = u[0]*v[1]-u[1]*v[0]; double l = Norm(N[f]); if (l > 0) for (int k = 0; k < 3; k++) N[f][k] /= l; }
  double tetMin = 180.0;
  for (int f = 0; f < 4; f++) for (int g = f+1; g < 4; g++) { double c = std::max(-1.0, std::min(1.0, -Dot(N[f], N[g]))); tetMin = std::min(tetMin, std::acos(c)*180.0/M_PI); }
  return tetMin;
}

static void DihedralStats(const std::vector<double> &pts, const std::vector<ll> &tets, const char *label)
{
  ll n = (ll)(tets.size()/4), under5 = 0, under10 = 0, under15 = 0; double smallest = 180.0;
  for (ll t = 0; t < n; t++) { double d = MinDihedral(pts, &tets[4*t]); smallest = std::min(smallest, d); if (d < 5) under5++; if (d < 10) under10++; if (d < 15) under15++; }
  printf("  %s: %lld tetrahedra, smallest dihedral %.2f degrees, under 5: %lld (%.2f%%), under 10: %lld (%.2f%%), under 15: %lld (%.2f%%)\n", label, n, smallest, under5, n ? 100.0*under5/n : 0.0, under10, n ? 100.0*under10/n : 0.0, under15, n ? 100.0*under15/n : 0.0);
}

// The junction zone filled and joined to the prism layers: the whole wall as one tetrahedral mesh, over the given trimmed levels.
static void FillHybridWithLevels(const Interface &iface, int numLayers, const std::vector<Surface> &surfs, bool innerLevels = true, const BandRelaxation &relaxation = BandRelaxation())
{
  std::vector<double> fractions; for (int k = 1; k <= numLayers; k++) fractions.push_back((double)k/numLayers);
  std::string err;
  OffsetField field; Report fr; field.Build(iface, fr, err);
  std::vector<unsigned char> structured; ZoneReport zr; ZoneOptions zo;
  ClassifyPrismZone(iface, field, numLayers, zo, structured, zr, err);
  { long long nr = 0, nw = 0; double lw = 0.0; std::vector<CoveredRegion> regions;
    const double maxRadiusOverThickness = getenv("SV_HYBRID_WRAP_RADIUS") != nullptr ? atof(getenv("SV_HYBRID_WRAP_RADIUS")) : 12.0;
    const double wrapFactor = getenv("SV_HYBRID_WRAP_FACTOR") != nullptr ? atof(getenv("SV_HYBRID_WRAP_FACTOR")) : 0.0;
    ll w = WrapCoveredRegions(iface, field, zr.coveredPoints, zr.coveredBy, maxRadiusOverThickness, structured, nr, nw, lw, &regions, wrapFactor);
    if (nr > 0) printf("  %lld junction regions with a covered point, %lld wrapped around their vessels (the longest %.3g, %lld triangles)\n", nr, nw, lw, w);
    if (getenv("SV_HYBRID_ZIPS") != nullptr)
      for (size_t r = 0; r < regions.size(); r++)
        printf("    region %zu: %lld triangles, %lld points around (%.2f, %.2f, %.2f), %.0f%% covered (facing %.2f, %.0f%% under -0.8), radius %.3g (by curvature %.3g), thickness %.3g (r/t %.2f), against %lld%s%s\n", r, regions[r].numTriangles,
            regions[r].numPoints, regions[r].center[0], regions[r].center[1], regions[r].center[2], 100.0*regions[r].coveredFraction, regions[r].facing, 100.0*regions[r].facingFraction, regions[r].radius, regions[r].curvatureRadius, regions[r].thickness,
            regions[r].thickness > 0.0 ? regions[r].radius/regions[r].thickness : 0.0, regions[r].numContacts,
            regions[r].wrapped ? (std::string(regions[r].wrapped == 2 ? ", wrapped with one by " : ", wrapped by ") + std::to_string(regions[r].numWidened) + " triangles").c_str() : "",
            regions[r].wrapped ? (std::string("; ") + std::to_string(regions[r].numLoops) + " loops, the longest " + std::to_string(regions[r].longestLoop) + " circumferences, turning over 90 deg at " + std::to_string((int)(100.0*regions[r].loopTurning)) + "% of their points").c_str() : "");
    if (getenv("SV_HYBRID_ZIPS") != nullptr)
      for (size_t r = 0; r < regions.size(); r++)
      {
        const std::vector<double> &d = regions[r].facingDots;
        long long n = 0, u8 = 0, u5 = 0, u2 = 0;
        for (size_t m = 0; m + 3 < d.size(); m += 4) { n++; if (d[m] < -0.8) u8++; if (d[m] < -0.5) u5++; if (d[m] < -0.2) u2++; }
        if (n > 0) printf("      region %zu: of %lld covered points with a covering triangle, dot under -0.8: %.0f%%, under -0.5: %.0f%%, under -0.2: %.0f%%\n", r, n, 100.0*u8/n, 100.0*u5/n, 100.0*u2/n);
      }
    if (getenv("SV_HYBRID_REGION_DETAIL") != nullptr)
    {
      // the covered points of one region: a histogram of the dots and their positions, ten bins of dot, the mean position per bin
      const size_t r = (size_t)atoi(getenv("SV_HYBRID_REGION_DETAIL"));
      if (r < regions.size())
      {
        const std::vector<double> &d = regions[r].facingDots;
        for (int bin = 0; bin < 10; bin++)
        {
          const double lo = -1.0 + 0.2*bin, hi = lo + 0.2;
          long long n = 0; double c[3] = {0.0, 0.0, 0.0}, lo3[3] = {1e300, 1e300, 1e300}, hi3[3] = {-1e300, -1e300, -1e300};
          for (size_t m = 0; m + 3 < d.size(); m += 4)
          {
            if (d[m] < lo || d[m] >= hi) continue;
            n++;
            for (int k = 0; k < 3; k++) { c[k] += d[m + 1 + k]; lo3[k] = std::min(lo3[k], d[m + 1 + k]); hi3[k] = std::max(hi3[k], d[m + 1 + k]); }
          }
          if (n == 0) continue;
          printf("      dot [%.1f, %.1f): %lld points, mean (%.2f, %.2f, %.2f), x %.1f..%.1f y %.1f..%.1f z %.1f..%.1f\n", lo, hi, n, c[0]/n, c[1]/n, c[2]/n, lo3[0], hi3[0], lo3[1], hi3[1], lo3[2], hi3[2]);
        }
      }
    }
    if (getenv("SV_HYBRID_CLASSIFY_ONLY") != nullptr) return; }
  PrismMesh pm; JunctionShell shell;
  // a loop whose piece has no chain of its own widens its junction region by two rings, up to three times (as the glue does)
  int chainTakenInARow = 0;   // as the glue: the rings double while a loop's chain is taken by another loop (two regions on two vessels against each other)
  int crossingAttempts = 0;   // as the glue: the rings around the loops whose strips cross double with each such attempt
  for (int attempt = 0; ; attempt++)
  {
    pm = PrismMesh(); shell = JunctionShell();
    if (BuildPrismLayers(iface, structured, numLayers, pm, err) != 0) { printf("  FAIL prisms: %s\n", err.c_str()); numFailed++; return; }
    // with SV_HYBRID_SLOT_SIZES set: the slot smoothing's reach, for a run without it (0) or with another
    BandRelaxation relaxUsed = relaxation;
    if (getenv("SV_HYBRID_SLOT_SIZES") != nullptr) relaxUsed.slotSizes = atof(getenv("SV_HYBRID_SLOT_SIZES"));
    const bool built = BuildJunctionShell(iface, field, structured, pm, surfs, fractions, 1, 2, shell, err, innerLevels, relaxUsed) == 0;
    // as the glue does: a built shell with triangles passing through one another or folded onto each other widens the junction zone at the prisms whose layer
    // surface passes through a piece (a ring) and around the loops whose strips are at fault (two rings: a strip hops between the sides of a piece too narrow for it);
    // a loop without a chain widens its region; a piece with a chain no loop takes (a hole) gives the prisms owning or passing through it to the junction zone
    if (built && shell.numCrossingTriangles == 0 && shell.numFoldedEdges == 0) break;
    // with SV_HYBRID_ZIPS set: the strips at fault, the worst first, to see where and what they are
    if (built && getenv("SV_HYBRID_ZIPS") != nullptr)
    {
      std::vector<size_t> order;
      for (size_t z = 0; z < shell.zips.size(); z++) if (shell.zips[z].crossingTriangles > 0 || shell.zips[z].foldedEdges > 0) order.push_back(z);
      std::sort(order.begin(), order.end(), [&](size_t p, size_t q) { return shell.zips[p].crossingTriangles > shell.zips[q].crossingTriangles; });
      ll perLevel[8] = {0, 0, 0, 0, 0, 0, 0, 0}, loopsPerLevel[8] = {0, 0, 0, 0, 0, 0, 0, 0};
      for (size_t i = 0; i < order.size(); i++) { const ZipReport &z = shell.zips[order[i]]; if (z.level >= 0 && z.level < 8) { perLevel[z.level] += z.crossingTriangles; loopsPerLevel[z.level]++; } }
      printf("    strips at fault: %zu of %zu; crossings per level:", order.size(), shell.zips.size());
      for (int l = 0; l < 8; l++) if (loopsPerLevel[l] > 0) printf(" level %d: %lld in %lld strips;", l, perLevel[l], loopsPerLevel[l]);
      printf("\n");
      for (size_t i = 0; i < order.size() && i < 12; i++)
      {
        const ZipReport &z = shell.zips[order[i]];
        printf("    level %d loop %lld pts, chain %lld pts, strip %lld tris, %lld crossing, %lld folded, mean distance %.3f, agreement %.2f, around (%.2f, %.2f, %.2f)\n",
            z.level, z.loopPoints, z.chainPoints, z.stripTriangles, z.crossingTriangles, z.foldedEdges, z.meanDistance, z.agreement, z.centre[0], z.centre[1], z.centre[2]);
      }
    }
    // with SV_HYBRID_DUMP_ATTEMPT set: the first shell built with faults goes to that file (points, triangles, markers, prism point ids) and the test stops
    if (built && getenv("SV_HYBRID_DUMP_ATTEMPT") != nullptr && (getenv("SV_HYBRID_DUMP_ATTEMPT_AT") == nullptr || atoi(getenv("SV_HYBRID_DUMP_ATTEMPT_AT")) == attempt + 1))
    {
      FILE *f = fopen(getenv("SV_HYBRID_DUMP_ATTEMPT"), "wb");
      if (f != nullptr)
      {
        auto writeVec = [&](const void *data, size_t bytes, ll count) { fwrite(&count, sizeof(ll), 1, f); fwrite(data, 1, bytes, f); };
        writeVec(shell.points.data(), shell.points.size()*sizeof(double), (ll)shell.points.size());
        writeVec(shell.triangles.data(), shell.triangles.size()*sizeof(ll), (ll)shell.triangles.size());
        std::vector<ll> markers(shell.markers.begin(), shell.markers.end());
        writeVec(markers.data(), markers.size()*sizeof(ll), (ll)markers.size());
        writeVec(shell.prismPoint.data(), shell.prismPoint.size()*sizeof(ll), (ll)shell.prismPoint.size());
        fclose(f);
        printf("  the faulted shell of attempt %d written to %s%s\n", attempt + 1, getenv("SV_HYBRID_DUMP_ATTEMPT"), getenv("SV_HYBRID_DUMP_ATTEMPT_GO") != nullptr ? "" : "; stopping here");
      }
      if (getenv("SV_HYBRID_DUMP_ATTEMPT_GO") == nullptr) return;
    }
    std::vector<ll> seedPoints, seedTriangles; int rings = 1;
    std::vector<std::pair<std::vector<ll>, double> > wraps;   // loops wrapped round their vessels, each by its own distance
    // as the glue: a loop is wrapped round its vessel by the radius the lumen chords give (the curvature's is a quarter of it at a crotch),
    // unless the vessel is over twelve walls in radius (a band round the aorta would take most of it; rings there)
    auto wrapOrSeed = [&](const std::vector<ll> &loop, double fallbackRadius, bool wrap)
    {
      double radius = wrap ? ChordRadius(iface, field, loop) : 0.0;
      if (wrap && !(radius > 0.0)) radius = fallbackRadius;
      std::vector<double> th;
      for (size_t m = 0; m < loop.size(); m++) if (loop[m] >= 0 && (size_t)loop[m] < iface.thickness.size()) th.push_back(iface.thickness[(size_t)loop[m]]);
      std::sort(th.begin(), th.end());
      const double thickness = th.empty() ? 0.0 : th[th.size()/2];
      if (radius > 0.0 && radius <= 12.0*thickness) wraps.push_back(std::make_pair(loop, WrapDistance(radius)));
      else seedPoints.insert(seedPoints.end(), loop.begin(), loop.end());
    };
    if (built)
    {
      seedTriangles = shell.crossingStructuredTriangles;
      for (size_t z = 0; z < shell.zips.size(); z++)
      {
        const ZipReport &zip = shell.zips[z];
        if (zip.crossingTriangles == 0 && zip.foldedEdges == 0) continue;
        // a loop far from a third of its chain is on a vessel against another, and a loop whose strip is still at fault after a widening by
        // rings is at a crotch the rings do not mend (thin vessels on the refined model: strips over the short zone walls): wrapped
        wrapOrSeed(zip.loop, zip.loopRadius, zip.farFraction > 0.3 || crossingAttempts >= 1);
      }
      rings = std::min(16, 2 << crossingAttempts);
      crossingAttempts++;
    }
    else if (!shell.failedLoopPoints.empty())
    {
      chainTakenInARow = shell.failedLoopChainTaken ? chainTakenInARow + 1 : 0;
      wrapOrSeed(shell.failedLoopPoints, shell.failedLoopRadius, shell.failedLoopChainTaken);
      rings = wraps.empty() ? std::min(16, 2 << chainTakenInARow) : 2;
    }
    else if (!shell.failedTriangles.empty()) { seedTriangles = shell.failedTriangles; rings = 1; }
    const int maxAttempts = getenv("SV_HYBRID_ATTEMPTS") ? atoi(getenv("SV_HYBRID_ATTEMPTS")) : 16;
    if ((seedPoints.empty() && seedTriangles.empty() && wraps.empty()) || attempt >= maxAttempts) { if (built) break; printf("  FAIL shell: %s\n", err.c_str()); numFailed++; return; }
    ll widened = 0;
    const std::vector<unsigned char> zoneBefore = structured;
    if (!seedPoints.empty() || !seedTriangles.empty()) widened += WidenJunctionZone(iface, seedPoints, seedTriangles, rings, structured);
    for (size_t w = 0; w < wraps.size(); w++)
    {
      std::vector<ll> none; ll wrapped = WidenJunctionZone(iface, wraps[w].first, none, 0, structured, nullptr, wraps[w].second); widened += wrapped;
      printf("  attempt %d: a loop of %zu points wrapped around its vessel, %.3g along the surface (%lld triangles)\n", attempt + 1, wraps[w].first.size(), wraps[w].second, wrapped);
    }
    // as the glue: a widening that changed nothing would build the same shell again
    if (structured == zoneBefore) { printf("  attempt %d: the widening changed no triangle of the zone\n", attempt + 1); if (built) break; printf("  FAIL shell: %s\n", err.c_str()); numFailed++; return; }
    if (built) printf("  attempt %d: the shell has %lld crossing triangles and %lld folded edges; the junction zone widened by %lld triangles at %zu prisms and around %zu loop points (%d rings)\n", attempt + 1, shell.numCrossingTriangles, shell.numFoldedEdges, widened, seedTriangles.size(), seedPoints.size(), rings);
    else printf("  attempt %d: %s; the junction zone widened by %lld triangles (%d rings)\n", attempt + 1, err.c_str(), widened, rings);
  }
  ll nb, nn, nm; CountEdges(shell.triangles, nb, nn, nm);
  std::vector<unsigned char> cr; double at[3]; ll crossings = svenvelope::CountCrossingTriangles(shell.points, shell.triangles, cr, at);
  std::map<int, ll> byMarker; for (size_t i = 0; i < shell.markers.size(); i++) byMarker[shell.markers[i]]++;
  { ll thin = 0, total = 0; for (size_t z = 0; z < shell.zips.size(); z++) { thin += shell.zips[z].trianglesUnder15; total += shell.zips[z].stripTriangles; }
    printf("  zipper strips: %lld triangles, %lld with an angle under 15 degrees (%.1f%%)\n", total, thin, total > 0 ? 100.0*thin/total : 0.0); }
  printf("  band relaxed: %lld flips, %lld moves; under 10 degrees %lld -> %lld of %lld, smallest %.2f -> %.2f; relaxed shell had %lld crossing, %lld folded; %lld triangles undone locally%s\n", shell.numBandFlips, shell.numBandPointsMoved, shell.numBandUnder10Before, shell.numBandUnder10After, shell.numBandTriangles, shell.bandSmallestAngleBefore, shell.bandSmallestAngleAfter, shell.relaxationCrossingsAfter, shell.relaxationFoldsAfter, shell.numRelaxationUndone, shell.relaxationReverted ? " (reverted whole)" : "");
  printf("  creases: %lld band points held on crease edges; band triangles chording off the level %lld -> %lld; %lld slot points smoothed off the level\n", shell.numBandCreasePoints, shell.numBandChordsBefore, shell.numBandChordsAfter, shell.numBandSlotPoints);
  printf("  band edges collapsed %lld, free points on three triangles removed %lld (%lld triangles gone)%s\n", shell.numBandCollapses, shell.numBandVerticesRemoved, shell.numBandTrianglesRemoved, shell.collapseReverted ? " (undone: the shell crossed or folded)" : "");
  {
    // the band triangles still under 10 degrees, by how many of their corners are prism points (fixed by the relaxation): 0 piece interior, 1-2 strip or on the loop, 3 all fixed
    ll byFixed[4] = {0, 0, 0, 0}, worstTri = -1; double worst = 180.0; int worstFixed = 0;
    for (size_t i = 0; i < shell.markers.size(); i++)
    {
      if (shell.markers[i] == 1 || shell.markers[i] >= 300) continue;
      const ll *tri = &shell.triangles[3*i];
      double best = 180.0;
      for (int j = 0; j < 3; j++)
      {
        double ab[3], ac[3];
        for (int m = 0; m < 3; m++) { ab[m] = shell.points[3*tri[(j + 1)%3] + m] - shell.points[3*tri[j] + m]; ac[m] = shell.points[3*tri[(j + 2)%3] + m] - shell.points[3*tri[j] + m]; }
        const double la = std::sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]), lc = std::sqrt(ac[0]*ac[0] + ac[1]*ac[1] + ac[2]*ac[2]);
        if (!(la > 0.0) || !(lc > 0.0)) { best = 0.0; break; }
        best = std::min(best, std::acos(std::max(-1.0, std::min(1.0, (ab[0]*ac[0] + ab[1]*ac[1] + ab[2]*ac[2])/(la*lc))))*180.0/M_PI);
      }
      if (best >= 10.0) continue;
      int fixedCorners = 0;
      for (int j = 0; j < 3; j++) if (shell.prismPoint[(size_t)tri[j]] >= 0) fixedCorners++;
      byFixed[fixedCorners]++;
      if (best < worst) { worst = best; worstTri = (ll)i; worstFixed = fixedCorners; }
    }
    printf("  band triangles under 10 degrees by fixed corners: 0:%lld 1:%lld 2:%lld 3:%lld", byFixed[0], byFixed[1], byFixed[2], byFixed[3]);
    if (worstTri >= 0)
    {
      const ll *tri = &shell.triangles[3*worstTri];
      double c[3] = {0, 0, 0};
      for (int j = 0; j < 3; j++) for (int m = 0; m < 3; m++) c[m] += shell.points[3*tri[j] + m]/3.0;
      printf("; the worst %.2f degrees, marker %d, %d fixed corners, at (%.3f, %.3f, %.3f)", worst, shell.markers[worstTri], worstFixed, c[0], c[1], c[2]);
    }
    printf("\n");
  }
  printf("  junction shell: %lld points, %zu triangles (zone interface %lld, zipper %lld); by marker:", (ll)(shell.points.size()/3), shell.triangles.size()/3, shell.numJunctionTriangles, shell.numZipperTriangles);
  for (std::map<int, ll>::iterator it = byMarker.begin(); it != byMarker.end(); ++it) printf(" %d:%lld", it->first, it->second);
  printf("; boundary %lld, non-manifold %lld, miswound %lld, crossing %lld\n", nb, nn, nm, crossings);
  // the rings of the zone walls at the layers inside the wall are on three triangles: the wall below, the wall above and the layer's strip
  const ll ringEdges = innerLevels ? (ll)(numLayers - 1)*(ll)(pm.zoneBoundaryEdges.size()/2) : 0;   // without the inner layer surfaces the rings below the outer one lie on two wall triangles only
  Check(nb == 0 && nn == ringEdges && nm == 0 && crossings == 0, "the junction shell is closed, wound consistently and free of crossings (the layer rings on three triangles apart)");
  if (nb || nn != ringEdges || nm || crossings) return;
  // TetGen as the fill flow calls it
  tetgenio in, out; in.firstnumber = 0; in.numberofpoints = (int)(shell.points.size()/3); in.pointlist = new REAL[shell.points.size()]; for (size_t i = 0; i < shell.points.size(); i++) in.pointlist[i] = shell.points[i];
  in.numberoffacets = (int)(shell.triangles.size()/3); in.facetlist = new tetgenio::facet[in.numberoffacets]; in.facetmarkerlist = new int[in.numberoffacets]();
  for (int i = 0; i < in.numberoffacets; i++) { tetgenio::facet *f = &in.facetlist[i]; f->numberofpolygons = 1; f->polygonlist = new tetgenio::polygon[1]; f->numberofholes = 0; f->holelist = nullptr; tetgenio::polygon *pg = &f->polygonlist[0]; pg->numberofvertices = 3; pg->vertexlist = new int[3]; for (int j = 0; j < 3; j++) pg->vertexlist[j] = (int)shell.triangles[3*i+j]; in.facetmarkerlist[i] = shell.markers[i]; }
  tetgenbehavior b; b.plc = 1; b.nobisect = 1; b.nojettison = 1; b.quality = 1; b.minratio = 1.414; b.mindihedral = 10.0; b.neighout = 2; b.quiet = 1;
  bool accepted = true;
  try { tetrahedralize(&b, &in, &out); } catch (int r) { printf("  TetGen error %d\n", r); accepted = false; }
  Check(accepted, "TetGen fills the junction zone");
  if (!accepted) return;
  // the whole wall over one point array: the prism mesh's points, then the shell's own (the pieces'), then the zone's Steiner points
  std::vector<double> pts = pm.points; const ll nIn = in.numberofpoints;
  std::vector<ll> shellToAll(nIn, -1);
  for (ll i = 0; i < nIn; i++) { if (shell.prismPoint[i] >= 0) shellToAll[i] = shell.prismPoint[i]; else { shellToAll[i] = (ll)(pts.size()/3); for (int k = 0; k < 3; k++) pts.push_back(shell.points[3*i+k]); } }
  std::vector<ll> steinerToAll;
  for (int i = nIn; i < out.numberofpoints; i++) { steinerToAll.push_back((ll)(pts.size()/3)); for (int k = 0; k < 3; k++) pts.push_back(out.pointlist[3*i+k]); }
  std::vector<ll> zoneTets; for (int t = 0; t < out.numberoftetrahedra; t++) for (int m = 0; m < 4; m++) { ll v = out.tetrahedronlist[4*t+m]; zoneTets.push_back(v < nIn ? shellToAll[v] : steinerToAll[v - nIn]); }
  // the slivers of the fill eased as the fill flow does it: flips, then the free corners moved a little - the
  // pieces' own points and the mesher's own; the prism mesh's points (the interface's and the zone walls') not
  {
    std::vector<unsigned char> fixedPts(pts.size()/3, 0);
    for (size_t v = 0; v < pm.points.size()/3 && v < fixedPts.size(); v++) fixedPts[v] = 1;
    // the shell's triangles are the facets, and every shell point not on the outer surface alone stays
    std::vector<ll> facets;
    for (size_t i = 0; i + 2 < shell.triangles.size(); i += 3)
    {
      for (int j = 0; j < 3; j++) { facets.push_back(shellToAll[shell.triangles[i + j]]); if (shell.markers[i/3] != 2) fixedPts[(size_t)shellToAll[shell.triangles[i + j]]] = 1; }
    }
    double volumeBefore = 0.0; ll turnedBefore = 0;
    for (size_t t = 0; t + 3 < zoneTets.size(); t += 4) { const double v = TetVolume6(pts, &zoneTets[t]); volumeBefore += std::fabs(v); if (v < 0.0) turnedBefore++; }
    SliverReport sliver;
    ImproveSlivers(pts, zoneTets, fixedPts, facets, b.mindihedral, 5, 0.35, sliver);
    double volumeAfter = 0.0; ll turnedAfter = 0;
    for (size_t t = 0; t + 3 < zoneTets.size(); t += 4) { const double v = TetVolume6(pts, &zoneTets[t]); volumeAfter += std::fabs(v); if (v < 0.0) turnedAfter++; }
    printf("  slivers eased: %lld flips across a face, %lld around an edge (%lld tetrahedra -> %lld), %lld moves of %lld points; under %.0f degrees %lld -> %lld, under 1 degree %lld -> %lld, smallest %.3f -> %.3f; volume %.6g -> %.6g, turned %lld -> %lld\n",
        sliver.numFlips23, sliver.numFlips32, sliver.numTetrahedra, sliver.numTetrahedraAfter, sliver.numMoves, sliver.numPointsMoved, b.mindihedral, sliver.numUnderBefore, sliver.numUnderAfter, sliver.numUnder1Before, sliver.numUnder1After, sliver.smallestBefore, sliver.smallestAfter, volumeBefore/6, volumeAfter/6, turnedBefore, turnedAfter);
    Check(sliver.numUnderAfter < sliver.numUnderBefore, "easing the slivers leaves fewer of them");
    Check(turnedAfter == 0 && std::fabs(volumeAfter - volumeBefore) < 1.0e-3*volumeBefore, "easing the slivers turns no tetrahedron and keeps the zones' volume to a thousandth");
  }
  std::vector<ll> all = pm.tetrahedra; all.insert(all.end(), zoneTets.begin(), zoneTets.end());
  std::map<std::array<ll,3>, int> count; FaceCensus(all, count);
  ll numBoundary = 0, numOver = 0; for (std::map<std::array<ll,3>, int>::const_iterator it = count.begin(); it != count.end(); ++it) { if (it->second == 1) numBoundary++; if (it->second > 2) numOver++; }
  // the boundary of the wall: every interface triangle, the tops, the rim sides, the pieces and strips of the outer level
  ll interfaceOnBoundary = 0; for (size_t i = 0; i + 2 < iface.triangles.size(); i += 3) if (IsBoundaryFace(count, &iface.triangles[i])) interfaceOnBoundary++;
  ll topsOnBoundary = 0; for (size_t i = 0; i + 2 < pm.topTriangles.size(); i += 3) if (IsBoundaryFace(count, &pm.topTriangles[i])) topsOnBoundary++;
  ll rimOnBoundary = 0; for (size_t i = 0; i + 2 < pm.rimTriangles.size(); i += 3) if (IsBoundaryFace(count, &pm.rimTriangles[i])) rimOnBoundary++;
  ll outerShell = 0, outerOnBoundary = 0, wallsShared = 0, walls = 0, layersInside = 0, layersTotal = 0;
  for (size_t i = 0; i < shell.markers.size(); i++)
  {
    const ll *tri = &shell.triangles[3*i];
    std::array<ll,3> key = {shellToAll[tri[0]], shellToAll[tri[1]], shellToAll[tri[2]]}; std::sort(key.begin(), key.end());
    int c = count.count(key) ? count[key] : 0;
    if (shell.markers[i] == 2) { outerShell++; if (c == 1) outerOnBoundary++; }
    else if (shell.markers[i] >= 300) { walls++; if (c == 2) wallsShared++; }
    else if (shell.markers[i] >= 100) { layersTotal++; if (c == 2) layersInside++; }
  }
  double volume = 0.0; for (size_t t = 0; t + 3 < all.size(); t += 4) volume += TetVolume(pts, &all[t]);
  double area = 0.0, expected = 0.0;
  for (size_t i = 0; i + 2 < iface.triangles.size(); i += 3) { const ll *tt = &iface.triangles[i]; double e1[3], e2[3]; Sub(&iface.points[3*tt[1]], &iface.points[3*tt[0]], e1); Sub(&iface.points[3*tt[2]], &iface.points[3*tt[0]], e2); double crs[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]}; double a = 0.5*Norm(crs); area += a; expected += a*(iface.thickness[tt[0]] + iface.thickness[tt[1]] + iface.thickness[tt[2]])/3.0; }
  printf("  filled: prism %zu + zone %zu = %zu tetrahedra on %lld points (%d Steiner); boundary faces %lld: interface %lld/%zu, tops %lld/%zu, rim sides %lld/%zu, outer pieces and strips %lld/%lld; zone walls shared by two tetrahedra %lld/%lld, layer pieces and strips inside %lld/%lld; faces on more than two %lld; volume %.3f against area times thickness %.3f (%.3f)\n",
      pm.tetrahedra.size()/4, zoneTets.size()/4, all.size()/4, (ll)(pts.size()/3), out.numberofpoints - (int)nIn, numBoundary, interfaceOnBoundary, iface.triangles.size()/3, topsOnBoundary, pm.topTriangles.size()/3, rimOnBoundary, pm.rimTriangles.size()/3, outerOnBoundary, outerShell, wallsShared, walls, layersInside, layersTotal, numOver, volume, expected, volume/expected);
  Check(numOver == 0, "no face is shared by more than two tetrahedra");
  Check(interfaceOnBoundary == (ll)(iface.triangles.size()/3), "every interface triangle is a boundary face of the wall");
  Check(topsOnBoundary == (ll)(pm.topTriangles.size()/3) && rimOnBoundary == (ll)(pm.rimTriangles.size()/3) && outerOnBoundary == outerShell, "the tops, the rim sides and the outer pieces with their strips are the rest of the boundary");
  Check(numBoundary == interfaceOnBoundary + topsOnBoundary + rimOnBoundary + outerOnBoundary, "and nothing else is");
  Check(wallsShared == walls, "every zone wall triangle is shared by a prism tetrahedron and a zone tetrahedron");
  Check(layersInside == layersTotal, "every layer piece and strip triangle lies inside the wall between two tetrahedra");
  Check(volume > 1.0*expected && volume < 1.4*expected, "the wall's volume is the interface's area times the thickness, allowing for the curvature");
  DihedralStats(pts, pm.tetrahedra, "prism layers");
  DihedralStats(pts, zoneTets, "junction zone");
  DihedralStats(pts, all, "whole wall");
  {
    // the worst zone tetrahedron: where, and whose points (the prism mesh's, the pieces', or Steiner points of the fill)
    const ll numPrism = (ll)(pm.points.size()/3), numShellOwn = (ll)(pts.size()/3) - (ll)steinerToAll.size();
    ll worst = -1; double worstAngle = 180.0;
    for (size_t t = 0; t < zoneTets.size()/4; t++) { const double d = MinDihedral(pts, &zoneTets[4*t]); if (d < worstAngle) { worstAngle = d; worst = (ll)t; } }
    if (worst >= 0)
    {
      int prism = 0, piece = 0, steiner = 0; double c[3] = {0, 0, 0};
      for (int m = 0; m < 4; m++) { const ll v = zoneTets[4*worst + m]; if (v < numPrism) prism++; else if (v < numShellOwn) piece++; else steiner++; for (int k = 0; k < 3; k++) c[k] += pts[3*v + k]/4.0; }
      printf("  the worst zone tetrahedron: %.4f degrees at (%.3f, %.3f, %.3f), corners: %d prism, %d piece, %d Steiner\n", worstAngle, c[0], c[1], c[2], prism, piece, steiner);
    }
    // with SV_HYBRID_DUMP set, the shell and the zone's tetrahedra go to that file for a look afterwards
    if (const char *dump = getenv("SV_HYBRID_DUMP"))
    {
      FILE *f = fopen(dump, "wb");
      if (f != nullptr)
      {
        auto writeVec = [&](const void *data, size_t bytes, ll count) { fwrite(&count, sizeof(ll), 1, f); fwrite(data, 1, bytes, f); };
        writeVec(shell.points.data(), shell.points.size()*sizeof(double), (ll)shell.points.size());
        writeVec(shell.triangles.data(), shell.triangles.size()*sizeof(ll), (ll)shell.triangles.size());
        std::vector<ll> markers(shell.markers.begin(), shell.markers.end());
        writeVec(markers.data(), markers.size()*sizeof(ll), (ll)markers.size());
        writeVec(shell.prismPoint.data(), shell.prismPoint.size()*sizeof(ll), (ll)shell.prismPoint.size());
        writeVec(pts.data(), pts.size()*sizeof(double), (ll)pts.size());
        writeVec(zoneTets.data(), zoneTets.size()*sizeof(ll), (ll)zoneTets.size());
        fclose(f);
        printf("  the shell and the zone's tetrahedra written to %s\n", dump);
      }
    }
  }
}

// The same from the interface alone: the offset surfaces built and trimmed at the caps first.
static void FillHybrid(const Interface &iface, int numLayers)
{
  std::vector<double> fractions; for (int k = 1; k <= numLayers; k++) fractions.push_back((double)k/numLayers);
  std::vector<Surface> surfs; std::vector<Report> reps; std::string err; Options options;
  if (getenv("SV_CREASE_DEGREES") != nullptr) options.creaseDegrees = atof(getenv("SV_CREASE_DEGREES"));   // 0 leaves the valleys as marched, for comparison
  if (BuildOffsetSurfaces(iface, options, fractions, Delaunay, nullptr, surfs, reps, err) != 0) { printf("  FAIL build: %s\n", err.c_str()); numFailed++; return; }
  for (size_t k = 0; k < reps.size(); k++)
    printf("  level %zu valleys: %lld edges chorded one as marched, %lld split at the valley, %lld corners moved onto it, %lld left\n", k + 1, reps[k].numCreaseEdges, reps[k].numCreasePointsInserted, reps[k].numCreasePointsSnapped, reps[k].numCreaseCutsRefused);
  printf("  levels nested: %lld point moves pulling a level from under the valleys of the one over it\n", reps.back().numLevelPointsPulledIn);
  for (size_t k = 1; k < surfs.size(); k++)
  {
    std::vector<double> pts(surfs[k-1].points); std::vector<ll> tris(surfs[k-1].triangles);
    const ll base = (ll)(pts.size()/3);
    pts.insert(pts.end(), surfs[k].points.begin(), surfs[k].points.end());
    for (size_t i = 0; i < surfs[k].triangles.size(); i++) tris.push_back(surfs[k].triangles[i] + base);
    std::vector<unsigned char> cr; double at[3];
    const ll n = svenvelope::CountCrossingTriangles(pts, tris, cr, at);
    printf("  levels %zu and %zu together: %lld crossing triangles\n", k, k + 1, n);
    Check(n == 0, "consecutive levels do not cross each other");
  }
  // The valley where the branch meets its parent, cut into the outer level
  // (Options::creaseDegrees) against the contour as marched: triangles whose
  // centre stands off the level by over a tenth of their edge are chords
  // over the valley, the teeth along a branch root in the user's picture of
  // 2026-10-20.
  if (options.creaseDegrees > 0.0)
  {
    OffsetField field; Report fieldReport; std::string fieldError;
    if (field.Build(iface, fieldReport, fieldError, options.chordTolerance) == 0)
    {
      auto chords = [&](const Surface &surf) -> ll
      {
        ll n = 0;
        for (size_t i = 0; i + 2 < surf.triangles.size(); i += 3)
        {
          double c[3] = {0.0, 0.0, 0.0}, edge = 0.0;
          for (int j = 0; j < 3; j++)
          {
            const double *p = &surf.points[(size_t)3*surf.triangles[i + j]], *q = &surf.points[(size_t)3*surf.triangles[i + (j + 1)%3]];
            for (int k = 0; k < 3; k++) c[k] += p[k]/3.0;
            edge += Dist(p, q)/3.0;
          }
          double size, thickness; ll rim;
          field.Local(c, size, thickness, rim);
          if (rim >= 0) continue;   // the collars and their domes past the caps are trimmed away later
          if (std::abs(field.Evaluate(c, fractions.back())) > 0.1*edge) n++;
        }
        return n;
      };
      Options asMarched = options; asMarched.creaseDegrees = 0.0;
      std::vector<Surface> surfs0; std::vector<Report> reps0; std::string err0;
      if (BuildOffsetSurfaces(iface, asMarched, fractions, Delaunay, nullptr, surfs0, reps0, err0) == 0)
      {
        const ll before = chords(surfs0.back()), after = chords(surfs.back());
        printf("  valley chords on the outer level: %lld as marched, %lld with the valleys cut\n", before, after);
        Check(before > 0, "the contour as marched chords over the branch root's valley");
        Check(after*2 <= before, "the cut leaves at most half of the chords");
      }
    }
  }
  std::vector<CapPlane> planes = PlanesFromRims(iface, surfs.back());
  for (size_t k = 0; k < surfs.size(); k++) { TrimReport tr; if (TrimSurfaceAtCaps(surfs[k], planes, 0.1, tr, err) != 0) { printf("  FAIL trim at caps: %s\n", err.c_str()); numFailed++; return; } }
  FillHybridWithLevels(iface, numLayers, surfs);
}


// 9. The widening stops at the cap rims: seeds next to a rim grow the zone
// everywhere but onto the triangles with a point on the rim, and an island
// at the rim is kept.
static void TestWideningAtRims()
{
  printf("test 9: the zone widening leaves the band at the cap rims\n");
  Interface iface; AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3);
  const ll nt = (ll)(iface.triangles.size()/3), np = (ll)(iface.points.size()/3);
  // the rim points: on boundary edges
  std::map<std::pair<ll, ll>, int> edgeCount;
  for (ll t = 0; t < nt; t++) for (int j = 0; j < 3; j++) { ll a = iface.triangles[3*t + j], b = iface.triangles[3*t + (j + 1)%3]; edgeCount[std::make_pair(std::min(a, b), std::max(a, b))]++; }
  std::vector<unsigned char> rimPoint(np, 0);
  for (std::map<std::pair<ll, ll>, int>::const_iterator it = edgeCount.begin(); it != edgeCount.end(); ++it) if (it->second == 1) { rimPoint[it->first.first] = 1; rimPoint[it->first.second] = 1; }
  // seeds: the points of the third ring from the z = 0 rim (a loop that failed there)
  std::vector<ll> seeds;
  for (ll p = 0; p < np; p++) if (std::abs(iface.points[3*p + 2] - 0.75) < 0.2) seeds.push_back(p);
  std::vector<unsigned char> structured(nt, 1);
  // the seeds' own region: the triangles at the seed points are junction already
  for (ll t = 0; t < nt; t++) for (int j = 0; j < 3; j++) if (std::abs(iface.points[3*iface.triangles[3*t + j] + 2] - 0.75) < 0.2) structured[t] = 0;
  std::vector<ll> none;
  ll leftAtRims = 0;
  ll widened = WidenJunctionZone(iface, seeds, none, 16, structured, &leftAtRims);
  ll rimJunction = 0, rimTriangles = 0;
  for (ll t = 0; t < nt; t++)
  {
    bool atRim = false;
    for (int j = 0; j < 3; j++) if (rimPoint[iface.triangles[3*t + j]]) atRim = true;
    if (!atRim) continue;
    rimTriangles++;
    if (!structured[t]) rimJunction++;
  }
  printf("  %lld triangles widened with 16 rings, %lld left at the rims; %lld of the %lld rim triangles are junction\n", widened, leftAtRims, rimJunction, rimTriangles);
  Check(widened > 0, "the zone widened");
  Check(leftAtRims > 0, "triangles at the rim were held back");
  Check(rimJunction == 0, "no triangle with a point on a rim became junction");
  // seed triangles on the rim (a prism crossing the shell there, or owning a
  // hole) stay prisms too, and the pinches the widening opens leave the rim
  // band alone: seeded at the triangles of both rims and of a ring between
  std::vector<unsigned char> structured2(nt, 1);
  std::vector<ll> rimSeeds;
  for (ll t = 0; t < nt; t++)
  {
    bool atRim = false;
    for (int j = 0; j < 3; j++) if (rimPoint[iface.triangles[3*t + j]]) atRim = true;
    if (atRim || std::abs(iface.points[3*iface.triangles[3*t] + 2] - 4.0) < 0.2) rimSeeds.push_back(t);
  }
  std::vector<ll> noPoints;
  ll leftAtRims2 = 0;
  ll widened2 = WidenJunctionZone(iface, noPoints, rimSeeds, 2, structured2, &leftAtRims2);
  ll rimJunction2 = 0;
  for (ll t = 0; t < nt; t++)
  {
    bool atRim = false;
    for (int j = 0; j < 3; j++) if (rimPoint[iface.triangles[3*t + j]]) atRim = true;
    if (atRim && !structured2[t]) rimJunction2++;
  }
  printf("  seeded at %zu triangles on the rims and a ring between: %lld widened, %lld held back, %lld rim triangles junction\n", rimSeeds.size(), widened2, leftAtRims2, rimJunction2);
  Check(widened2 > 0 && leftAtRims2 > 0, "the ring widened and the rim seeds were held back");
  Check(rimJunction2 == 0, "no seed triangle on a rim became junction");
  // the boundary smoothing works on this widening's front alone: notches in
  // a zone far from the seeds (a band at z 5-7 with every fourth quad of the
  // row below it taken, whole quads so that no pinch is made) are left as
  // they are by a widening at z 2
  std::vector<unsigned char> structured3(nt, 1);
  auto zOf = [&](ll t, int j) { return iface.points[3*iface.triangles[3*t + j] + 2]; };
  ll inRow = 0;
  for (ll t = 0; t < nt; t++)
  {
    const double lo = std::min(zOf(t, 0), std::min(zOf(t, 1), zOf(t, 2))), hi = std::max(zOf(t, 0), std::max(zOf(t, 1), zOf(t, 2)));
    if (lo > 4.99 && hi < 7.01) structured3[t] = 0;
    else if (hi < 5.01 && hi > 4.99 && (inRow++ % 8) < 2) structured3[t] = 0;
  }
  std::map<std::pair<ll, ll>, std::vector<ll> > edgeTriangles;
  for (ll t = 0; t < nt; t++) for (int j = 0; j < 3; j++) { ll a = iface.triangles[3*t + j], b = iface.triangles[3*t + (j + 1)%3]; edgeTriangles[std::make_pair(std::min(a, b), std::max(a, b))].push_back(t); }
  std::vector<ll> notchesFar;
  for (ll t = 0; t < nt; t++)
  {
    if (!structured3[t]) continue;
    int junction = 0;
    for (int j = 0; j < 3; j++)
    {
      ll a = iface.triangles[3*t + j], b = iface.triangles[3*t + (j + 1)%3];
      const std::vector<ll> &on = edgeTriangles[std::make_pair(std::min(a, b), std::max(a, b))];
      for (size_t q = 0; q < on.size(); q++) if (on[q] != t && !structured3[on[q]]) junction++;
    }
    if (junction >= 2) notchesFar.push_back(t);
  }
  std::vector<ll> seeds3;
  for (ll p = 0; p < np; p++) if (std::abs(iface.points[3*p + 2] - 2.0) < 0.1) seeds3.push_back(p);
  for (ll t = 0; t < nt; t++) for (int j = 0; j < 3; j++) if (std::abs(zOf(t, j) - 2.0) < 0.1) structured3[t] = 0;
  WidenJunctionZone(iface, seeds3, none, 1, structured3);
  ll notchesTaken = 0;
  for (size_t m = 0; m < notchesFar.size(); m++) if (!structured3[notchesFar[m]]) notchesTaken++;
  printf("  %zu notches in a zone far from the seeds, %lld of them taken by the widening\n", notchesFar.size(), notchesTaken);
  Check(!notchesFar.empty() && notchesTaken == 0, "the smoothing left the zone away from the widening's front alone");
}

int main()
{
  {
    printf("test 1: a tube of radius 1 with a wall of 0.3, three layers - all of it prisms\n");
    Interface iface; AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3);
    std::vector<unsigned char> structured; ZoneReport zr;
    if (Classify(iface, 3, structured, zr))
    {
      Check(zr.numStructured == zr.numTriangles, "every triangle is structured");
      Check(zr.numPointsCovered == 0 && zr.numTrianglesInverted == 0 && zr.numTrianglesCrossing == 0, "nothing is covered, inverted or crossing");
    }
  }
  {
    printf("test 2: a thin branch (r 0.3, wall 0.1) on a thick parent (r 1, wall 0.5), three layers - the junction alone is tetrahedra\n");
    Interface iface; MakeJunction(0.5, 0.1, 0.3, 48, 20, iface);
    std::vector<unsigned char> structured; ZoneReport zr;
    if (Classify(iface, 3, structured, zr))
    {
      const double ostium[3] = {1.0, 0.0, 6.0};
      double far = FarthestJunction(iface, structured, ostium);
      printf("  the junction zone reaches %.2f from the ostium\n", far);
      Check(zr.numJunctionRegions == 1, "one junction region");
      Check(zr.numStructured > 0.6*zr.numTriangles, "most of the wall is prisms");
      Check(far < 2.0, "the junction zone stays within 2 of the ostium");
      Check(zr.numPointsCovered > 0, "the branch root's layer points are covered by the parent's wall");
    }
  }
  {
    printf("test 3: two tubes whose walls overlap between them (r 1 wall 0.3 and r 0.5 wall 0.1, surfaces 0.3 apart), one layer - a strip on each faces the other\n");
    Interface iface; AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3); AddTube(iface, 1.8, 0.0, 0.5, 0.0, 8.0, 20, 32, 0.1);
    std::vector<unsigned char> structured; ZoneReport zr;
    if (Classify(iface, 1, structured, zr))
    {
      // every junction triangle faces the other tube: its centre lies between the axes (x between 0 and 1.8;
      // the margin rings carry the zone around each tube, but not past its equator)
      ll off = 0, zone = 0;
      for (size_t t = 0; t < structured.size(); t++)
      {
        if (structured[t]) continue;
        zone++;
        const ll *tt = &iface.triangles[3*t];
        double x = (iface.points[3*tt[0]] + iface.points[3*tt[1]] + iface.points[3*tt[2]])/3.0;
        if (x < 0.0 || x > 1.8) off++;
      }
      printf("  junction zone %lld triangles, %lld of them not between the axes\n", zone, off);
      Check(zone > 0 && off == 0, "the junction zone is the strips facing the other tube");
      Check(zr.numStructured > 0.5*zr.numTriangles, "most of both tubes is prisms");
    }
    printf("test 3b: the same strips wrapped round their tubes before the shell (WrapCoveredRegions)\n");
    {
      OffsetField field; Report fr; std::string err;
      if (field.Build(iface, fr, err) != 0) { printf("  FAIL field: %s\n", err.c_str()); numFailed++; }
      else
      {
        long long nr = 0, nw = 0; double lw = 0.0; std::vector<CoveredRegion> regions;
        const ll w = WrapCoveredRegions(iface, field, zr.coveredPoints, zr.coveredBy, 12.0, structured, nr, nw, lw, &regions);
        printf("  %lld regions with a covered point, %lld wrapped (%lld triangles)\n", nr, nw, w);
        for (size_t r = 0; r < regions.size(); r++)
          printf("    region %zu around (%.2f, %.2f, %.2f): facing %.2f (%.0f%% under -0.8), radius %.3f by the chord (%.3f by curvature), against %lld, wrapped %d; %lld loops, the longest %.2f circumferences, turning over 90 deg at %.0f%%\n",
              r, regions[r].center[0], regions[r].center[1], regions[r].center[2], regions[r].facing, 100.0*regions[r].facingFraction, regions[r].radius, regions[r].curvatureRadius,
              regions[r].numContacts, (int)regions[r].wrapped, regions[r].numLoops, regions[r].longestLoop, 100.0*regions[r].loopTurning);
        Check(nr == 2 && nw == 2, "both strips face each other and are wrapped");
        bool radii = regions.size() == 2, facing = regions.size() == 2, bands = regions.size() == 2;
        for (size_t r = 0; r < regions.size(); r++)
        {
          const double expected = regions[r].center[0] < 0.9 ? 1.0 : 0.5;   // the tube of radius 1 at x 0, the one of radius 0.5 at x 1.8
          radii = radii && std::fabs(regions[r].radius - expected) < 0.1*expected;
          facing = facing && regions[r].facing < -0.8 && regions[r].numContacts == 1;
          bands = bands && regions[r].wrapped && regions[r].numLoops == 2 && regions[r].longestLoop > 0.8 && regions[r].longestLoop < 2.0 && regions[r].loopTurning < 0.1;
        }
        Check(radii, "the chord gives each tube's radius within a tenth");
        Check(facing, "the covering triangles face the covered points (dot under -0.8) and each strip is against the other");
        Check(bands, "each band is cut square round its tube: two end loops about a circumference long, few sharp turns");
        ll zone = 0;
        for (size_t t = 0; t < structured.size(); t++) if (!structured[t]) zone++;
        printf("  junction zone now %lld of %lld triangles\n", zone, (ll)structured.size());
      }
    }
  }
  {
    printf("test 4: prism layers over the whole tube (r 1, wall 0.3, three layers)\n");
    Interface iface; AddTube(iface, 0.0, 0.0, 1.0, 0.0, 8.0, 32, 32, 0.3);
    std::vector<unsigned char> structured(iface.triangles.size()/3, 1);
    PrismMesh pm; CheckPrisms(iface, structured, 3, pm);
    double volume = 0.0; for (size_t t = 0; t + 3 < pm.tetrahedra.size(); t += 4) volume += TetVolume(pm.points, &pm.tetrahedra[t]);
    double annulus = M_PI*(1.3*1.3 - 1.0)*8.0;
    printf("  volume %.4f against the annulus %.4f (%.3f)\n", volume, annulus, volume/annulus);
    Check(volume > 0.95*annulus && volume < 1.02*annulus, "the prisms fill the annulus between the tube and its offset");
    Check(pm.sideTriangles.empty() && pm.rimTriangles.size()/3 == (size_t)(2*3*64), "no zone walls, and the rim sides close the two ends");
  }
  {
    printf("test 5: prism layers over the structured zone of the junction, the junction left open toward its tetrahedra\n");
    Interface iface; MakeJunction(0.5, 0.1, 0.3, 48, 20, iface);
    std::vector<unsigned char> structured; ZoneReport zr;
    if (Classify(iface, 3, structured, zr))
    {
      PrismMesh pm; CheckPrisms(iface, structured, 3, pm);
    }
  }
  {
    printf("test 6: the three layer surfaces of the junction, trimmed to the junction zone and zipped to the prism zone's layers\n");
    Interface iface; MakeJunction(0.5, 0.1, 0.3, 48, 20, iface);
    std::vector<double> fractions = {1.0/3.0, 2.0/3.0, 1.0};
    std::vector<Surface> surfs; std::vector<Report> reps; std::string err; Options options;
    if (BuildOffsetSurfaces(iface, options, fractions, Delaunay, nullptr, surfs, reps, err) != 0) { printf("  FAIL build: %s\n", err.c_str()); numFailed++; }
    else
    {
      std::vector<CapPlane> planes = PlanesFromRims(iface, surfs.back());
      for (size_t k = 0; k < surfs.size(); k++) { TrimReport tr; if (TrimSurfaceAtCaps(surfs[k], planes, 0.1, tr, err) != 0) { printf("  FAIL trim at caps: %s\n", err.c_str()); numFailed++; } }
      OffsetField field; Report fr; field.Build(iface, fr, err);
      std::vector<unsigned char> structured; ZoneReport zr; ZoneOptions zo;
      ClassifyPrismZone(iface, field, 3, zo, structured, zr, err);
      PrismMesh pm; BuildPrismLayers(iface, structured, 3, pm, err);
      for (int k = 1; k <= 3; k++) CheckLevelPiece(iface, field, surfs[(size_t)k-1], fractions[(size_t)k-1], k, structured, pm);
    }
  }
  {
    printf("test 7: the whole wall of the junction as prism layers plus the junction zone's tetrahedra, three layers\n");
    Interface iface; MakeJunction(0.5, 0.1, 0.3, 48, 20, iface);
    FillHybrid(iface, 3);
  }
  {
    printf("test 8: the same with one layer\n");
    Interface iface; MakeJunction(0.5, 0.1, 0.3, 48, 20, iface);
    FillHybrid(iface, 1);
  }
  TestWideningAtRims();
  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
