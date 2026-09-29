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
// Build and run:
//   g++ -O2 -std=c++17 -I. Testing/test_hybrid.cxx sv_tetgenmesh_offset.cxx sv_tetgenmesh_envelope.cxx -o test_hybrid && ./test_hybrid
#include "sv_tetgenmesh_offset.h"
#include "sv_tetgenmesh_envelope.h"
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
      // every junction triangle faces the other tube: its centre lies between the axes (x between 0 and 1.8)
      ll off = 0, zone = 0;
      for (size_t t = 0; t < structured.size(); t++)
      {
        if (structured[t]) continue;
        zone++;
        const ll *tt = &iface.triangles[3*t];
        double x = (iface.points[3*tt[0]] + iface.points[3*tt[1]] + iface.points[3*tt[2]])/3.0;
        if (x < 0.3 || x > 1.6) off++;
      }
      printf("  junction zone %lld triangles, %lld of them not between the axes\n", zone, off);
      Check(zone > 0 && off == 0, "the junction zone is the strips facing the other tube");
      Check(zr.numStructured > 0.5*zr.numTriangles, "most of both tubes is prisms");
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
  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
