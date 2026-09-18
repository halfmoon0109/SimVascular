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


/** @file sv_tetgenmesh_envelope.cxx
 *  @brief The outer envelope of a closed, self-intersecting triangle surface.
 *  @details See sv_tetgenmesh_envelope.h for what is computed and why. The
 *  file has no dependency beyond the standard library so that it can be built
 *  and run on its own against synthetic surfaces.
 */

#include "sv_tetgenmesh_envelope.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <iterator>
#include <map>
#include <utility>
#include <vector>

namespace svenvelope {

namespace {

typedef long long ll;

//---------------------
// Vector arithmetic
//---------------------

inline void Sub(const double a[3], const double b[3], double r[3])
{
  r[0] = a[0] - b[0];
  r[1] = a[1] - b[1];
  r[2] = a[2] - b[2];
}

inline void Cross(const double a[3], const double b[3], double r[3])
{
  r[0] = a[1]*b[2] - a[2]*b[1];
  r[1] = a[2]*b[0] - a[0]*b[2];
  r[2] = a[0]*b[1] - a[1]*b[0];
}

inline double Dot(const double a[3], const double b[3])
{
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline double Norm(const double a[3])
{
  return std::sqrt(Dot(a, a));
}

inline double Distance(const double a[3], const double b[3])
{
  double d[3];
  Sub(a, b, d);
  return Norm(d);
}

inline bool Normalize(double a[3])
{
  double n = Norm(a);
  if (!(n > 0.0))
  {
    return false;
  }
  a[0] /= n;
  a[1] /= n;
  a[2] /= n;
  return true;
}

//---------------------
// Keys
//---------------------

struct EdgeKey
{
  ll a, b;
  EdgeKey(ll p, ll q) : a(std::min(p, q)), b(std::max(p, q)) {}
  bool operator<(const EdgeKey &o) const
  {
    return (a != o.a) ? (a < o.a) : (b < o.b);
  }
};

struct EdgeTriangleKey
{
  ll a, b, tri;
  EdgeTriangleKey(ll p, ll q, ll t) : a(std::min(p, q)), b(std::max(p, q)), tri(t) {}
  bool operator<(const EdgeTriangleKey &o) const
  {
    if (a != o.a) return a < o.a;
    if (b != o.b) return b < o.b;
    return tri < o.tri;
  }
};

struct TripleKey
{
  ll a, b, c;
  TripleKey(ll p, ll q, ll r)
  {
    ll s[3] = {p, q, r};
    std::sort(s, s + 3);
    a = s[0];
    b = s[1];
    c = s[2];
  }
  bool operator<(const TripleKey &o) const
  {
    if (a != o.a) return a < o.a;
    if (b != o.b) return b < o.b;
    return c < o.c;
  }
};

//---------------------
// Grid
//---------------------
// Triangles binned on a uniform grid by their bounding boxes, so that the
// triangles near a point, a box or a ray are the ones in the bins it covers.

struct Grid
{
  double origin[3] = {0.0, 0.0, 0.0};
  double size = 1.0;
  int n[3] = {1, 1, 1};
  std::vector<ll> start;
  std::vector<ll> cells;

  int Bin(double x, int axis) const
  {
    double f = std::floor((x - origin[axis])/size);
    if (f < 0.0) return 0;
    if (f >= (double)n[axis]) return n[axis] - 1;
    return (int)f;
  }

  ll Index(int i, int j, int k) const
  {
    return ((ll)k*n[1] + j)*n[0] + i;
  }

  // boxes: six values per triangle. typicalExtent: the size a bin is a
  // couple of. maxBins: the bin count the grid is kept under, at the cost of
  // more triangles per bin.
  void Build(const std::vector<double> &boxes, ll numTris, double typicalExtent, ll maxBins)
  {
    double lo[3], hi[3];
    for (int k = 0; k < 3; k++)
    {
      lo[k] = std::numeric_limits<double>::max();
      hi[k] = -std::numeric_limits<double>::max();
    }
    for (ll t = 0; t < numTris; t++)
    {
      for (int k = 0; k < 3; k++)
      {
        lo[k] = std::min(lo[k], boxes[(size_t)6*t + 2*k]);
        hi[k] = std::max(hi[k], boxes[(size_t)6*t + 2*k + 1]);
      }
    }
    if (numTris == 0)
    {
      for (int k = 0; k < 3; k++)
      {
        lo[k] = 0.0;
        hi[k] = 1.0;
      }
    }
    double extent[3];
    double volume = 1.0;
    for (int k = 0; k < 3; k++)
    {
      extent[k] = std::max(hi[k] - lo[k], 1.0e-12);
      volume *= extent[k];
    }
    size = std::max(2.0*typicalExtent, 1.0e-12);
    double floorSize = std::cbrt(volume/(double)maxBins);
    size = std::max(size, floorSize);
    // A little slack so that nothing sits on the far boundary.
    for (int k = 0; k < 3; k++)
    {
      origin[k] = lo[k] - 0.5*size;
      double count = std::ceil((extent[k] + size)/size);
      n[k] = std::max(1, (int)count);
    }
    ll numBins = (ll)n[0]*n[1]*n[2];
    start.assign((size_t)numBins + 1, 0);
    for (ll t = 0; t < numTris; t++)
    {
      const double *b = &boxes[(size_t)6*t];
      int i0 = Bin(b[0], 0), i1 = Bin(b[1], 0);
      int j0 = Bin(b[2], 1), j1 = Bin(b[3], 1);
      int k0 = Bin(b[4], 2), k1 = Bin(b[5], 2);
      for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++)
          for (int i = i0; i <= i1; i++)
            start[(size_t)Index(i, j, k) + 1]++;
    }
    for (ll b = 0; b < numBins; b++)
    {
      start[(size_t)b + 1] += start[(size_t)b];
    }
    cells.assign((size_t)start[(size_t)numBins], -1);
    std::vector<ll> cursor(start.begin(), start.end() - 1);
    for (ll t = 0; t < numTris; t++)
    {
      const double *b = &boxes[(size_t)6*t];
      int i0 = Bin(b[0], 0), i1 = Bin(b[1], 0);
      int j0 = Bin(b[2], 1), j1 = Bin(b[3], 1);
      int k0 = Bin(b[4], 2), k1 = Bin(b[5], 2);
      for (int k = k0; k <= k1; k++)
        for (int j = j0; j <= j1; j++)
          for (int i = i0; i <= i1; i++)
            cells[(size_t)cursor[(size_t)Index(i, j, k)]++] = t;
    }
  }
};

inline bool BoxesOverlap(const double *a, const double *b)
{
  for (int k = 0; k < 3; k++)
  {
    if (a[2*k] > b[2*k+1] || b[2*k] > a[2*k+1])
    {
      return false;
    }
  }
  return true;
}

//---------------------
// Triangle geometry
//---------------------

struct TriangleGeometry
{
  std::vector<double> normal;   // three per triangle, unit; zero for a degenerate triangle
  std::vector<double> box;      // six per triangle
  std::vector<double> extent;   // largest box side per triangle
  std::vector<unsigned char> degenerate;

  void Build(const std::vector<double> &points, const std::vector<ll> &tris, ll numTris)
  {
    normal.assign((size_t)3*numTris, 0.0);
    box.assign((size_t)6*numTris, 0.0);
    extent.assign((size_t)numTris, 0.0);
    degenerate.assign((size_t)numTris, 0);
    for (ll t = 0; t < numTris; t++)
    {
      const double *p0 = &points[(size_t)3*tris[(size_t)3*t]];
      const double *p1 = &points[(size_t)3*tris[(size_t)3*t + 1]];
      const double *p2 = &points[(size_t)3*tris[(size_t)3*t + 2]];
      double e1[3], e2[3], n[3];
      Sub(p1, p0, e1);
      Sub(p2, p0, e2);
      Cross(e1, e2, n);
      double area2 = Norm(n);
      double longest2 = std::max(Dot(e1, e1), Dot(e2, e2));
      if (!(area2 > 1.0e-14*longest2))
      {
        degenerate[(size_t)t] = 1;
      }
      else
      {
        for (int k = 0; k < 3; k++)
        {
          normal[(size_t)3*t + k] = n[k]/area2;
        }
      }
      double *b = &box[(size_t)6*t];
      for (int k = 0; k < 3; k++)
      {
        b[2*k] = std::min(p0[k], std::min(p1[k], p2[k]));
        b[2*k+1] = std::max(p0[k], std::max(p1[k], p2[k]));
        extent[(size_t)t] = std::max(extent[(size_t)t], b[2*k+1] - b[2*k]);
      }
    }
  }
};

// Barycentric weights of a point against a triangle, in the triangle's own
// plane (the point is projected along the normal by the arithmetic).
inline void Barycentric(const double p[3], const double a[3], const double b[3],
    const double c[3], double w[3])
{
  double v0[3], v1[3], v2[3];
  Sub(b, a, v0);
  Sub(c, a, v1);
  Sub(p, a, v2);
  double d00 = Dot(v0, v0), d01 = Dot(v0, v1), d11 = Dot(v1, v1);
  double d20 = Dot(v2, v0), d21 = Dot(v2, v1);
  double denom = d00*d11 - d01*d01;
  if (!(std::abs(denom) > 0.0))
  {
    w[0] = w[1] = w[2] = -1.0;
    return;
  }
  w[1] = (d11*d20 - d01*d21)/denom;
  w[2] = (d00*d21 - d01*d20)/denom;
  w[0] = 1.0 - w[1] - w[2];
}

// One end of the segment two triangles cross along.
struct CrossingEnd
{
  int kind;        // 0 a shared corner, 1 on an edge of the first triangle, 2 on an edge of the second
  ll v0, v1;       // the edge, v0 < v1; for a corner both are the corner
  double u;        // along the edge from v0
  double p[3];
};

// The ends of the segment along which two triangles cross, as the places
// where an edge of either passes through the other. Two triangles that share
// a corner meet there by construction; the corner is one end, and only the
// edges opposite it are tried for the other. Two that share an edge are not
// crossed (they can only overlap, coplanar) and give nothing.
// Returns the number of ends found, at most seven; the caller picks the pair.
int CrossingEnds(const std::vector<double> &points, const ll ta[3], const double na[3],
    const ll tb[3], const double nb[3], double baryTol, CrossingEnd ends[7], int &numShared)
{
  int found = 0;
  numShared = 0;
  ll shared = -1;
  for (int m = 0; m < 3; m++)
  {
    for (int n = 0; n < 3; n++)
    {
      if (ta[m] == tb[n])
      {
        numShared++;
        shared = ta[m];
      }
    }
  }
  if (numShared >= 2)
  {
    return 0;
  }
  if (numShared == 1)
  {
    CrossingEnd &e = ends[found++];
    e.kind = 0;
    e.v0 = e.v1 = shared;
    e.u = 0.0;
    for (int k = 0; k < 3; k++)
    {
      e.p[k] = points[(size_t)3*shared + k];
    }
  }
  for (int side = 0; side < 2; side++)
  {
    const ll *edgesOf = (side == 0) ? ta : tb;
    const ll *other = (side == 0) ? tb : ta;
    const double *otherNormal = (side == 0) ? nb : na;
    const double *o0 = &points[(size_t)3*other[0]];
    const double *o1 = &points[(size_t)3*other[1]];
    const double *o2 = &points[(size_t)3*other[2]];
    for (int j = 0; j < 3; j++)
    {
      ll p = edgesOf[j], q = edgesOf[(j+1)%3];
      if (p == shared || q == shared)
      {
        continue;
      }
      ll v0 = std::min(p, q), v1 = std::max(p, q);
      const double *a = &points[(size_t)3*v0];
      const double *b = &points[(size_t)3*v1];
      double ra[3], rb[3];
      Sub(a, o0, ra);
      Sub(b, o0, rb);
      double da = Dot(otherNormal, ra);
      double db = Dot(otherNormal, rb);
      if (!((da > 0.0 && db < 0.0) || (da < 0.0 && db > 0.0)))
      {
        continue;
      }
      double u = da/(da - db);
      double x[3];
      for (int k = 0; k < 3; k++)
      {
        x[k] = a[k] + u*(b[k] - a[k]);
      }
      double w[3];
      Barycentric(x, o0, o1, o2, w);
      if (w[0] < -baryTol || w[1] < -baryTol || w[2] < -baryTol)
      {
        continue;
      }
      CrossingEnd &e = ends[found++];
      e.kind = 1 + side;
      e.v0 = v0;
      e.v1 = v1;
      e.u = u;
      for (int k = 0; k < 3; k++)
      {
        e.p[k] = x[k];
      }
    }
  }
  return found;
}

// Picks the two ends of the crossing segment out of the candidates: with a
// shared corner, the corner and the candidate nearest it; otherwise the two
// farthest apart, since a candidate doubled near an edge of the other
// triangle lies between the true ends. Returns false when there is no
// segment of any length.
bool PickCrossing(const CrossingEnd ends[7], int found, int numShared, double minLength,
    CrossingEnd out[2])
{
  if (numShared == 1)
  {
    if (found < 2)
    {
      return false;
    }
    int best = -1;
    double bestDistance = std::numeric_limits<double>::max();
    for (int e = 1; e < found; e++)
    {
      double d = Distance(ends[0].p, ends[e].p);
      if (d < bestDistance)
      {
        bestDistance = d;
        best = e;
      }
    }
    if (best < 0 || bestDistance <= minLength)
    {
      return false;
    }
    out[0] = ends[0];
    out[1] = ends[best];
    return true;
  }
  if (found < 2)
  {
    return false;
  }
  int bi = -1, bj = -1;
  double bestDistance = -1.0;
  for (int i = 0; i < found; i++)
  {
    for (int j = i + 1; j < found; j++)
    {
      double d = Distance(ends[i].p, ends[j].p);
      if (d > bestDistance)
      {
        bestDistance = d;
        bi = i;
        bj = j;
      }
    }
  }
  if (bestDistance <= minLength)
  {
    return false;
  }
  out[0] = ends[bi];
  out[1] = ends[bj];
  return true;
}

//---------------------
// Point registry
//---------------------
// The crease points, keyed so that every triangle that has one gets the same
// id: a point on an edge by the edge and the triangle crossing it, a point
// where two creases cross inside a triangle by the three triangles meeting
// there. Points closer than a snap tolerance along an edge are one point, and
// a point within it of an edge's end is that end.

struct PointRegistry
{
  std::vector<double> &points;
  std::vector<int> &kind;
  std::map<EdgeTriangleKey, ll> edgeTriangle;
  std::map<EdgeKey, std::vector<std::pair<double, ll> > > onEdge;
  std::map<TripleKey, ll> triple;
  double snapU;
  double coincidence;
  double cell;
  std::map<std::pair<std::pair<ll, ll>, ll>, std::vector<ll> > hash;
  ll numInput;

  PointRegistry(std::vector<double> &pts, std::vector<int> &knd, double snap, double coincide)
    : points(pts), kind(knd), snapU(snap), coincidence(coincide), cell(std::max(4.0*coincide, 1.0e-300))
  {
    numInput = (ll)(points.size()/3);
    for (ll id = 0; id < numInput; id++)
    {
      hash[Cell(&points[(size_t)3*id])].push_back(id);
    }
  }

  std::pair<std::pair<ll, ll>, ll> Cell(const double p[3]) const
  {
    return std::make_pair(std::make_pair((ll)std::floor(p[0]/cell), (ll)std::floor(p[1]/cell)),
        (ll)std::floor(p[2]/cell));
  }

  // The point already registered within the coincidence distance of a
  // position, or -1. Two crossing points that land on one another - an edge
  // of one triangle through an edge of another, a crease through a corner,
  // two creases through the same point - are one point, whichever way each
  // was found, and this is what makes them so.
  ll Near(const double p[3])
  {
    std::pair<std::pair<ll, ll>, ll> c = Cell(p);
    ll best = -1;
    double bestDistance = coincidence;
    for (ll dk = -1; dk <= 1; dk++)
    {
      for (ll dj = -1; dj <= 1; dj++)
      {
        for (ll di = -1; di <= 1; di++)
        {
          std::pair<std::pair<ll, ll>, ll> key = std::make_pair(
              std::make_pair(c.first.first + di, c.first.second + dj), c.second + dk);
          std::map<std::pair<std::pair<ll, ll>, ll>, std::vector<ll> >::iterator found = hash.find(key);
          if (found == hash.end())
          {
            continue;
          }
          for (size_t m = 0; m < found->second.size(); m++)
          {
            double d = Distance(p, &points[(size_t)3*found->second[m]]);
            if (d <= bestDistance)
            {
              bestDistance = d;
              best = found->second[m];
            }
          }
        }
      }
    }
    return best;
  }

  ll FindOrNew(const double p[3], int k)
  {
    ll near = Near(p);
    if (near >= 0)
    {
      return near;
    }
    ll id = (ll)(points.size()/3);
    points.push_back(p[0]);
    points.push_back(p[1]);
    points.push_back(p[2]);
    kind.push_back(k);
    hash[Cell(p)].push_back(id);
    return id;
  }

  ll EdgePoint(ll v0, ll v1, ll tri, double u, const double p[3])
  {
    if (u < snapU)
    {
      return v0;
    }
    if (u > 1.0 - snapU)
    {
      return v1;
    }
    EdgeTriangleKey key(v0, v1, tri);
    std::map<EdgeTriangleKey, ll>::iterator found = edgeTriangle.find(key);
    if (found != edgeTriangle.end())
    {
      return found->second;
    }
    std::vector<std::pair<double, ll> > &list = onEdge[EdgeKey(v0, v1)];
    for (size_t m = 0; m < list.size(); m++)
    {
      if (std::abs(list[m].first - u) < snapU)
      {
        edgeTriangle[key] = list[m].second;
        return list[m].second;
      }
    }
    ll id = FindOrNew(p, 1);
    list.push_back(std::make_pair(u, id));
    edgeTriangle[key] = id;
    return id;
  }

  ll TriplePoint(ll a, ll b, ll c, const double p[3], bool &created)
  {
    TripleKey key(a, b, c);
    std::map<TripleKey, ll>::iterator found = triple.find(key);
    created = false;
    if (found != triple.end())
    {
      return found->second;
    }
    ll before = (ll)(points.size()/3);
    ll id = FindOrNew(p, 2);
    created = id == before;
    triple[key] = id;
    return id;
  }
};

//---------------------
// Ray casting
//---------------------

struct RayHit
{
  int status;   // 0 miss, 1 hit, 2 grazed an edge - the ray is not to be trusted
  int sign;     // +1 leaving through the front of the triangle, -1 entering through its back
  double t;
};

inline RayHit RayTriangle(const double o[3], const double d[3], const double a[3],
    const double b[3], const double c[3], double edgeTol)
{
  RayHit hit = {0, 0, 0.0};
  double e1[3], e2[3], pvec[3], tvec[3], qvec[3];
  Sub(b, a, e1);
  Sub(c, a, e2);
  Cross(d, e2, pvec);
  double det = Dot(e1, pvec);
  double scale = std::sqrt(Dot(e1, e1)*Dot(e2, e2));
  if (!(std::abs(det) > 1.0e-14*scale))
  {
    return hit;
  }
  double inv = 1.0/det;
  Sub(o, a, tvec);
  double u = Dot(tvec, pvec)*inv;
  if (u < -edgeTol || u > 1.0 + edgeTol)
  {
    return hit;
  }
  Cross(tvec, e1, qvec);
  double v = Dot(d, qvec)*inv;
  if (v < -edgeTol || u + v > 1.0 + edgeTol)
  {
    return hit;
  }
  double t = Dot(e2, qvec)*inv;
  if (!(t > 0.0))
  {
    return hit;
  }
  if (u < edgeTol || v < edgeTol || u + v > 1.0 - edgeTol)
  {
    hit.status = 2;
    hit.t = t;
    return hit;
  }
  hit.status = 1;
  // det = e1 . (d x e2) = -d . (e1 x e2): negative when the ray runs along
  // the triangle's normal, which is leaving the solid through its front.
  hit.sign = (det < 0.0) ? 1 : -1;
  hit.t = t;
  return hit;
}

// A small deterministic pseudo-random sequence for perturbing ray directions.
struct Scramble
{
  unsigned long long state;
  explicit Scramble(unsigned long long seed) : state(seed*6364136223846793005ULL + 1442695040888963407ULL) {}
  double Next()
  {
    state = state*6364136223846793005ULL + 1442695040888963407ULL;
    return ((double)(state >> 11))/9007199254740992.0*2.0 - 1.0;
  }
};

// The winding number of the closed surface just off a point in a direction:
// the signed count of the triangles a ray from the point crosses, leaving out
// one triangle (the one the point is on). Returns false when the ray grazed
// an edge and the count cannot be trusted.
bool WindingAlongRay(const double origin[3], const double direction[3], ll skipTriangle,
    const std::vector<double> &points, const std::vector<ll> &tris,
    const TriangleGeometry &geometry, const Grid &grid, std::vector<ll> &stamp, ll &stampValue,
    double edgeTol, int &winding, const std::vector<unsigned char> *skipSet = nullptr)
{
  winding = 0;
  stampValue++;
  // Amanatides & Woo over the bins the ray passes through.
  int cell[3], step[3];
  double tMax[3], tDelta[3];
  for (int k = 0; k < 3; k++)
  {
    cell[k] = grid.Bin(origin[k], k);
    if (direction[k] > 0.0)
    {
      step[k] = 1;
      double boundary = grid.origin[k] + (double)(cell[k] + 1)*grid.size;
      tMax[k] = (boundary - origin[k])/direction[k];
      tDelta[k] = grid.size/direction[k];
    }
    else if (direction[k] < 0.0)
    {
      step[k] = -1;
      double boundary = grid.origin[k] + (double)cell[k]*grid.size;
      tMax[k] = (boundary - origin[k])/direction[k];
      tDelta[k] = -grid.size/direction[k];
    }
    else
    {
      step[k] = 0;
      tMax[k] = std::numeric_limits<double>::max();
      tDelta[k] = std::numeric_limits<double>::max();
    }
  }
  while (true)
  {
    ll bin = grid.Index(cell[0], cell[1], cell[2]);
    for (ll c = grid.start[(size_t)bin]; c < grid.start[(size_t)bin + 1]; c++)
    {
      ll t = grid.cells[(size_t)c];
      if (t == skipTriangle || stamp[(size_t)t] == stampValue)
      {
        continue;
      }
      stamp[(size_t)t] = stampValue;
      if (geometry.degenerate[(size_t)t] || (skipSet != nullptr && (*skipSet)[(size_t)t]))
      {
        continue;
      }
      // The ray misses the triangle's box entirely if it never reaches it or
      // starts beyond it along every axis it moves in; a cheap reject.
      const double *box = &geometry.box[(size_t)6*t];
      bool miss = false;
      for (int k = 0; k < 3 && !miss; k++)
      {
        if (direction[k] > 0.0)
        {
          miss = origin[k] > box[2*k+1];
        }
        else if (direction[k] < 0.0)
        {
          miss = origin[k] < box[2*k];
        }
        else
        {
          miss = origin[k] < box[2*k] || origin[k] > box[2*k+1];
        }
      }
      if (miss)
      {
        continue;
      }
      const double *a = &points[(size_t)3*tris[(size_t)3*t]];
      const double *b = &points[(size_t)3*tris[(size_t)3*t + 1]];
      const double *cc = &points[(size_t)3*tris[(size_t)3*t + 2]];
      RayHit hit = RayTriangle(origin, direction, a, b, cc, edgeTol);
      if (hit.status == 2)
      {
        return false;
      }
      if (hit.status == 1)
      {
        winding += hit.sign;
      }
    }
    int axis = 0;
    if (tMax[1] < tMax[axis]) axis = 1;
    if (tMax[2] < tMax[axis]) axis = 2;
    if (step[axis] == 0)
    {
      break;
    }
    cell[axis] += step[axis];
    if (cell[axis] < 0 || cell[axis] >= grid.n[axis])
    {
      break;
    }
    tMax[axis] += tDelta[axis];
  }
  return true;
}

//---------------------
// The planar graph inside one triangle
//---------------------

struct LocalGraph
{
  std::vector<ll> id;
  std::vector<double> x, y;
  std::vector<std::vector<int> > adj;
  std::map<ll, int> local;
  // Which edges of the graph are crease segments (1) rather than pieces of
  // the triangle's own boundary (0), by sorted vertex pair.
  std::map<std::pair<int, int>, unsigned char> tag;
  double origin[3], axisU[3], axisV[3];

  void Frame(const double p0[3], const double p1[3], const double normal[3])
  {
    for (int k = 0; k < 3; k++)
    {
      origin[k] = p0[k];
      axisU[k] = p1[k] - p0[k];
    }
    Normalize(axisU);
    Cross(normal, axisU, axisV);
    Normalize(axisV);
  }

  int Vertex(ll gid, const double p[3])
  {
    std::map<ll, int>::iterator found = local.find(gid);
    if (found != local.end())
    {
      return found->second;
    }
    double r[3];
    Sub(p, origin, r);
    int v = (int)id.size();
    id.push_back(gid);
    x.push_back(Dot(r, axisU));
    y.push_back(Dot(r, axisV));
    adj.push_back(std::vector<int>());
    local[gid] = v;
    return v;
  }

  void Edge(int a, int b, unsigned char isCrease = 0)
  {
    if (a == b)
    {
      return;
    }
    std::pair<int, int> key(std::min(a, b), std::max(a, b));
    if (isCrease)
    {
      tag[key] = 1;
    }
    else if (tag.find(key) == tag.end())
    {
      tag[key] = 0;
    }
    for (size_t m = 0; m < adj[(size_t)a].size(); m++)
    {
      if (adj[(size_t)a][m] == b)
      {
        return;
      }
    }
    adj[(size_t)a].push_back(b);
    adj[(size_t)b].push_back(a);
  }

  unsigned char IsCrease(int a, int b) const
  {
    std::map<std::pair<int, int>, unsigned char>::const_iterator it = tag.find(std::make_pair(std::min(a, b), std::max(a, b)));
    return (it == tag.end()) ? 0 : it->second;
  }

  void SortAround()
  {
    for (size_t v = 0; v < adj.size(); v++)
    {
      std::vector<int> &around = adj[v];
      double cx = x[v], cy = y[v];
      std::sort(around.begin(), around.end(), [&](int p, int q)
      {
        return std::atan2(y[(size_t)p] - cy, x[(size_t)p] - cx) <
            std::atan2(y[(size_t)q] - cy, x[(size_t)q] - cx);
      });
    }
  }

  // The faces with the graph's edges on their left, as vertex cycles.
  void Faces(std::vector<std::vector<int> > &faces)
  {
    faces.clear();
    std::vector<std::vector<char> > used(adj.size());
    for (size_t v = 0; v < adj.size(); v++)
    {
      used[v].assign(adj[v].size(), 0);
    }
    for (size_t v = 0; v < adj.size(); v++)
    {
      for (size_t m = 0; m < adj[v].size(); m++)
      {
        if (used[v][m])
        {
          continue;
        }
        std::vector<int> cycle;
        int from = (int)v;
        int to = adj[v][m];
        used[v][m] = 1;
        cycle.push_back(from);
        size_t guard = 0;
        while (true)
        {
          // Arriving at 'to' from 'from': leave along the edge just
          // clockwise of the one back to 'from', which keeps the face on
          // the left.
          std::vector<int> &around = adj[(size_t)to];
          size_t pos = 0;
          for (; pos < around.size(); pos++)
          {
            if (around[pos] == from)
            {
              break;
            }
          }
          size_t nextPos = (pos + around.size() - 1) % around.size();
          int next = around[nextPos];
          if (to == (int)v && next == adj[v][m])
          {
            break;
          }
          if (used[(size_t)to][nextPos])
          {
            break;
          }
          used[(size_t)to][nextPos] = 1;
          cycle.push_back(to);
          from = to;
          to = next;
          if (++guard > 4*adj.size() + 8)
          {
            break;
          }
        }
        faces.push_back(cycle);
      }
    }
  }

  double SignedArea(const std::vector<int> &cycle) const
  {
    double area = 0.0;
    for (size_t m = 0; m < cycle.size(); m++)
    {
      size_t n = (m + 1) % cycle.size();
      area += x[(size_t)cycle[m]]*y[(size_t)cycle[n]] - x[(size_t)cycle[n]]*y[(size_t)cycle[m]];
    }
    return 0.5*area;
  }

  // Whether a point of the plane lies inside a polygon (crossing number).
  bool PointInPolygon(double px, double py, const std::vector<int> &poly) const
  {
    bool inside = false;
    size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++)
    {
      double xi = x[(size_t)poly[i]], yi = y[(size_t)poly[i]];
      double xj = x[(size_t)poly[j]], yj = y[(size_t)poly[j]];
      if ((yi > py) != (yj > py) && px < (xj - xi)*(py - yi)/(yj - yi) + xi)
      {
        inside = !inside;
      }
    }
    return inside;
  }

  // Whether two segments of the plane cross in their interiors; segments
  // that share an end never do.
  bool SegmentsCross(int a, int b, int c, int d) const
  {
    if (a == c || a == d || b == c || b == d)
    {
      return false;
    }
    auto orient = [&](int p, int q, int r)
    {
      return (x[(size_t)q] - x[(size_t)p])*(y[(size_t)r] - y[(size_t)p]) -
          (y[(size_t)q] - y[(size_t)p])*(x[(size_t)r] - x[(size_t)p]);
    };
    double o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b);
    return ((o1 > 0.0 && o2 < 0.0) || (o1 < 0.0 && o2 > 0.0)) &&
        ((o3 > 0.0 && o4 < 0.0) || (o3 < 0.0 && o4 > 0.0));
  }

  // Joins a hole, a clockwise cycle lying inside a counter-clockwise polygon,
  // into that polygon by a bridge edge run both ways, so that the result is
  // one weakly simple polygon the ear clipping can take (Eberly's method: the
  // hole's rightmost vertex is joined to the polygon vertex it can see along
  // and above the ray to its right). The bridge may cross no edge of the
  // polygon, of the hole, or of the other holes still to be joined. Returns
  // false if no bridge could be found; the polygon is then left as it was.
  bool BridgeHole(std::vector<int> &outer, const std::vector<int> &hole,
      const std::vector<std::vector<int> > &otherHoles) const
  {
    size_t hs = hole.size();
    if (hs < 3 || outer.size() < 3)
    {
      return false;
    }
    size_t mi = 0;
    for (size_t i = 1; i < hs; i++)
    {
      if (x[(size_t)hole[i]] > x[(size_t)hole[mi]] ||
          (x[(size_t)hole[i]] == x[(size_t)hole[mi]] && y[(size_t)hole[i]] > y[(size_t)hole[mi]]))
      {
        mi = i;
      }
    }
    int M = hole[mi];
    double mx = x[(size_t)M], my = y[(size_t)M];
    size_t n = outer.size();
    // The nearest point of the polygon the ray to the right of M meets.
    double bestX = std::numeric_limits<double>::max();
    int hitVertex = -1;
    size_t hitEdge = n;
    for (size_t i = 0; i < n; i++)
    {
      int a = outer[i], b = outer[(i + 1) % n];
      double ax = x[(size_t)a], ay = y[(size_t)a], bx = x[(size_t)b], by = y[(size_t)b];
      if (ay == my && ax >= mx && ax < bestX)
      {
        bestX = ax;
        hitVertex = a;
        hitEdge = n;
      }
      if ((ay > my) != (by > my))
      {
        double t = (my - ay)/(by - ay);
        double ix = ax + t*(bx - ax);
        if (ix >= mx && ix < bestX)
        {
          bestX = ix;
          hitEdge = i;
          hitVertex = -1;
        }
      }
    }
    if (hitVertex < 0 && hitEdge == n)
    {
      return false;
    }
    std::vector<int> candidates;
    if (hitVertex >= 0)
    {
      candidates.push_back(hitVertex);
    }
    else
    {
      int a = outer[hitEdge], b = outer[(hitEdge + 1) % n];
      int P = (x[(size_t)a] > x[(size_t)b]) ? a : b;
      double ix = bestX, iy = my;
      // Reflex vertices of the polygon inside the triangle M, I, P block the
      // view of P; the one nearest the ray in angle is seen instead.
      std::vector<std::pair<std::pair<double, double>, int> > blockers;
      for (size_t i = 0; i < n; i++)
      {
        int r = outer[i];
        if (r == P || r == M)
        {
          continue;
        }
        int prev = outer[(i + n - 1) % n], next = outer[(i + 1) % n];
        double cross = (x[(size_t)r] - x[(size_t)prev])*(y[(size_t)next] - y[(size_t)r]) -
            (y[(size_t)r] - y[(size_t)prev])*(x[(size_t)next] - x[(size_t)r]);
        if (!(cross < 0.0))
        {
          continue;
        }
        double rx = x[(size_t)r], ry = y[(size_t)r];
        double s1 = (ix - mx)*(ry - my) - (iy - my)*(rx - mx);
        double s2 = (x[(size_t)P] - ix)*(ry - iy) - (y[(size_t)P] - iy)*(rx - ix);
        double s3 = (mx - x[(size_t)P])*(ry - y[(size_t)P]) - (my - y[(size_t)P])*(rx - x[(size_t)P]);
        bool inside = (s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0) || (s1 <= 0.0 && s2 <= 0.0 && s3 <= 0.0);
        if (!inside)
        {
          continue;
        }
        double angle = std::atan2(std::abs(ry - my), rx - mx);
        double d2 = (rx - mx)*(rx - mx) + (ry - my)*(ry - my);
        blockers.push_back(std::make_pair(std::make_pair(angle, d2), r));
      }
      std::sort(blockers.begin(), blockers.end());
      for (size_t i = 0; i < blockers.size(); i++)
      {
        candidates.push_back(blockers[i].second);
      }
      candidates.push_back(P);
    }
    // Every other vertex to the right of M, nearest first, as a fallback.
    {
      std::vector<std::pair<double, int> > rest;
      for (size_t i = 0; i < n; i++)
      {
        int r = outer[i];
        if (r == M || !(x[(size_t)r] >= mx))
        {
          continue;
        }
        double d2 = (x[(size_t)r] - mx)*(x[(size_t)r] - mx) + (y[(size_t)r] - my)*(y[(size_t)r] - my);
        rest.push_back(std::make_pair(d2, r));
      }
      std::sort(rest.begin(), rest.end());
      for (size_t i = 0; i < rest.size(); i++)
      {
        candidates.push_back(rest[i].second);
      }
    }
    for (size_t c = 0; c < candidates.size(); c++)
    {
      int P = candidates[c];
      if (P == M)
      {
        continue;
      }
      bool blocked = false;
      for (size_t i = 0; i < n && !blocked; i++)
      {
        blocked = SegmentsCross(M, P, outer[i], outer[(i + 1) % n]);
      }
      for (size_t i = 0; i < hs && !blocked; i++)
      {
        blocked = SegmentsCross(M, P, hole[i], hole[(i + 1) % hs]);
      }
      for (size_t o = 0; o < otherHoles.size() && !blocked; o++)
      {
        const std::vector<int> &other = otherHoles[o];
        for (size_t i = 0; i < other.size() && !blocked; i++)
        {
          blocked = SegmentsCross(M, P, other[i], other[(i + 1) % other.size()]);
        }
      }
      if (blocked)
      {
        continue;
      }
      // The first occurrence of P; a vertex a previous bridge doubled is the
      // same point either way.
      size_t idx = 0;
      for (; idx < n; idx++)
      {
        if (outer[idx] == P)
        {
          break;
        }
      }
      std::vector<int> merged;
      merged.reserve(n + hs + 2);
      for (size_t i = 0; i <= idx; i++)
      {
        merged.push_back(outer[i]);
      }
      for (size_t k = 0; k <= hs; k++)
      {
        merged.push_back(hole[(mi + k) % hs]);
      }
      for (size_t i = idx; i < n; i++)
      {
        merged.push_back(outer[i]);
      }
      outer.swap(merged);
      return true;
    }
    return false;
  }

  // Drops the spikes a dangling edge leaves in a traced cycle (a vertex
  // whose neighbours in the cycle are the same vertex) and repeated
  // vertices, until the cycle is a plain polygon.
  void Simplify(std::vector<int> &cycle) const
  {
    bool changed = true;
    while (changed && cycle.size() >= 3)
    {
      changed = false;
      for (size_t m = 0; m < cycle.size(); m++)
      {
        size_t prev = (m + cycle.size() - 1) % cycle.size();
        size_t next = (m + 1) % cycle.size();
        if (cycle[prev] == cycle[next] || cycle[m] == cycle[next])
        {
          size_t drop = (cycle[m] == cycle[next]) ? next : m;
          cycle.erase(cycle.begin() + (long)drop);
          changed = true;
          break;
        }
      }
    }
  }

  // Ear clipping of a counter-clockwise polygon. Returns false if no ear
  // could be found at some step, in which case the remainder was fanned.
  bool Triangulate(const std::vector<int> &polygon, std::vector<int> &triangles) const
  {
    std::vector<int> ring(polygon);
    bool clean = true;
    while (ring.size() > 3)
    {
      size_t n = ring.size();
      bool clipped = false;
      for (size_t m = 0; m < n && !clipped; m++)
      {
        int a = ring[(m + n - 1) % n], b = ring[m], c = ring[(m + 1) % n];
        double ax = x[(size_t)a], ay = y[(size_t)a];
        double bx = x[(size_t)b], by = y[(size_t)b];
        double cx = x[(size_t)c], cy = y[(size_t)c];
        double cross = (bx - ax)*(cy - ay) - (by - ay)*(cx - ax);
        // An ear of no area - three points in a line, which a crease
        // running straight across a triangle leaves in every polygon it
        // bounds - is not an ear: rounding can make its turn look convex,
        // and clipping it would leave a sliver whose third edge is on
        // nothing.
        double span = (bx - ax)*(bx - ax) + (by - ay)*(by - ay) + (cx - ax)*(cx - ax) + (cy - ay)*(cy - ay);
        if (!(cross > 1.0e-12*span))
        {
          continue;
        }
        // A vertex on the ear's closing edge - the next point of a crease
        // running straight on - blocks it as surely as one inside it, and
        // rounding must not be allowed to say otherwise, so the test is
        // widened by a little.
        double blockTol = 1.0e-9*span;
        bool blocked = false;
        for (size_t o = 0; o < n && !blocked; o++)
        {
          int p = ring[o];
          if (p == a || p == b || p == c)
          {
            continue;
          }
          double px = x[(size_t)p], py = y[(size_t)p];
          double s1 = (bx - ax)*(py - ay) - (by - ay)*(px - ax);
          double s2 = (cx - bx)*(py - by) - (cy - by)*(px - bx);
          double s3 = (ax - cx)*(py - cy) - (ay - cy)*(px - cx);
          blocked = s1 >= -blockTol && s2 >= -blockTol && s3 >= -blockTol;
        }
        if (blocked)
        {
          continue;
        }
        triangles.push_back(a);
        triangles.push_back(b);
        triangles.push_back(c);
        ring.erase(ring.begin() + (long)m);
        clipped = true;
      }
      if (!clipped)
      {
        clean = false;
        for (size_t m = 1; m + 1 < ring.size(); m++)
        {
          triangles.push_back(ring[0]);
          triangles.push_back(ring[m]);
          triangles.push_back(ring[m + 1]);
        }
        return clean;
      }
    }
    if (ring.size() == 3)
    {
      triangles.push_back(ring[0]);
      triangles.push_back(ring[1]);
      triangles.push_back(ring[2]);
    }
    return clean;
  }
};

struct Segment
{
  ll p, q;      // point ids
  ll other;     // the triangle crossed along it
};

// The point three planes meet at, or false if they are nearly parallel.
bool ThreePlanes(const double n0[3], const double p0[3], const double n1[3], const double p1[3],
    const double n2[3], const double p2[3], double x[3])
{
  double c12[3], c20[3], c01[3];
  Cross(n1, n2, c12);
  Cross(n2, n0, c20);
  Cross(n0, n1, c01);
  double det = Dot(n0, c12);
  if (!(std::abs(det) > 1.0e-9))
  {
    return false;
  }
  double d0 = Dot(n0, p0), d1 = Dot(n1, p1), d2 = Dot(n2, p2);
  for (int k = 0; k < 3; k++)
  {
    x[k] = (d0*c12[k] + d1*c20[k] + d2*c01[k])/det;
  }
  return true;
}

void NoteFault(Report &report, const char *what, const double at[3])
{
  report.numArrangementFaults++;
  if (report.firstFault.empty())
  {
    report.firstFault = what;
    for (int k = 0; k < 3; k++)
    {
      report.firstFaultAt[k] = at[k];
    }
  }
}

} // namespace

//---------------------
// BuildOuterEnvelope
//---------------------

int BuildOuterEnvelope(const Surface &surface, Envelope &envelope, Report &report,
    std::string &error)
{
  auto startTime = std::chrono::steady_clock::now();
  report = Report();
  envelope = Envelope();

  ll numPts = (ll)(surface.points.size()/3);
  ll numTris = (ll)(surface.triangles.size()/3);
  ll numSheet = surface.numSheetTriangles;
  if (surface.points.size() != (size_t)3*numPts || surface.triangles.size() != (size_t)3*numTris)
  {
    error = "the surface arrays are not multiples of three";
    return 1;
  }
  if (numSheet < 0 || numSheet > numTris)
  {
    error = "the sheet triangle count is outside the triangle list";
    return 1;
  }
  for (size_t m = 0; m < surface.triangles.size(); m++)
  {
    if (surface.triangles[m] < 0 || surface.triangles[m] >= numPts)
    {
      error = "a triangle refers to a point outside the point list";
      return 1;
    }
  }
  if (numSheet == 0)
  {
    error = "there is no sheet to take the envelope of";
    return 1;
  }
  // The winding number needs a consistently wound closed surface: every edge
  // traversed once each way. Which way round the surface is wound - normals
  // out of the solid or into it - is not assumed but read off its signed
  // volume, and the rays below are shot to the outer side accordingly, so
  // the classification does not depend on the convention the caller used.
  {
    std::map<EdgeKey, std::pair<int, int> > use;
    for (ll t = 0; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = surface.triangles[(size_t)3*t + j], b = surface.triangles[(size_t)3*t + (j+1)%3];
        std::pair<int, int> &u = use[EdgeKey(a, b)];
        u.first++;
        if (a < b)
        {
          u.second++;
        }
      }
    }
    ll numOpen = 0, numMiswound = 0;
    for (std::map<EdgeKey, std::pair<int, int> >::iterator it = use.begin(); it != use.end(); ++it)
    {
      if (it->second.first != 2)
      {
        numOpen++;
      }
      else if (it->second.second != 1)
      {
        numMiswound++;
      }
    }
    if (numOpen > 0 || numMiswound > 0)
    {
      char what[160];
      snprintf(what, sizeof(what), "the closed surface is not a consistently wound closed surface: %lld edges are not on exactly two triangles and %lld are traversed the same way twice",
          numOpen, numMiswound);
      error = what;
      return 1;
    }
    double volume = 0.0;
    for (ll t = 0; t < numTris; t++)
    {
      const double *a = &surface.points[(size_t)3*surface.triangles[(size_t)3*t]];
      const double *b = &surface.points[(size_t)3*surface.triangles[(size_t)3*t + 1]];
      const double *c = &surface.points[(size_t)3*surface.triangles[(size_t)3*t + 2]];
      double bc[3];
      Cross(b, c, bc);
      volume += Dot(a, bc)/6.0;
    }
    if (!(std::abs(volume) > 0.0))
    {
      error = "the closed surface encloses no volume, so which side is out cannot be told";
      return 1;
    }
    report.windingSense = (volume > 0.0) ? 1 : -1;
    report.signedVolume = volume;
  }

  // The points grow with the crease points; the triangles do not.
  envelope.points = surface.points;
  envelope.pointKind.assign((size_t)numPts, 0);
  const std::vector<ll> &tris = surface.triangles;
  std::vector<double> &points = envelope.points;

  TriangleGeometry geometry;
  geometry.Build(points, tris, numTris);
  double meanExtent = 0.0;
  ll numDegenerate = 0;
  for (ll t = 0; t < numSheet; t++)
  {
    meanExtent += geometry.extent[(size_t)t];
    if (geometry.degenerate[(size_t)t])
    {
      numDegenerate++;
    }
  }
  meanExtent /= (double)numSheet;
  if (numDegenerate > 0)
  {
    const ll *tp = nullptr;
    for (ll t = 0; t < numSheet; t++)
    {
      if (geometry.degenerate[(size_t)t])
      {
        tp = &tris[(size_t)3*t];
        break;
      }
    }
    double at[3];
    for (int k = 0; k < 3; k++)
    {
      at[k] = (points[(size_t)3*tp[0] + k] + points[(size_t)3*tp[1] + k] + points[(size_t)3*tp[2] + k])/3.0;
    }
    char what[128];
    snprintf(what, sizeof(what), "%lld sheet triangles have no area", numDegenerate);
    NoteFault(report, what, at);
  }

  const ll maxBins = 4000000;
  Grid grid;
  grid.Build(geometry.box, numTris, meanExtent, maxBins);

  // Intersect every pair of sheet triangles whose boxes overlap and that do
  // not share an edge.
  const double baryTol = 1.0e-9;
  const double snapU = 1.0e-7;
  PointRegistry registry(envelope.points, envelope.pointKind, snapU, snapU*meanExtent);
  std::vector<std::vector<Segment> > segmentsOf((size_t)numSheet);
  std::vector<ll> stamp((size_t)numTris, 0);
  ll stampValue = 0;
  for (ll a = 0; a < numSheet; a++)
  {
    if (geometry.degenerate[(size_t)a])
    {
      continue;
    }
    const ll *ta = &tris[(size_t)3*a];
    const double *na = &geometry.normal[(size_t)3*a];
    const double *boxA = &geometry.box[(size_t)6*a];
    stampValue++;
    int i0 = grid.Bin(boxA[0], 0), i1 = grid.Bin(boxA[1], 0);
    int j0 = grid.Bin(boxA[2], 1), j1 = grid.Bin(boxA[3], 1);
    int k0 = grid.Bin(boxA[4], 2), k1 = grid.Bin(boxA[5], 2);
    for (int k = k0; k <= k1; k++)
    {
      for (int j = j0; j <= j1; j++)
      {
        for (int i = i0; i <= i1; i++)
        {
          ll bin = grid.Index(i, j, k);
          for (ll c = grid.start[(size_t)bin]; c < grid.start[(size_t)bin + 1]; c++)
          {
            ll b = grid.cells[(size_t)c];
            if (b <= a || b >= numSheet || stamp[(size_t)b] == stampValue)
            {
              continue;
            }
            stamp[(size_t)b] = stampValue;
            if (geometry.degenerate[(size_t)b])
            {
              continue;
            }
            if (!BoxesOverlap(boxA, &geometry.box[(size_t)6*b]))
            {
              continue;
            }
            report.numPairsTested++;
            const ll *tb = &tris[(size_t)3*b];
            const double *nb = &geometry.normal[(size_t)3*b];
            CrossingEnd ends[7];
            int numShared = 0;
            int found = CrossingEnds(points, ta, na, tb, nb, baryTol, ends, numShared);
            if (numShared >= 2 || found == 0)
            {
              continue;
            }
            double minLength = 1.0e-9*std::min(geometry.extent[(size_t)a], geometry.extent[(size_t)b]);
            // Every end found is a real crossing of an edge with a triangle,
            // and the triangle on the other side of that edge finds it too,
            // so every one is registered; the registry makes two that land
            // on one another one point.
            ll endId[7];
            for (int e = 0; e < found; e++)
            {
              const CrossingEnd &end = ends[e];
              if (end.kind == 0)
              {
                endId[e] = end.v0;
              }
              else
              {
                ll crossed = (end.kind == 1) ? b : a;
                endId[e] = registry.EdgePoint(end.v0, end.v1, crossed, end.u, end.p);
              }
            }
            CrossingEnd pair[2];
            if (!PickCrossing(ends, found, numShared, minLength, pair))
            {
              if (numShared == 0 && found == 1)
              {
                NoteFault(report, "a crossing with only one end found (a tolerance miss at an edge)", ends[0].p);
              }
              continue;
            }
            ll ids[2];
            for (int e = 0; e < 2; e++)
            {
              for (int f = 0; f < found; f++)
              {
                if (ends[f].kind == pair[e].kind && ends[f].v0 == pair[e].v0 && ends[f].v1 == pair[e].v1 &&
                    ends[f].u == pair[e].u)
                {
                  ids[e] = endId[f];
                  break;
                }
              }
            }
            if (ids[0] == ids[1])
            {
              continue;
            }
            report.numPairsCrossing++;
            Segment sa = {ids[0], ids[1], b};
            Segment sb = {ids[0], ids[1], a};
            segmentsOf[(size_t)a].push_back(sa);
            segmentsOf[(size_t)b].push_back(sb);
          }
        }
      }
    }
  }

  for (std::map<EdgeKey, std::vector<std::pair<double, ll> > >::iterator it = registry.onEdge.begin();
      it != registry.onEdge.end(); ++it)
  {
    std::sort(it->second.begin(), it->second.end());
  }
  report.numCreasePoints = 0;
  for (size_t p = (size_t)numPts; p < envelope.pointKind.size(); p++)
  {
    if (envelope.pointKind[p] == 1)
    {
      report.numCreasePoints++;
    }
  }

  // Split every crossed triangle into pieces along its segments, and gather
  // every piece of the sheet: the uncut triangles as they are.
  std::vector<ll> pieceTris;
  std::vector<ll> pieceSource;
  std::vector<unsigned char> pieceWhole;
  std::vector<unsigned char> pieceCreaseEdge;   // three per piece: edge j lies on a crease
  // The pieces of one face of the arrangement - one polygon of a split
  // triangle, or an uncut triangle - lie on one side of every crease, so
  // they are classified together, by a ray from the largest of them; a
  // sliver of a face along a crease could not be classified on its own.
  std::vector<ll> pieceFace;
  std::vector<ll> faceRepresentative;
  std::vector<double> faceRepresentativeArea;
  std::vector<std::vector<int> > faces;
  std::vector<int> polygon, fan;
  for (ll a = 0; a < numSheet; a++)
  {
    if (geometry.degenerate[(size_t)a])
    {
      continue;
    }
    const ll *ta = &tris[(size_t)3*a];
    std::vector<Segment> &segments = segmentsOf[(size_t)a];
    bool anyOnEdge = false;
    for (int j = 0; j < 3 && !anyOnEdge; j++)
    {
      anyOnEdge = registry.onEdge.find(EdgeKey(ta[j], ta[(j+1)%3])) != registry.onEdge.end();
    }
    if (segments.empty() && !anyOnEdge)
    {
      for (int j = 0; j < 3; j++)
      {
        pieceTris.push_back(ta[j]);
        pieceCreaseEdge.push_back(0);
      }
      pieceSource.push_back(a);
      pieceWhole.push_back(1);
      pieceFace.push_back((ll)faceRepresentative.size());
      faceRepresentative.push_back((ll)pieceSource.size() - 1);
      faceRepresentativeArea.push_back(1.0);
      continue;
    }
    report.numTrianglesSplit++;
    const double *na = &geometry.normal[(size_t)3*a];
    LocalGraph graph;
    graph.Frame(&points[(size_t)3*ta[0]], &points[(size_t)3*ta[1]], na);
    int corner[3];
    for (int j = 0; j < 3; j++)
    {
      corner[j] = graph.Vertex(ta[j], &points[(size_t)3*ta[j]]);
    }
    // The boundary, split at the points on each edge in order along it.
    for (int j = 0; j < 3; j++)
    {
      ll p = ta[j], q = ta[(j+1)%3];
      std::vector<int> chain;
      chain.push_back(corner[j]);
      std::map<EdgeKey, std::vector<std::pair<double, ll> > >::iterator on =
          registry.onEdge.find(EdgeKey(p, q));
      if (on != registry.onEdge.end())
      {
        const std::vector<std::pair<double, ll> > &list = on->second;
        // The list runs from the lower id; the edge runs p -> q.
        if (p < q)
        {
          for (size_t m = 0; m < list.size(); m++)
          {
            chain.push_back(graph.Vertex(list[m].second, &points[(size_t)3*list[m].second]));
          }
        }
        else
        {
          for (size_t m = list.size(); m-- > 0;)
          {
            chain.push_back(graph.Vertex(list[m].second, &points[(size_t)3*list[m].second]));
          }
        }
      }
      chain.push_back(corner[(j+1)%3]);
      for (size_t m = 0; m + 1 < chain.size(); m++)
      {
        graph.Edge(chain[m], chain[m+1]);
      }
    }
    // The segments, split where two from different sheets cross inside the
    // triangle. Such a crossing is a point three triangles meet at, and it
    // is keyed on the three so that each of them splits at the same point.
    struct Split
    {
      double s;
      int vertex;
    };
    std::vector<std::vector<Split> > splits(segments.size());
    std::vector<int> segP(segments.size()), segQ(segments.size());
    for (size_t s = 0; s < segments.size(); s++)
    {
      segP[s] = graph.Vertex(segments[s].p, &points[(size_t)3*segments[s].p]);
      segQ[s] = graph.Vertex(segments[s].q, &points[(size_t)3*segments[s].q]);
    }
    const double crossTol = 1.0e-9;
    for (size_t s = 0; s < segments.size(); s++)
    {
      for (size_t r = s + 1; r < segments.size(); r++)
      {
        if (segments[s].other == segments[r].other)
        {
          continue;
        }
        if (segP[s] == segP[r] || segP[s] == segQ[r] || segQ[s] == segP[r] || segQ[s] == segQ[r])
        {
          continue;
        }
        double p1x = graph.x[(size_t)segP[s]], p1y = graph.y[(size_t)segP[s]];
        double d1x = graph.x[(size_t)segQ[s]] - p1x, d1y = graph.y[(size_t)segQ[s]] - p1y;
        double p2x = graph.x[(size_t)segP[r]], p2y = graph.y[(size_t)segP[r]];
        double d2x = graph.x[(size_t)segQ[r]] - p2x, d2y = graph.y[(size_t)segQ[r]] - p2y;
        double denom = d1x*d2y - d1y*d2x;
        double scale = std::sqrt((d1x*d1x + d1y*d1y)*(d2x*d2x + d2y*d2y));
        if (!(std::abs(denom) > 1.0e-12*scale))
        {
          continue;
        }
        double rx = p2x - p1x, ry = p2y - p1y;
        double u = (rx*d2y - ry*d2x)/denom;
        double v = (rx*d1y - ry*d1x)/denom;
        if (u < -crossTol || u > 1.0 + crossTol || v < -crossTol || v > 1.0 + crossTol)
        {
          continue;
        }
        bool uAtEnd = u < crossTol || u > 1.0 - crossTol;
        bool vAtEnd = v < crossTol || v > 1.0 - crossTol;
        if (uAtEnd && vAtEnd)
        {
          continue;
        }
        if (uAtEnd)
        {
          Split sp = {v, (u < 0.5) ? segP[s] : segQ[s]};
          splits[r].push_back(sp);
          continue;
        }
        if (vAtEnd)
        {
          Split sp = {u, (v < 0.5) ? segP[r] : segQ[r]};
          splits[s].push_back(sp);
          continue;
        }
        ll b = segments[s].other, c = segments[r].other;
        double x3[3];
        if (!ThreePlanes(na, &points[(size_t)3*ta[0]],
              &geometry.normal[(size_t)3*b], &points[(size_t)3*tris[(size_t)3*b]],
              &geometry.normal[(size_t)3*c], &points[(size_t)3*tris[(size_t)3*c]], x3))
        {
          // Nearly parallel planes: lift the crossing found in this plane.
          double lx = p1x + u*d1x, ly = p1y + u*d1y;
          for (int k = 0; k < 3; k++)
          {
            x3[k] = graph.origin[k] + lx*graph.axisU[k] + ly*graph.axisV[k];
          }
        }
        bool created = false;
        ll id = registry.TriplePoint(a, b, c, x3, created);
        if (created)
        {
          report.numTriplePoints++;
        }
        int lv = graph.Vertex(id, &points[(size_t)3*id]);
        Split su = {u, lv};
        Split sv = {v, lv};
        splits[s].push_back(su);
        splits[r].push_back(sv);
      }
    }
    for (size_t s = 0; s < segments.size(); s++)
    {
      std::vector<Split> &sp = splits[s];
      std::sort(sp.begin(), sp.end(), [](const Split &l, const Split &r) { return l.s < r.s; });
      int last = segP[s];
      for (size_t m = 0; m < sp.size(); m++)
      {
        graph.Edge(last, sp[m].vertex, 1);
        last = sp[m].vertex;
      }
      graph.Edge(last, segQ[s], 1);
    }
    graph.SortAround();
    graph.Faces(faces);
    double triangleArea = 0.0;
    {
      std::vector<int> outline;
      outline.push_back(corner[0]);
      outline.push_back(corner[1]);
      outline.push_back(corner[2]);
      triangleArea = graph.SignedArea(outline);
    }
    double coveredArea = 0.0;
    size_t firstPiece = pieceSource.size();
    bool faulted = false;
    bool holeFaulted = false;
    // A face of positive area is a piece. A face of negative area is the
    // outside of one connected component of the graph: for the component
    // that holds the triangle's boundary it is the triangle run the wrong
    // way round, and for any other component - a crease loop closed inside
    // the triangle, left where a vertex of the other sheet pokes through -
    // it is the boundary of a hole in the piece around it. Such a hole is
    // joined to its piece by a bridge so that the piece can be triangulated;
    // the loop's own inside is a piece of its own like any other.
    std::vector<std::vector<int> > positive, holes;
    std::vector<double> positiveArea;
    for (size_t f = 0; f < faces.size(); f++)
    {
      polygon = faces[f];
      graph.Simplify(polygon);
      if (polygon.size() < 3)
      {
        continue;
      }
      double area = graph.SignedArea(polygon);
      if (area > 1.0e-12*triangleArea)
      {
        positive.push_back(polygon);
        positiveArea.push_back(area);
        coveredArea += area;
      }
      else if (area < -1.0e-12*triangleArea)
      {
        bool onBoundary = false;
        for (size_t m = 0; m < polygon.size() && !onBoundary; m++)
        {
          onBoundary = polygon[m] == corner[0] || polygon[m] == corner[1] || polygon[m] == corner[2];
        }
        if (!onBoundary)
        {
          holes.push_back(polygon);
          coveredArea += area;
        }
      }
    }
    std::vector<std::vector<size_t> > holesOf(positive.size());
    for (size_t h = 0; h < holes.size(); h++)
    {
      // The smallest piece the hole lies in is the one around it. The
      // loop's own inside faces share its vertices and are not around it,
      // whatever a point test on the loop's boundary says.
      int best = -1;
      double hx = graph.x[(size_t)holes[h][0]], hy = graph.y[(size_t)holes[h][0]];
      for (size_t f = 0; f < positive.size(); f++)
      {
        bool shares = false;
        for (size_t m = 0; m < positive[f].size() && !shares; m++)
        {
          for (size_t k = 0; k < holes[h].size() && !shares; k++)
          {
            shares = positive[f][m] == holes[h][k];
          }
        }
        if (shares)
        {
          continue;
        }
        if (graph.PointInPolygon(hx, hy, positive[f]) && (best < 0 || positiveArea[f] < positiveArea[(size_t)best]))
        {
          best = (int)f;
        }
      }
      if (best < 0)
      {
        holeFaulted = true;
        continue;
      }
      holesOf[(size_t)best].push_back(h);
    }
    for (size_t f = 0; f < positive.size(); f++)
    {
      polygon = positive[f];
      for (size_t k = 0; k < holesOf[f].size(); k++)
      {
        std::vector<std::vector<int> > others;
        for (size_t k2 = k + 1; k2 < holesOf[f].size(); k2++)
        {
          others.push_back(holes[holesOf[f][k2]]);
        }
        if (!graph.BridgeHole(polygon, holes[holesOf[f][k]], others))
        {
          holeFaulted = true;
        }
      }
      fan.clear();
      if (!graph.Triangulate(polygon, fan))
      {
        faulted = true;
      }
      ll face = (ll)faceRepresentative.size();
      faceRepresentative.push_back(-1);
      faceRepresentativeArea.push_back(-1.0);
      for (size_t m = 0; m + 2 < fan.size(); m += 3)
      {
        pieceTris.push_back(graph.id[(size_t)fan[m]]);
        pieceTris.push_back(graph.id[(size_t)fan[m+1]]);
        pieceTris.push_back(graph.id[(size_t)fan[m+2]]);
        for (int j = 0; j < 3; j++)
        {
          pieceCreaseEdge.push_back(graph.IsCrease(fan[m + j], fan[m + (j+1)%3]));
        }
        pieceSource.push_back(a);
        pieceWhole.push_back(0);
        pieceFace.push_back(face);
        std::vector<int> one(fan.begin() + (long)m, fan.begin() + (long)m + 3);
        double pieceArea = graph.SignedArea(one);
        if (pieceArea > faceRepresentativeArea[(size_t)face])
        {
          faceRepresentativeArea[(size_t)face] = pieceArea;
          faceRepresentative[(size_t)face] = (ll)pieceSource.size() - 1;
        }
      }
    }
    report.numPieces += (ll)(pieceSource.size() - firstPiece);
    double centre[3];
    for (int k = 0; k < 3; k++)
    {
      centre[k] = (points[(size_t)3*ta[0] + k] + points[(size_t)3*ta[1] + k] + points[(size_t)3*ta[2] + k])/3.0;
    }
    if (faulted)
    {
      NoteFault(report, "a piece of a crossed triangle could not be ear-clipped and was fanned", centre);
    }
    if (holeFaulted)
    {
      NoteFault(report, "a crease loop closed inside a triangle could not be joined to the piece around it", centre);
    }
    if (std::abs(coveredArea - triangleArea) > 1.0e-6*triangleArea)
    {
      NoteFault(report, "the pieces of a crossed triangle do not add up to it (a crease segment dangling inside it)", centre);
    }
  }

  // The count of turned-over triangles, for the log: where the extrusion
  // runs past the centre of curvature of a concave region the sheet turns
  // over, and its triangles there face the way they came from. With the
  // direction each triangle faced before the extrusion, which the caller
  // can supply, that is read exactly; without it, a triangle facing against
  // two of its edge-neighbours is taken for turned over. It is only a
  // count: what the fold leaves on the kept surface is caught below, as
  // pockets and shreds, because a turned-over triangle at the edge of the
  // sheet or at the end of a fold can still carry envelope.
  {
    if (surface.sheetNormal.size() == (size_t)3*numSheet)
    {
      for (ll t = 0; t < numSheet; t++)
      {
        double ref[3] = {surface.sheetNormal[(size_t)3*t], surface.sheetNormal[(size_t)3*t + 1], surface.sheetNormal[(size_t)3*t + 2]};
        if (!Normalize(ref))
        {
          continue;
        }
        if (geometry.degenerate[(size_t)t] || Dot(&geometry.normal[(size_t)3*t], ref) < 0.0)
        {
          report.numInverted++;
        }
      }
    }
    else
    {
      std::map<EdgeKey, std::vector<ll> > owners;
      for (ll t = 0; t < numSheet; t++)
      {
        for (int j = 0; j < 3; j++)
        {
          owners[EdgeKey(tris[(size_t)3*t + j], tris[(size_t)3*t + (j+1)%3])].push_back(t);
        }
      }
      std::vector<int> against((size_t)numSheet, 0);
      for (std::map<EdgeKey, std::vector<ll> >::iterator it = owners.begin(); it != owners.end(); ++it)
      {
        if (it->second.size() != 2)
        {
          continue;
        }
        ll t0 = it->second[0], t1 = it->second[1];
        if (Dot(&geometry.normal[(size_t)3*t0], &geometry.normal[(size_t)3*t1]) < 0.0)
        {
          against[(size_t)t0]++;
          against[(size_t)t1]++;
        }
      }
      for (ll t = 0; t < numSheet; t++)
      {
        if (against[(size_t)t] >= 2)
        {
          report.numInverted++;
        }
      }
    }
  }

  // Keep each piece whose outer side has winding number zero: nothing of the
  // solid lies on top of it. The ray goes roughly along the piece's normal,
  // scrambled a little so that it does not run along mesh edges, and is shot
  // again from a new direction when it grazes an edge.
  ll numPieceTris = (ll)pieceSource.size();
  ll numFaces = (ll)faceRepresentative.size();
  std::vector<int> windingOf((size_t)numPieceTris, 0);
  std::vector<unsigned char> keep((size_t)numPieceTris, 0);
  std::vector<unsigned char> decided((size_t)numPieceTris, 0);
  std::vector<int> faceWinding((size_t)numFaces, 0);
  std::vector<unsigned char> faceDecided((size_t)numFaces, 0);
  const double rayEdgeTol = 1.0e-9;
  const int maxRayTries = 24;
  std::map<int, ll> histogram;
  for (ll f = 0; f < numFaces; f++)
  {
    ll p = faceRepresentative[(size_t)f];
    if (p < 0)
    {
      continue;
    }
    const ll src = pieceSource[(size_t)p];
    const ll *pp = &pieceTris[(size_t)3*p];
    double centre[3];
    for (int k = 0; k < 3; k++)
    {
      centre[k] = (points[(size_t)3*pp[0] + k] + points[(size_t)3*pp[1] + k] + points[(size_t)3*pp[2] + k])/3.0;
    }
    // The outer side of the piece is along its geometric normal when the
    // surface is wound outward, and against it when wound inward.
    const double *n = &geometry.normal[(size_t)3*src];
    const double sense = (double)report.windingSense;
    int winding = 0;
    bool ok = false;
    for (int attempt = 0; attempt < maxRayTries && !ok; attempt++)
    {
      Scramble scramble((unsigned long long)p*1000003ULL + (unsigned long long)attempt*7919ULL + 17ULL);
      double spread = 0.05 + 0.1*attempt;
      double d[3];
      for (int k = 0; k < 3; k++)
      {
        d[k] = sense*n[k] + spread*scramble.Next();
      }
      if (!Normalize(d) || !(sense*Dot(d, n) > 0.2))
      {
        continue;
      }
      report.numRays++;
      if (attempt > 0)
      {
        report.numRayRetries++;
      }
      ok = WindingAlongRay(centre, d, pieceSource[(size_t)p], points, tris, geometry,
          grid, stamp, stampValue, rayEdgeTol, winding);
    }
    if (!ok)
    {
      continue;
    }
    faceDecided[(size_t)f] = 1;
    faceWinding[(size_t)f] = winding;
  }
  for (ll p = 0; p < numPieceTris; p++)
  {
    ll f = pieceFace[(size_t)p];
    if (!faceDecided[(size_t)f])
    {
      report.numUndecided++;
      continue;
    }
    decided[(size_t)p] = 1;
    windingOf[(size_t)p] = faceWinding[(size_t)f];
    histogram[windingOf[(size_t)p]]++;
    keep[(size_t)p] = (windingOf[(size_t)p] == 0) ? 1 : 0;
  }
  if (!histogram.empty())
  {
    int lo = histogram.begin()->first, hi = histogram.rbegin()->first;
    report.windingOffset = -lo;
    report.windingHistogram.assign((size_t)(hi - lo + 1), 0);
    for (std::map<int, ll>::iterator it = histogram.begin(); it != histogram.end(); ++it)
    {
      report.windingHistogram[(size_t)(it->first - lo)] = it->second;
    }
  }

  // Pockets and shreds. Where the extrusion folds over, at the cusps of a
  // swallowtail, the discrete sheet crumples: its outgoing and returning
  // layers overshoot the fold and cross, and between them lies a slit with
  // nothing of the solid in it. The pieces on the walls of a slit have
  // winding number zero on their outer side and pass as envelope, though
  // they lie deep inside the wall. What gives them away is that they are
  // cut off and enclosed: everything around them is inside the solid and
  // dropped, so they come out as a component of the kept surface on their
  // own, touching no cap rim - a closed pocket or an open shred - and the
  // rest of the closed surface winds around them. So each component without
  // a rim edge is put to the winding number once more, from its largest
  // piece, with the sheet triangles the component itself came from left
  // out of the count: what is left winds around a slit and around nothing
  // that is envelope. A closed surface of its own, however small or coarse,
  // has nothing else around it and stays; a wall with a hole in it stays
  // too, and the hole is reported.
  std::vector<unsigned char> pocket((size_t)numPieceTris, 0);
  {
    std::vector<ll> parent((size_t)numPieceTris);
    for (ll p = 0; p < numPieceTris; p++)
    {
      parent[(size_t)p] = p;
    }
    std::function<ll(ll)> find = [&](ll p) -> ll
    {
      while (parent[(size_t)p] != p)
      {
        parent[(size_t)p] = parent[(size_t)parent[(size_t)p]];
        p = parent[(size_t)p];
      }
      return p;
    };
    std::map<EdgeKey, std::vector<ll> > keptOn;
    for (ll p = 0; p < numPieceTris; p++)
    {
      if (!keep[(size_t)p])
      {
        continue;
      }
      for (int j = 0; j < 3; j++)
      {
        keptOn[EdgeKey(pieceTris[(size_t)3*p + j], pieceTris[(size_t)3*p + (j+1)%3])].push_back(p);
      }
    }
    // The rim edges are the edges the closing fans share with the sheet.
    std::map<EdgeKey, unsigned char> rimEdge;
    for (ll t = numSheet; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        EdgeKey key(tris[(size_t)3*t + j], tris[(size_t)3*t + (j+1)%3]);
        if (keptOn.count(key))
        {
          rimEdge[key] = 1;
        }
      }
    }
    for (std::map<EdgeKey, std::vector<ll> >::iterator it = keptOn.begin(); it != keptOn.end(); ++it)
    {
      const std::vector<ll> &on = it->second;
      for (size_t i = 1; i < on.size(); i++)
      {
        ll a = find(on[0]), b = find(on[i]);
        if (a != b)
        {
          parent[(size_t)b] = a;
        }
      }
    }
    struct Component { unsigned char rimmed; ll largestPiece; double largestArea; std::vector<ll> pieces; };
    std::map<ll, Component> components;
    for (ll p = 0; p < numPieceTris; p++)
    {
      if (!keep[(size_t)p])
      {
        continue;
      }
      const ll *pp = &pieceTris[(size_t)3*p];
      const double *a = &points[(size_t)3*pp[0]], *b = &points[(size_t)3*pp[1]], *c = &points[(size_t)3*pp[2]];
      double e1[3], e2[3], n[3];
      Sub(b, a, e1);
      Sub(c, a, e2);
      Cross(e1, e2, n);
      double area = 0.5*Norm(n);
      Component &comp = components[find(p)];
      if (comp.pieces.empty() || area > comp.largestArea)
      {
        comp.largestArea = area;
        comp.largestPiece = p;
      }
      comp.pieces.push_back(p);
    }
    for (std::map<EdgeKey, std::vector<ll> >::iterator it = keptOn.begin(); it != keptOn.end(); ++it)
    {
      if (it->second.size() == 1 && rimEdge.count(it->first))
      {
        components[find(it->second[0])].rimmed = 1;
      }
    }
    std::vector<unsigned char> ownSource((size_t)numTris, 0);
    for (std::map<ll, Component>::iterator it = components.begin(); it != components.end(); ++it)
    {
      const Component &comp = it->second;
      if (comp.rimmed)
      {
        continue;
      }
      for (size_t i = 0; i < comp.pieces.size(); i++)
      {
        ownSource[(size_t)pieceSource[(size_t)comp.pieces[i]]] = 1;
      }
      ll p = comp.largestPiece;
      const ll *pp = &pieceTris[(size_t)3*p];
      double centre[3];
      for (int k = 0; k < 3; k++)
      {
        centre[k] = (points[(size_t)3*pp[0] + k] + points[(size_t)3*pp[1] + k] + points[(size_t)3*pp[2] + k])/3.0;
      }
      const double *n = &geometry.normal[(size_t)3*pieceSource[(size_t)p]];
      const double sense = (double)report.windingSense;
      int winding = 0;
      bool ok = false;
      for (int attempt = 0; attempt < maxRayTries && !ok; attempt++)
      {
        Scramble scramble((unsigned long long)p*1000003ULL + (unsigned long long)attempt*7919ULL + 4243ULL);
        double spread = 0.05 + 0.1*attempt;
        double d[3];
        for (int k = 0; k < 3; k++)
        {
          d[k] = sense*n[k] + spread*scramble.Next();
        }
        if (!Normalize(d) || !(sense*Dot(d, n) > 0.2))
        {
          continue;
        }
        report.numRays++;
        ok = WindingAlongRay(centre, d, pieceSource[(size_t)p], points, tris, geometry,
            grid, stamp, stampValue, rayEdgeTol, winding, &ownSource);
      }
      for (size_t i = 0; i < comp.pieces.size(); i++)
      {
        ownSource[(size_t)pieceSource[(size_t)comp.pieces[i]]] = 0;
      }
      // Undecided is left as it is: nothing is dropped on a guess.
      if (!ok || winding == 0)
      {
        continue;
      }
      report.numPockets++;
      for (size_t i = 0; i < comp.pieces.size(); i++)
      {
        pocket[(size_t)comp.pieces[i]] = 1;
        keep[(size_t)comp.pieces[i]] = 0;
        report.numPocketPieces++;
      }
    }
  }

  // The result, and how it hangs together: every edge of it should be on
  // exactly two pieces, traversed once each way, except the rims the caller
  // closes.
  std::map<EdgeKey, std::pair<int, int> > edgeUse;   // (count, forward count)
  for (ll p = 0; p < numPieceTris; p++)
  {
    if (!keep[(size_t)p])
    {
      if (pieceWhole[(size_t)p] && decided[(size_t)p])
      {
        report.numWholeDropped++;
      }
      for (int j = 0; j < 3; j++)
      {
        envelope.droppedTriangles.push_back(pieceTris[(size_t)3*p + j]);
      }
      envelope.droppedSource.push_back(pieceSource[(size_t)p]);
      envelope.droppedWinding.push_back(pocket[(size_t)p] ? -998 : (decided[(size_t)p] ? windingOf[(size_t)p] : -999));
      continue;
    }
    const ll *pp = &pieceTris[(size_t)3*p];
    for (int j = 0; j < 3; j++)
    {
      envelope.triangles.push_back(pp[j]);
      envelope.creaseEdge.push_back(pieceCreaseEdge[(size_t)3*p + j]);
      ll a = pp[j], b = pp[(j+1)%3];
      std::pair<int, int> &use = edgeUse[EdgeKey(a, b)];
      use.first++;
      if (a < b)
      {
        use.second++;
      }
    }
    envelope.source.push_back(pieceSource[(size_t)p]);
    envelope.whole.push_back(pieceWhole[(size_t)p]);
    report.numPiecesKept++;
    if (pieceWhole[(size_t)p])
    {
      report.numWholeKept++;
    }
  }
  for (std::map<EdgeKey, std::pair<int, int> >::iterator it = edgeUse.begin(); it != edgeUse.end(); ++it)
  {
    int count = it->second.first, forward = it->second.second;
    if (count == 1)
    {
      report.numBoundaryEdges++;
    }
    else if (count > 2)
    {
      report.numNonManifoldEdges++;
      if (report.firstFault.empty())
      {
        report.firstFault = "an edge of the envelope on more than two pieces";
        for (int k = 0; k < 3; k++)
        {
          report.firstFaultAt[k] = 0.5*(points[(size_t)3*it->first.a + k] + points[(size_t)3*it->first.b + k]);
        }
      }
    }
    else if (forward != 1)
    {
      report.numMiswoundEdges++;
      if (report.firstFault.empty())
      {
        report.firstFault = "an edge of the envelope traversed the same way by both its pieces";
        for (int k = 0; k < 3; k++)
        {
          report.firstFaultAt[k] = 0.5*(points[(size_t)3*it->first.a + k] + points[(size_t)3*it->first.b + k]);
        }
      }
    }
  }
  report.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
  return 0;
}

//---------------------
// CleanEnvelopeSlivers
//---------------------

namespace {

// The aspect ratio of a triangle, longest edge over smallest altitude scaled
// so that an equilateral triangle is 1, the convention of the mesh quality
// reports. A triangle with no area is infinitely bad. Edge j runs from
// corner j to corner j+1; the shortest and longest are returned by index.
double TriangleAspect(const double *c[3], int &shortestEdge, int &longestEdge)
{
  double length[3];
  shortestEdge = 0;
  longestEdge = 0;
  for (int j = 0; j < 3; j++)
  {
    length[j] = Distance(c[j], c[(j+1)%3]);
    if (length[j] < length[shortestEdge]) shortestEdge = j;
    if (length[j] > length[longestEdge]) longestEdge = j;
  }
  double e1[3], e2[3], n[3];
  Sub(c[1], c[0], e1);
  Sub(c[2], c[0], e2);
  Cross(e1, e2, n);
  double twiceArea = Norm(n);
  if (!(twiceArea > 0.0) || !(length[longestEdge] > 0.0))
  {
    return std::numeric_limits<double>::infinity();
  }
  return length[longestEdge]*length[longestEdge]*std::sqrt(3.0)/(2.0*twiceArea);
}

// The distance from a point to a triangle (Ericson's closest-point cases).
double PointTriangleDistance(const double p[3], const double a[3], const double b[3], const double c[3])
{
  double ab[3], ac[3], ap[3];
  Sub(b, a, ab);
  Sub(c, a, ac);
  Sub(p, a, ap);
  double d1 = Dot(ab, ap), d2 = Dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return Distance(p, a);
  double bp[3];
  Sub(p, b, bp);
  double d3 = Dot(ab, bp), d4 = Dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return Distance(p, b);
  double vc = d1*d4 - d3*d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    double v = (d1 - d3 != 0.0) ? d1/(d1 - d3) : 0.0;
    double q[3] = {a[0] + v*ab[0], a[1] + v*ab[1], a[2] + v*ab[2]};
    return Distance(p, q);
  }
  double cp[3];
  Sub(p, c, cp);
  double d5 = Dot(ab, cp), d6 = Dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return Distance(p, c);
  double vb = d5*d2 - d1*d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    double w = (d2 - d6 != 0.0) ? d2/(d2 - d6) : 0.0;
    double q[3] = {a[0] + w*ac[0], a[1] + w*ac[1], a[2] + w*ac[2]};
    return Distance(p, q);
  }
  double va = d3*d6 - d5*d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    double w = ((d4 - d3) + (d5 - d6) != 0.0) ? (d4 - d3)/((d4 - d3) + (d5 - d6)) : 0.0;
    double q[3] = {b[0] + w*(c[0] - b[0]), b[1] + w*(c[1] - b[1]), b[2] + w*(c[2] - b[2])};
    return Distance(p, q);
  }
  double denom = va + vb + vc;
  if (!(denom != 0.0)) return Distance(p, a);
  double v = vb/denom, w = vc/denom;
  double q[3] = {a[0] + ab[0]*v + ac[0]*w, a[1] + ab[1]*v + ac[1]*w, a[2] + ab[2]*v + ac[2]*w};
  return Distance(p, q);
}

// The sliver cleanup's view of the surface: the triangles as they are being
// changed, the star of each point, and which points and edges are not to be
// touched. A star may hold triangles that no longer contain the point (a
// flip moves a triangle off a point) or that are gone, so every walk over
// a star checks both.
struct SliverMesh
{
  std::vector<double> &points;
  std::vector<ll> tris;
  std::vector<unsigned char> alive;
  std::vector<unsigned char> reshaped;
  std::vector<ll> source;
  std::vector<unsigned char> creaseFlag;   // three per triangle: edge j lies on a crease
  std::vector<std::vector<ll> > star;
  std::vector<unsigned char> fixed;    // per point: a rim point or one the caller fixed; never moves
  std::vector<unsigned char> crease;   // per point: a crease point; moves only onto another crease point, and then only within the deviation bound
  double aspectLimit;
  // The triangles binned where they started; a reshaped triangle stays
  // within an edge of where it was, so a query widened by a bin finds it.
  Grid grid;
  double meanExtent = 0.0;
  mutable std::vector<ll> stamp;
  mutable ll stampValue = 0;
  // Triangles made after the grid was built (by splitting), which the grid
  // does not hold and every crossing query has to look at as well.
  std::vector<ll> extraTris;
  // The kinds of the envelope's points, to mark a point that was moved onto
  // a crease as a crease point.
  std::vector<int> *kind = nullptr;

  SliverMesh(std::vector<double> &pts, double limit) : points(pts), aspectLimit(limit) {}

  void BuildGrid()
  {
    ll numTris = (ll)(tris.size()/3);
    TriangleGeometry geometry;
    geometry.Build(points, tris, numTris);
    meanExtent = 0.0;
    for (ll t = 0; t < numTris; t++)
    {
      meanExtent += geometry.extent[(size_t)t];
    }
    meanExtent /= std::max<ll>(numTris, 1);
    grid.Build(geometry.box, numTris, meanExtent, 4000000);
    stamp.assign((size_t)numTris, 0);
    stampValue = 0;
  }

  static void BoxOf(const double *c[3], double box[6])
  {
    for (int k = 0; k < 3; k++)
    {
      box[2*k] = std::min(c[0][k], std::min(c[1][k], c[2][k]));
      box[2*k + 1] = std::max(c[0][k], std::max(c[1][k], c[2][k]));
    }
  }

  // Whether a triangle on these corners would pass through any triangle of
  // the surface as it will stand: the live ones as they are, except those
  // being changed, which are taken as they will be, and those being taken
  // away. A crossing of some length only, as the volume mesher counts them.
  bool WouldCross(const ll tri[3], const std::vector<ll> &changedIds,
      const std::vector<std::array<ll, 3> > &changedTris, const std::vector<ll> &removedIds) const
  {
    const double *c[3] = {At(tri[0]), At(tri[1]), At(tri[2])};
    double n[3];
    NormalOf(c, n);
    double e1[3], e2[3];
    Sub(c[1], c[0], e1);
    Sub(c[2], c[0], e2);
    double longest2 = std::max(Dot(e1, e1), Dot(e2, e2));
    double area2 = Norm(n);
    if (!(area2 > 1.0e-14*longest2)) return true;
    for (int k = 0; k < 3; k++) n[k] /= area2;
    double box[6];
    BoxOf(c, box);
    double extent = 0.0;
    for (int k = 0; k < 3; k++) extent = std::max(extent, box[2*k + 1] - box[2*k]);
    const double baryTol = 1.0e-9;
    auto test = [&](const ll other[3]) -> bool
    {
      if ((other[0] == tri[0] || other[0] == tri[1] || other[0] == tri[2]) &&
          (other[1] == tri[0] || other[1] == tri[1] || other[1] == tri[2]) &&
          (other[2] == tri[0] || other[2] == tri[1] || other[2] == tri[2])) return false;
      const double *oc[3] = {At(other[0]), At(other[1]), At(other[2])};
      double on[3];
      NormalOf(oc, on);
      double oe1[3], oe2[3];
      Sub(oc[1], oc[0], oe1);
      Sub(oc[2], oc[0], oe2);
      double olongest2 = std::max(Dot(oe1, oe1), Dot(oe2, oe2));
      double oarea2 = Norm(on);
      if (!(oarea2 > 1.0e-14*olongest2)) return false;
      for (int k = 0; k < 3; k++) on[k] /= oarea2;
      double obox[6];
      BoxOf(oc, obox);
      if (!BoxesOverlap(box, obox)) return false;
      double oextent = 0.0;
      for (int k = 0; k < 3; k++) oextent = std::max(oextent, obox[2*k + 1] - obox[2*k]);
      CrossingEnd ends[7];
      int numShared = 0;
      int found = CrossingEnds(points, tri, n, other, on, baryTol, ends, numShared);
      if (numShared >= 2 || found == 0) return false;
      CrossingEnd pair[2];
      return PickCrossing(ends, found, numShared, 1.0e-7*std::min(extent, oextent), pair);
    };
    // The other triangles being changed, as they will be.
    for (size_t i = 0; i < changedTris.size(); i++)
    {
      if (test(changedTris[i].data())) return true;
    }
    // Everything else within reach, widened by a bin for what has moved.
    stampValue++;
    int i0 = std::max(0, grid.Bin(box[0], 0) - 1), i1 = std::min(grid.n[0] - 1, grid.Bin(box[1], 0) + 1);
    int j0 = std::max(0, grid.Bin(box[2], 1) - 1), j1 = std::min(grid.n[1] - 1, grid.Bin(box[3], 1) + 1);
    int k0 = std::max(0, grid.Bin(box[4], 2) - 1), k1 = std::min(grid.n[2] - 1, grid.Bin(box[5], 2) + 1);
    for (int k = k0; k <= k1; k++)
    {
      for (int j = j0; j <= j1; j++)
      {
        for (int i = i0; i <= i1; i++)
        {
          ll bin = grid.Index(i, j, k);
          for (ll cc = grid.start[(size_t)bin]; cc < grid.start[(size_t)bin + 1]; cc++)
          {
            ll o = grid.cells[(size_t)cc];
            if (stamp[(size_t)o] == stampValue) continue;
            stamp[(size_t)o] = stampValue;
            if (!alive[(size_t)o]) continue;
            bool skip = false;
            for (size_t r = 0; r < removedIds.size() && !skip; r++) if (removedIds[r] == o) skip = true;
            for (size_t r = 0; r < changedIds.size() && !skip; r++) if (changedIds[r] == o) skip = true;
            if (skip) continue;
            if (test(&tris[(size_t)3*o])) return true;
          }
        }
      }
    }
    for (size_t e = 0; e < extraTris.size(); e++)
    {
      ll o = extraTris[e];
      if (!alive[(size_t)o]) continue;
      bool skip = false;
      for (size_t r = 0; r < removedIds.size() && !skip; r++) if (removedIds[r] == o) skip = true;
      for (size_t r = 0; r < changedIds.size() && !skip; r++) if (changedIds[r] == o) skip = true;
      if (skip) continue;
      if (!BoxesOverlapTri(&tris[(size_t)3*o], box)) continue;
      if (test(&tris[(size_t)3*o])) return true;
    }
    return false;
  }

  bool BoxesOverlapTri(const ll other[3], const double box[6]) const
  {
    const double *oc[3] = {At(other[0]), At(other[1]), At(other[2])};
    double obox[6];
    BoxOf(oc, obox);
    return BoxesOverlap(box, obox);
  }

  const double *At(ll p) const { return &points[(size_t)3*p]; }

  bool Contains(ll t, ll p) const
  {
    return tris[(size_t)3*t] == p || tris[(size_t)3*t + 1] == p || tris[(size_t)3*t + 2] == p;
  }

  int CornerOf(ll t, ll p) const
  {
    for (int j = 0; j < 3; j++)
    {
      if (tris[(size_t)3*t + j] == p) return j;
    }
    return -1;
  }

  // The triangles on edge (a, b) that are alive, at most maxCount of them.
  int OnEdge(ll a, ll b, ll out[3], int maxCount) const
  {
    int n = 0;
    const std::vector<ll> &s = star[(size_t)a];
    for (size_t i = 0; i < s.size(); i++)
    {
      ll t = s[i];
      if (!alive[(size_t)t] || !Contains(t, a) || !Contains(t, b)) continue;
      bool seen = false;
      for (int k = 0; k < n; k++) if (out[k] == t) seen = true;
      if (seen) continue;
      if (n < maxCount) out[n] = t;
      n++;
    }
    return n;
  }

  // The flag of edge (p, q) of triangle t, or 0 if t has no such edge.
  unsigned char FlagOf(ll t, ll p, ll q) const
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3];
      if ((a == p && b == q) || (a == q && b == p)) return creaseFlag[(size_t)3*t + j];
    }
    return 0;
  }

  void SetFlag(ll t, ll p, ll q, unsigned char f)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3];
      if ((a == p && b == q) || (a == q && b == p)) creaseFlag[(size_t)3*t + j] = f;
    }
  }

  // An edge is a crease when the envelope marked it as one - a segment along
  // which two sheets cross - on any of the pieces on it. The dihedral across
  // it is whatever the sheets make, so it is never flipped and never read as
  // a fold. Two crease points joined by a piece of a sheet's own edge, or by
  // a diagonal of the splitting, are an ordinary edge.
  bool IsCreaseEdge(ll a, ll b) const
  {
    ll on[3];
    int n = OnEdge(a, b, on, 3);
    for (int k = 0; k < n && k < 3; k++)
    {
      if (FlagOf(on[k], a, b)) return true;
    }
    return false;
  }

  // The mean length of the edges of a set of triangles, the scale a move
  // in their midst is judged against.
  double MeanEdge(const std::vector<ll> &set) const
  {
    double sum = 0.0;
    int n = 0;
    for (size_t i = 0; i < set.size(); i++)
    {
      for (int j = 0; j < 3; j++)
      {
        sum += Distance(At(tris[(size_t)3*set[i] + j]), At(tris[(size_t)3*set[i] + (j+1)%3]));
        n++;
      }
    }
    return (n > 0) ? sum/n : 0.0;
  }

  void CornersAfter(ll t, ll from, ll to, const double *c[3], ll ids[3]) const
  {
    for (int j = 0; j < 3; j++)
    {
      ids[j] = tris[(size_t)3*t + j];
      if (ids[j] == from) ids[j] = to;
      c[j] = At(ids[j]);
    }
  }

  static void NormalOf(const double *c[3], double n[3])
  {
    double e1[3], e2[3];
    Sub(c[1], c[0], e1);
    Sub(c[2], c[0], e2);
    Cross(e1, e2, n);
  }

  double Aspect(ll t, int &shortest, int &longest) const
  {
    const double *c[3];
    ll ids[3];
    CornersAfter(t, -1, -1, c, ids);
    return TriangleAspect(c, shortest, longest);
  }

  // Two triangles that share an edge and are consistently wound have
  // normals to the same side unless one has folded over the edge, which
  // shows as a negative dot product.
  static bool Folded(const double n1[3], const double n2[3])
  {
    return !(Dot(n1, n2) > 0.0);
  }

  // Collapses point u onto point v when every check passes. The surface may
  // move by a twentieth of the reshaped triangles' longest edge, or by the
  // allowance the caller gives - the altitude below which a triangle is a
  // sliver, measured on the sliver being removed - whichever is more.
  bool Collapse(ll u, ll v, double allowance = 0.0)
  {
    if (u == v || fixed[(size_t)u]) return false;
    if (crease[(size_t)u] && !crease[(size_t)v]) return false;
    ll onEdge[3];
    if (OnEdge(u, v, onEdge, 3) != 2) return false;
    // The link condition: the points next to both u and v have to be
    // exactly the two apexes of the triangles on the edge, or the collapse
    // pinches the surface into a non-manifold.
    std::vector<ll> linkU, linkV;
    auto link = [&](ll p, std::vector<ll> &out)
    {
      const std::vector<ll> &s = star[(size_t)p];
      for (size_t i = 0; i < s.size(); i++)
      {
        ll t = s[i];
        if (!alive[(size_t)t] || !Contains(t, p)) continue;
        for (int j = 0; j < 3; j++)
        {
          ll q = tris[(size_t)3*t + j];
          if (q != p) out.push_back(q);
        }
      }
      std::sort(out.begin(), out.end());
      out.erase(std::unique(out.begin(), out.end()), out.end());
    };
    link(u, linkU);
    link(v, linkV);
    std::vector<ll> common;
    std::set_intersection(linkU.begin(), linkU.end(), linkV.begin(), linkV.end(), std::back_inserter(common));
    std::vector<ll> apex;
    for (int k = 0; k < 2; k++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll q = tris[(size_t)3*onEdge[k] + j];
        if (q != u && q != v) apex.push_back(q);
      }
    }
    std::sort(apex.begin(), apex.end());
    apex.erase(std::unique(apex.begin(), apex.end()), apex.end());
    if (common != apex || apex.size() != 2) return false;

    // The triangles that change: the rest of u's star, with u read as v.
    std::vector<ll> changed;
    const std::vector<ll> &su = star[(size_t)u];
    for (size_t i = 0; i < su.size(); i++)
    {
      ll t = su[i];
      if (!alive[(size_t)t] || !Contains(t, u) || t == onEdge[0] || t == onEdge[1]) continue;
      changed.push_back(t);
    }
    std::sort(changed.begin(), changed.end());
    changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
    if (changed.empty()) return false;
    const double moved = Distance(At(u), At(v));

    // The worst triangle touched, going in: the two on the edge and the ones
    // reshaped. No triangle may come out worse than that, or than the limit,
    // so the local worst never gets worse and the global worst never grows.
    double worstBefore = aspectLimit;
    for (int k = 0; k < 2; k++)
    {
      int s0, l0;
      worstBefore = std::max(worstBefore, Aspect(onEdge[k], s0, l0));
    }
    for (size_t i = 0; i < changed.size(); i++)
    {
      int s0, l0;
      worstBefore = std::max(worstBefore, Aspect(changed[i], s0, l0));
    }

    // The normal a triangle will have, read through the substitution and
    // with the two triangles on the edge gone.
    auto normalAfter = [&](ll t, double n[3])
    {
      const double *c[3];
      ll ids[3];
      CornersAfter(t, u, v, c, ids);
      NormalOf(c, n);
    };
    auto normalBefore = [&](ll t, double n[3])
    {
      const double *c[3];
      ll ids[3];
      CornersAfter(t, -1, -1, c, ids);
      NormalOf(c, n);
    };
    // The triangle on the other side of edge (x, y) after the collapse,
    // reading u as v; -1 for none, -2 for more than one.
    auto neighbourAfter = [&](ll t, ll x, ll y) -> ll
    {
      ll found = -1;
      for (int pass = 0; pass < 2; pass++)
      {
        ll p = (pass == 0) ? x : u;
        if (pass == 1 && x != v) break;
        const std::vector<ll> &s = star[(size_t)p];
        for (size_t i = 0; i < s.size(); i++)
        {
          ll m = s[i];
          if (m == t || m == found || !alive[(size_t)m] || m == onEdge[0] || m == onEdge[1]) continue;
          const double *c[3];
          ll ids[3];
          CornersAfter(m, u, v, c, ids);
          bool hasX = false, hasY = false;
          for (int j = 0; j < 3; j++)
          {
            if (ids[j] == x) hasX = true;
            if (ids[j] == y) hasY = true;
          }
          if (!hasX || !hasY) continue;
          if (found >= 0) return -2;
          found = m;
        }
      }
      return found;
    };

    for (size_t i = 0; i < changed.size(); i++)
    {
      ll t = changed[i];
      const double *before[3];
      ll idsBefore[3];
      CornersAfter(t, -1, -1, before, idsBefore);
      int s0, l0;
      double aspectBefore = TriangleAspect(before, s0, l0);
      double nBefore[3];
      NormalOf(before, nBefore);

      const double *after[3];
      ll ids[3];
      CornersAfter(t, u, v, after, ids);
      int s1, l1;
      double aspectAfter = TriangleAspect(after, s1, l1);
      if (!(aspectAfter <= worstBefore)) return false;
      double nAfter[3];
      NormalOf(after, nAfter);
      double lenAfter = Norm(nAfter);
      if (!(lenAfter > 0.0)) return false;
      // A sound triangle large against the move may not turn by more than
      // sixty degrees; a sliver has no direction to speak of, and a triangle
      // no bigger than the move is below what the move can be judged by, so
      // both are judged by their neighbours instead.
      if (aspectBefore <= aspectLimit && before[l0] != nullptr && Distance(before[l0], before[(l0+1)%3]) > 4.0*moved)
      {
        double lenBefore = Norm(nBefore);
        if (!(lenBefore > 0.0) || Dot(nBefore, nAfter) < 0.5*lenBefore*lenAfter) return false;
      }
      for (int j = 0; j < 3; j++)
      {
        ll x = ids[j], y = ids[(j+1)%3];
        if (x == y) return false;
        ll m = neighbourAfter(t, x, y);
        if (m == -2) return false;
        if (m < 0) continue;   // a boundary edge stays one
        if (IsCreaseEdge(x, y)) continue;
        // Two triangles no bigger than the move have no dihedral to speak
        // of; whether they fold is settled by the crossing check.
        if (!(Distance(after[l1], after[(l1+1)%3]) > 4.0*moved)) continue;
        {
          const double *cm[3];
          ll idsm[3];
          CornersAfter(m, u, v, cm, idsm);
          int sm, lm;
          TriangleAspect(cm, sm, lm);
          if (!(Distance(cm[lm], cm[(lm+1)%3]) > 4.0*moved)) continue;
        }
        double nm[3];
        normalAfter(m, nm);
        // A fold the move makes is refused; one that was there already, in
        // the crumple of a fold the extrusion left, is not this move's doing.
        if (Folded(nAfter, nm))
        {
          double nmBefore[3];
          normalBefore(m, nmBefore);
          if (!Folded(nBefore, nmBefore)) return false;
          // An existing fold may not deepen: the cosine across the edge may
          // not fall.
          double cosBefore = Dot(nBefore, nmBefore)/std::max(1.0e-300, Norm(nBefore)*Norm(nmBefore));
          double cosAfter = Dot(nAfter, nm)/std::max(1.0e-300, lenAfter*Norm(nm));
          if (cosAfter < cosBefore - 1.0e-12) return false;
        }
      }
    }

    // The surface may not move: the point taken away has to lie on the
    // reshaped triangles to within the bound, and may not have travelled
    // further than the bound along them either, or the collapse would
    // chamfer a ridge, pull a fold flat, or slide a point along a curved
    // surface. The bound is a twentieth of the mean edge around the point,
    // or the allowance, whichever is more; the mean rather than the longest
    // so that one long edge cannot buy a long move.
    {
      std::vector<ll> around(changed);
      around.push_back(onEdge[0]);
      around.push_back(onEdge[1]);
      double bound = std::max(0.05*MeanEdge(around), allowance);
      // A piece of the sheet itself may not slide either: along a curved
      // surface a long slide changes the shape though the point taken away
      // lies close to it. A piece a crease cut is judged by that distance
      // alone: where two creases run within the allowance of each other the
      // sheets there are as good as one, and a crease point may travel along
      // them to where they join.
      if (!(allowance > 0.0) && !(moved <= bound)) return false;
      double nearest = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < changed.size(); i++)
      {
        const double *c[3];
        ll ids[3];
        CornersAfter(changed[i], u, v, c, ids);
        nearest = std::min(nearest, PointTriangleDistance(At(u), c[0], c[1], c[2]));
      }
      if (!(nearest <= bound)) return false;
    }

    // Nothing reshaped may pass through the surface.
    {
      std::vector<std::array<ll, 3> > changedTris;
      std::vector<ll> removed;
      removed.push_back(onEdge[0]);
      removed.push_back(onEdge[1]);
      for (size_t i = 0; i < changed.size(); i++)
      {
        const double *c[3];
        std::array<ll, 3> ids;
        CornersAfter(changed[i], u, v, c, ids.data());
        changedTris.push_back(ids);
      }
      for (size_t i = 0; i < changed.size(); i++)
      {
        std::vector<std::array<ll, 3> > others;
        for (size_t k = 0; k < changedTris.size(); k++) if (k != i) others.push_back(changedTris[k]);
        if (WouldCross(changedTris[i].data(), changed, others, removed)) return false;
      }
    }

    // Apply. An edge that ran to u now runs to v and keeps its flag; where
    // two edges merge into one, the edge is a crease if either was.
    alive[(size_t)onEdge[0]] = 0;
    alive[(size_t)onEdge[1]] = 0;
    for (size_t i = 0; i < changed.size(); i++)
    {
      ll t = changed[i];
      int j = CornerOf(t, u);
      tris[(size_t)3*t + j] = v;
      reshaped[(size_t)t] = 1;
      star[(size_t)v].push_back(t);
    }
    star[(size_t)u].clear();
    for (size_t i = 0; i < changed.size(); i++)
    {
      ll t = changed[i];
      for (int j = 0; j < 3; j++)
      {
        ll p = tris[(size_t)3*t + j], q = tris[(size_t)3*t + (j+1)%3];
        if (p != v && q != v) continue;
        ll on[3];
        int n = OnEdge(p, q, on, 3);
        unsigned char f = 0;
        for (int k = 0; k < n && k < 3; k++) f |= FlagOf(on[k], p, q);
        for (int k = 0; k < n && k < 3; k++) SetFlag(on[k], p, q, f);
      }
    }
    return true;
  }

  // Moves the apex of a sliver onto the foot of its altitude on the long
  // edge, and splits the triangle across that edge there, when every check
  // passes. This is the move for a piece a crease cut nearly along one of its
  // edges: the apex is within the allowance of the long edge, so putting it
  // on the edge moves the surface by no more than that, the sliver becomes
  // nothing, and the triangle across the edge is split at the foot. A point
  // put on a crease edge becomes a crease point; a crease point may only be
  // put on a crease edge, so a crease is shortened or joined, never bent off
  // its line by more than the allowance.
  bool SnapToEdge(ll t, int le, double allowance)
  {
    ll v = tris[(size_t)3*t + le], w = tris[(size_t)3*t + (le+1)%3], u = tris[(size_t)3*t + (le+2)%3];
    if (fixed[(size_t)u]) return false;
    ll onEdge[3];
    if (OnEdge(v, w, onEdge, 3) != 2) return false;
    ll tp = (onEdge[0] == t) ? onEdge[1] : onEdge[0];
    if (tp == t || Contains(tp, u)) return false;
    int jw = CornerOf(tp, w);
    if (jw < 0 || tris[(size_t)3*tp + (jw+1)%3] != v) return false;   // tp has to run w -> v
    ll x = tris[(size_t)3*tp + (jw+2)%3];
    if (x == u) return false;
    ll onUX[3];
    if (OnEdge(u, x, onUX, 3) != 0) return false;   // u and x already joined
    bool creaseEdge = IsCreaseEdge(v, w);

    const double *pv = At(v), *pw = At(w);
    double vw[3];
    Sub(pw, pv, vw);
    double L2 = Dot(vw, vw);
    if (!(L2 > 0.0)) return false;
    double oldU[3] = {At(u)[0], At(u)[1], At(u)[2]};
    double vu[3];
    Sub(oldU, pv, vu);
    double s = Dot(vu, vw)/L2;
    const double sMin = 0.01;
    if (!(s > sMin && s < 1.0 - sMin)) return false;   // near an end the collapse is the move
    double f[3];
    for (int k = 0; k < 3; k++) f[k] = pv[k] + s*vw[k];
    // A crease point stays on a crease, unless the edge it would go to is so
    // close - within a tenth of the longest crease edge at the point - that
    // the crease is bent by less than its own resolution there.
    if (crease[(size_t)u] && !creaseEdge)
    {
      double longestCrease = 0.0;
      const std::vector<ll> &su0 = star[(size_t)u];
      for (size_t i = 0; i < su0.size(); i++)
      {
        ll m = su0[i];
        if (!alive[(size_t)m] || !Contains(m, u)) continue;
        for (int j = 0; j < 3; j++)
        {
          ll p = tris[(size_t)3*m + j], q = tris[(size_t)3*m + (j+1)%3];
          if ((p != u && q != u) || !creaseFlag[(size_t)3*m + j]) continue;
          longestCrease = std::max(longestCrease, Distance(At(p), At(q)));
        }
      }
      if (!(Distance(oldU, f) <= 0.1*longestCrease)) return false;
    }
    // The triangles that change: u's star but t, with u at the foot, and the
    // two halves of tp.
    std::vector<ll> changed;
    const std::vector<ll> &su = star[(size_t)u];
    for (size_t i = 0; i < su.size(); i++)
    {
      ll m = su[i];
      if (!alive[(size_t)m] || !Contains(m, u) || m == t) continue;
      changed.push_back(m);
    }
    std::sort(changed.begin(), changed.end());
    changed.erase(std::unique(changed.begin(), changed.end()), changed.end());
    if (changed.empty()) return false;
    // The move is the altitude; it may not exceed a twentieth of the mean
    // edge around the apex, or the allowance.
    {
      std::vector<ll> around(changed);
      around.push_back(t);
      around.push_back(tp);
      if (!(Distance(oldU, f) <= std::max(0.05*MeanEdge(around), allowance))) return false;
    }
    const double moved = Distance(oldU, f);
    std::vector<double> aspectBefore(changed.size()), longestBefore(changed.size());
    std::vector<std::array<double, 3> > normalBefore(changed.size());
    double worstBefore = aspectLimit;
    {
      int s0, l0;
      worstBefore = std::max(worstBefore, Aspect(t, s0, l0));
    }
    for (size_t i = 0; i < changed.size(); i++)
    {
      const double *c[3];
      ll ids[3];
      CornersAfter(changed[i], -1, -1, c, ids);
      int s0, l0;
      aspectBefore[i] = TriangleAspect(c, s0, l0);
      longestBefore[i] = Distance(c[l0], c[(l0+1)%3]);
      worstBefore = std::max(worstBefore, aspectBefore[i]);
      NormalOf(c, normalBefore[i].data());
    }
    const double *ctp[3];
    ll idsTp[3];
    CornersAfter(tp, -1, -1, ctp, idsTp);
    int s0, l0;
    double aspectTp = TriangleAspect(ctp, s0, l0);
    worstBefore = std::max(worstBefore, aspectTp);
    double nTp[3];
    NormalOf(ctp, nTp);

    // From here the point stands at the foot; it is put back on any failure.
    double *pu = &points[(size_t)3*u];
    for (int k = 0; k < 3; k++) pu[k] = f[k];
    auto restore = [&]()
    {
      for (int k = 0; k < 3; k++) pu[k] = oldU[k];
      return false;
    };

    std::vector<std::array<ll, 3> > newTris;
    for (size_t i = 0; i < changed.size(); i++)
    {
      std::array<ll, 3> ids = {{tris[(size_t)3*changed[i]], tris[(size_t)3*changed[i] + 1], tris[(size_t)3*changed[i] + 2]}};
      newTris.push_back(ids);
    }
    std::array<ll, 3> half1 = {{w, u, x}}, half2 = {{u, v, x}};
    newTris.push_back(half1);
    newTris.push_back(half2);
    std::vector<std::array<double, 3> > newNormal(newTris.size());
    for (size_t i = 0; i < newTris.size(); i++)
    {
      const double *c[3] = {At(newTris[i][0]), At(newTris[i][1]), At(newTris[i][2])};
      int s1, l1;
      double aspect = TriangleAspect(c, s1, l1);
      NormalOf(c, newNormal[i].data());
      double len = Norm(newNormal[i].data());
      if (!(len > 0.0)) return restore();
      if (!(aspect <= worstBefore)) return restore();
      if (i < changed.size())
      {
        // A sound triangle large against the move may not turn by more than
        // sixty degrees.
        if (aspectBefore[i] <= aspectLimit && longestBefore[i] > 4.0*moved)
        {
          double lenBefore = Norm(normalBefore[i].data());
          if (!(lenBefore > 0.0) || Dot(normalBefore[i].data(), newNormal[i].data()) < 0.5*lenBefore*len) return restore();
        }
      }
      else
      {
        if (!(Dot(nTp, newNormal[i].data()) > 0.0)) return restore();
      }
    }
    // The normal each new triangle had before the move, for telling a fold
    // the move makes from one that was there.
    auto oldNormalOf = [&](size_t i, double n[3])
    {
      if (i < changed.size())
      {
        for (int k = 0; k < 3; k++) n[k] = normalBefore[i][(size_t)k];
      }
      else
      {
        for (int k = 0; k < 3; k++) n[k] = nTp[k];
      }
    };
    // No fold across any edge of the new triangles, and no edge on more than
    // two of them. The edges along the crease, when the long edge is one,
    // are creases too.
    auto creaseLike = [&](ll p, ll q)
    {
      if (IsCreaseEdge(p, q)) return true;
      if (!creaseEdge) return false;
      return (p == u && (q == v || q == w)) || (q == u && (p == v || p == w));
    };
    auto longestOf = [&](const ll ids3[3])
    {
      double L = 0.0;
      for (int j = 0; j < 3; j++) L = std::max(L, Distance(At(ids3[j]), At(ids3[(j+1)%3])));
      return L;
    };
    for (size_t i = 0; i < newTris.size(); i++)
    {
      bool bigI = longestOf(newTris[i].data()) > 4.0*moved;
      for (int j = 0; j < 3; j++)
      {
        ll p = newTris[i][(size_t)j], q = newTris[i][(size_t)(j+1)%3];
        int numNeighbours = 0;
        bool folded = false;
        for (size_t k = 0; k < newTris.size(); k++)
        {
          if (k == i) continue;
          const std::array<ll, 3> &o = newTris[k];
          bool hasP = o[0] == p || o[1] == p || o[2] == p;
          bool hasQ = o[0] == q || o[1] == q || o[2] == q;
          if (!hasP || !hasQ) continue;
          numNeighbours++;
          if (!bigI || !(longestOf(o.data()) > 4.0*moved)) continue;
          if (!creaseLike(p, q) && Folded(newNormal[i].data(), newNormal[k].data()))
          {
            double oi[3], ok[3];
            oldNormalOf(i, oi);
            oldNormalOf(k, ok);
            if (!Folded(oi, ok)) folded = true;
            double cosBefore = Dot(oi, ok)/std::max(1.0e-300, Norm(oi)*Norm(ok));
            double cosAfter = Dot(newNormal[i].data(), newNormal[k].data())/std::max(1.0e-300, Norm(newNormal[i].data())*Norm(newNormal[k].data()));
            if (cosAfter < cosBefore - 1.0e-12) folded = true;
          }
        }
        ll on[3];
        int n = OnEdge(p, q, on, 3);
        for (int k = 0; k < n && k < 3; k++)
        {
          ll o = on[k];
          if (o == t || o == tp) continue;
          bool isChanged = false;
          for (size_t c = 0; c < changed.size() && !isChanged; c++) isChanged = changed[c] == o;
          if (isChanged) continue;
          numNeighbours++;
          const double *co[3];
          ll idso[3];
          CornersAfter(o, -1, -1, co, idso);
          if (!bigI || !(longestOf(idso) > 4.0*moved)) continue;
          double no[3];
          NormalOf(co, no);
          if (!creaseLike(p, q) && Folded(newNormal[i].data(), no))
          {
            double oi[3];
            oldNormalOf(i, oi);
            if (!Folded(oi, no)) folded = true;
            double cosBefore = Dot(oi, no)/std::max(1.0e-300, Norm(oi)*Norm(no));
            double cosAfter = Dot(newNormal[i].data(), no)/std::max(1.0e-300, Norm(newNormal[i].data())*Norm(no));
            if (cosAfter < cosBefore - 1.0e-12) folded = true;
          }
        }
        if (n > 3 || numNeighbours > 1 || folded) return restore();
      }
    }
    // Nothing new may pass through the surface.
    {
      std::vector<ll> changedIds(changed);
      changedIds.push_back(tp);
      std::vector<ll> removed;
      removed.push_back(t);
      removed.push_back(tp);
      for (size_t i = 0; i < newTris.size(); i++)
      {
        std::vector<std::array<ll, 3> > others;
        for (size_t k = 0; k < newTris.size(); k++) if (k != i) others.push_back(newTris[k]);
        if (WouldCross(newTris[i].data(), changedIds, others, removed)) return restore();
      }
    }

    // Apply. The halves of tp keep its flags on the edges they keep, the two
    // parts of the split edge take that edge's flag, and the new edge to the
    // foot is no crease.
    unsigned char flagVW = FlagOf(tp, v, w) | FlagOf(t, v, w);
    unsigned char flagWX = FlagOf(tp, w, x), flagXV = FlagOf(tp, x, v);
    alive[(size_t)t] = 0;
    tris[(size_t)3*tp] = w;
    tris[(size_t)3*tp + 1] = u;
    tris[(size_t)3*tp + 2] = x;
    creaseFlag[(size_t)3*tp] = flagVW;
    creaseFlag[(size_t)3*tp + 1] = 0;
    creaseFlag[(size_t)3*tp + 2] = flagWX;
    reshaped[(size_t)tp] = 1;
    ll fresh = (ll)(tris.size()/3);
    tris.push_back(u);
    tris.push_back(v);
    tris.push_back(x);
    creaseFlag.push_back(flagVW);
    creaseFlag.push_back(flagXV);
    creaseFlag.push_back(0);
    alive.push_back(1);
    reshaped.push_back(1);
    source.push_back(source[(size_t)tp]);
    extraTris.push_back(fresh);
    star[(size_t)u].push_back(tp);
    star[(size_t)u].push_back(fresh);
    star[(size_t)v].push_back(fresh);
    star[(size_t)x].push_back(fresh);
    if (creaseEdge)
    {
      crease[(size_t)u] = 1;
      if (kind != nullptr) (*kind)[(size_t)u] = 1;
    }
    return true;
  }

  // Takes a point out of the surface and triangulates the ring of its
  // neighbours afresh, when every check passes. This is the move for a point
  // a collapse cannot take away: a collapse joins the point to one
  // neighbour and needs the two links to meet only at the two apexes, which
  // in the crumple of a fold, where creases run into each other, they often
  // do not; the ring is triangulated without joining anything to anything
  // that is already joined, so it has no such condition. A crease point is
  // taken out only when exactly two crease edges meet at it, and the ring is
  // then triangulated on both sides of the edge joining those two, so the
  // crease runs straight through where the point was; a point where a
  // crease ends or three meet stays. The ring is triangulated in the plane
  // of the star's mean normal, taking at each step the ear of the best
  // shape.
  bool RemoveVertex(ll u, double allowance)
  {
    if (fixed[(size_t)u]) return false;
    // The star, and the ring edge each triangle has opposite u, in the order
    // the triangle runs (u, a, b): the ring runs a -> b.
    std::vector<ll> starTris;
    std::map<ll, ll> nextOf, triOf;
    const std::vector<ll> &su = star[(size_t)u];
    for (size_t i = 0; i < su.size(); i++)
    {
      ll m = su[i];
      if (!alive[(size_t)m] || !Contains(m, u)) continue;
      bool seen = false;
      for (size_t k = 0; k < starTris.size() && !seen; k++) seen = starTris[k] == m;
      if (seen) continue;
      starTris.push_back(m);
      int j = CornerOf(m, u);
      ll a = tris[(size_t)3*m + (j+1)%3], b = tris[(size_t)3*m + (j+2)%3];
      if (nextOf.count(a)) return false;   // not a simple ring
      nextOf[a] = b;
      triOf[a] = m;
    }
    size_t n = starTris.size();
    if (n < 3) return false;
    std::vector<ll> ring;
    ll start = nextOf.begin()->first, cur = start;
    for (size_t k = 0; k < n; k++)
    {
      ring.push_back(cur);
      std::map<ll, ll>::iterator it = nextOf.find(cur);
      if (it == nextOf.end()) return false;
      cur = it->second;
    }
    if (cur != start || ring.size() != n) return false;   // open or pinched
    for (size_t i = 0; i < n; i++) for (size_t k = i + 1; k < n; k++) if (ring[i] == ring[k]) return false;

    // The crease edges at u.
    std::vector<size_t> creaseAt;
    for (size_t i = 0; i < n; i++)
    {
      ll m = triOf[ring[i]];
      // edge (u, ring[i]) is on triangle m and on the previous triangle
      if (FlagOf(m, u, ring[i]) || FlagOf(triOf[ring[(i + n - 1) % n]], u, ring[i])) creaseAt.push_back(i);
    }
    if (creaseAt.size() != 0 && creaseAt.size() != 2) return false;
    bool creased = creaseAt.size() == 2;
    if (creased)
    {
      ll c1 = ring[creaseAt[0]], c2 = ring[creaseAt[1]];
      ll on[3];
      if ((creaseAt[1] + 1) % n != creaseAt[0] && (creaseAt[0] + 1) % n != creaseAt[1] && OnEdge(c1, c2, on, 3) != 0) return false;
    }

    // The frame: the star's mean normal. A point of the sheet itself (no
    // allowance) is taken out only where the sheet is smooth: across a
    // ridge the ring's triangles would cut the corner off, which no
    // distance bound on the point taken away catches.
    double meanN[3] = {0.0, 0.0, 0.0};
    double worstBefore = aspectLimit;
    std::vector<std::array<double, 3> > starNormal(n);
    for (size_t i = 0; i < n; i++)
    {
      const double *c[3];
      ll ids[3];
      CornersAfter(starTris[i], -1, -1, c, ids);
      NormalOf(c, starNormal[i].data());
      for (int k = 0; k < 3; k++) meanN[k] += starNormal[i][(size_t)k];
      int s0, l0;
      worstBefore = std::max(worstBefore, TriangleAspect(c, s0, l0));
    }
    if (!Normalize(meanN)) return false;
    if (!(allowance > 0.0))
    {
      for (size_t i = 0; i < n; i++)
      {
        double len = Norm(starNormal[i].data());
        if (!(len > 0.0)) continue;
        if (Dot(starNormal[i].data(), meanN) < 0.866*len) return false;   // more than thirty degrees off the mean: a ridge
      }
    }
    double axisU[3], axisV[3];
    {
      double r0[3];
      Sub(At(ring[0]), At(u), r0);
      double d = Dot(r0, meanN);
      for (int k = 0; k < 3; k++) axisU[k] = r0[k] - d*meanN[k];
      if (!Normalize(axisU)) return false;
      Cross(meanN, axisU, axisV);
    }
    std::vector<double> px(n), py(n);
    for (size_t i = 0; i < n; i++)
    {
      double r[3];
      Sub(At(ring[i]), At(u), r);
      px[i] = Dot(r, axisU);
      py[i] = Dot(r, axisV);
    }
    // The ring has to be a simple polygon in the frame, run the right way.
    {
      double area = 0.0;
      for (size_t i = 0; i < n; i++)
      {
        size_t j = (i + 1) % n;
        area += px[i]*py[j] - px[j]*py[i];
      }
      if (!(area > 0.0)) return false;
    }

    // Ear clipping of one chain of ring indices (a closed polygon), best ear
    // first; false if the polygon cannot be clipped.
    auto clip = [&](const std::vector<size_t> &poly, std::vector<std::array<size_t, 3> > &out) -> bool
    {
      std::vector<size_t> ringIdx(poly);
      while (ringIdx.size() > 3)
      {
        size_t m = ringIdx.size();
        double best = std::numeric_limits<double>::infinity();
        size_t bestAt = m;
        for (size_t i = 0; i < m; i++)
        {
          size_t a = ringIdx[(i + m - 1) % m], b = ringIdx[i], c = ringIdx[(i + 1) % m];
          double cross = (px[b] - px[a])*(py[c] - py[a]) - (py[b] - py[a])*(px[c] - px[a]);
          if (!(cross > 0.0)) continue;
          bool blocked = false;
          for (size_t o = 0; o < m && !blocked; o++)
          {
            size_t p = ringIdx[o];
            if (p == a || p == b || p == c) continue;
            double s1 = (px[b] - px[a])*(py[p] - py[a]) - (py[b] - py[a])*(px[p] - px[a]);
            double s2 = (px[c] - px[b])*(py[p] - py[b]) - (py[c] - py[b])*(px[p] - px[b]);
            double s3 = (px[a] - px[c])*(py[p] - py[c]) - (py[a] - py[c])*(px[p] - px[c]);
            blocked = s1 >= 0.0 && s2 >= 0.0 && s3 >= 0.0;
          }
          if (blocked) continue;
          const double *cc[3] = {At(ring[a]), At(ring[b]), At(ring[c])};
          int s0, l0;
          double aspect = TriangleAspect(cc, s0, l0);
          if (aspect < best)
          {
            best = aspect;
            bestAt = i;
          }
        }
        if (bestAt == m) return false;
        std::array<size_t, 3> tri = {{ringIdx[(bestAt + m - 1) % m], ringIdx[bestAt], ringIdx[(bestAt + 1) % m]}};
        out.push_back(tri);
        ringIdx.erase(ringIdx.begin() + (long)bestAt);
      }
      if (ringIdx.size() == 3)
      {
        std::array<size_t, 3> tri = {{ringIdx[0], ringIdx[1], ringIdx[2]}};
        out.push_back(tri);
      }
      return true;
    };
    std::vector<std::array<size_t, 3> > fans;
    if (creased)
    {
      // Two polygons, each running along the ring and closing across the
      // crease chord; each needs three vertices, else the chord is a ring
      // edge already and that side is nothing.
      size_t i1 = creaseAt[0], i2 = creaseAt[1];
      std::vector<size_t> side1, side2;
      for (size_t i = i1; ; i = (i + 1) % n) { side1.push_back(i); if (i == i2) break; }
      for (size_t i = i2; ; i = (i + 1) % n) { side2.push_back(i); if (i == i1) break; }
      if (side1.size() >= 3 && !clip(side1, fans)) return false;
      if (side2.size() >= 3 && !clip(side2, fans)) return false;
    }
    else
    {
      std::vector<size_t> all;
      for (size_t i = 0; i < n; i++) all.push_back(i);
      if (!clip(all, fans)) return false;
    }
    if (fans.empty()) return false;

    // The new triangles, their diagonals not already edges, none worse than
    // the worst of the star, no fold made or deepened across any edge,
    // and the point taken away close to them.
    std::vector<std::array<ll, 3> > newTris;
    std::vector<std::array<double, 3> > newNormal;
    for (size_t f = 0; f < fans.size(); f++)
    {
      std::array<ll, 3> ids = {{ring[fans[f][0]], ring[fans[f][1]], ring[fans[f][2]]}};
      const double *c[3] = {At(ids[0]), At(ids[1]), At(ids[2])};
      int s0, l0;
      double aspect = TriangleAspect(c, s0, l0);
      if (!(aspect <= worstBefore)) return false;
      std::array<double, 3> nn;
      NormalOf(c, nn.data());
      if (!(Norm(nn.data()) > 0.0)) return false;
      if (!(Dot(nn.data(), meanN) > 0.0)) return false;
      newTris.push_back(ids);
      newNormal.push_back(nn);
    }
    auto isRingEdge = [&](ll p, ll q) -> bool
    {
      for (size_t i = 0; i < n; i++)
      {
        ll a = ring[i], b = ring[(i + 1) % n];
        if ((a == p && b == q) || (a == q && b == p)) return true;
      }
      return false;
    };
    ll c1 = creased ? ring[creaseAt[0]] : -1, c2 = creased ? ring[creaseAt[1]] : -1;
    for (size_t i = 0; i < newTris.size(); i++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll p = newTris[i][(size_t)j], q = newTris[i][(size_t)(j+1)%3];
        bool ringEdge = isRingEdge(p, q);
        if (!ringEdge)
        {
          ll on[3];
          if (OnEdge(p, q, on, 3) != 0) return false;   // already an edge elsewhere
        }
        bool creaseLike = ringEdge ? IsCreaseEdge(p, q) : (creased && ((p == c1 && q == c2) || (p == c2 && q == c1)));
        // The neighbour across: another new triangle on a diagonal, the
        // outer triangle on a ring edge.
        int numNeighbours = 0;
        for (size_t k = 0; k < newTris.size(); k++)
        {
          if (k == i) continue;
          const std::array<ll, 3> &o = newTris[k];
          bool hasP = o[0] == p || o[1] == p || o[2] == p, hasQ = o[0] == q || o[1] == q || o[2] == q;
          if (!hasP || !hasQ) continue;
          numNeighbours++;
          if (creaseLike) continue;
          // Both new: a fold between them is one the move makes.
          if (Folded(newNormal[i].data(), newNormal[k].data())) return false;
        }
        if (ringEdge)
        {
          ll on[3];
          int cnt = OnEdge(p, q, on, 3);
          for (int k = 0; k < cnt && k < 3; k++)
          {
            ll o = on[k];
            bool inStar = false;
            for (size_t s = 0; s < starTris.size() && !inStar; s++) inStar = starTris[s] == o;
            if (inStar) continue;
            numNeighbours++;
            if (creaseLike) continue;
            const double *co[3];
            ll idso[3];
            CornersAfter(o, -1, -1, co, idso);
            double no[3];
            NormalOf(co, no);
            if (Folded(newNormal[i].data(), no))
            {
              // A fold the move makes is refused; one that was there between
              // the old star triangle on this edge and o may stay but not
              // deepen.
              ll starTri = -1;
              for (size_t s = 0; s < starTris.size(); s++)
              {
                if (Contains(starTris[s], p) && Contains(starTris[s], q)) starTri = starTris[s];
              }
              if (starTri < 0) return false;
              const double *cst[3];
              ll idst[3];
              CornersAfter(starTri, -1, -1, cst, idst);
              double nst[3];
              NormalOf(cst, nst);
              if (!Folded(nst, no)) return false;
              double cosBefore = Dot(nst, no)/std::max(1.0e-300, Norm(nst)*Norm(no));
              double cosAfter = Dot(newNormal[i].data(), no)/std::max(1.0e-300, Norm(newNormal[i].data())*Norm(no));
              if (cosAfter < cosBefore - 1.0e-12) return false;
            }
          }
        }
        if (numNeighbours > 1) return false;
      }
    }
    // The point taken away has to lie on the new triangles within the bound.
    {
      double bound = std::max(0.05*MeanEdge(starTris), allowance);
      double nearest = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < newTris.size(); i++)
      {
        nearest = std::min(nearest, PointTriangleDistance(At(u), At(newTris[i][0]), At(newTris[i][1]), At(newTris[i][2])));
      }
      if (!(nearest <= bound)) return false;
    }
    // Nothing new may pass through the surface.
    {
      std::vector<ll> none;
      for (size_t i = 0; i < newTris.size(); i++)
      {
        std::vector<std::array<ll, 3> > others;
        for (size_t k = 0; k < newTris.size(); k++) if (k != i) others.push_back(newTris[k]);
        if (WouldCross(newTris[i].data(), none, others, starTris)) return false;
      }
    }

    // Apply: the ring edges keep the flags of the star triangles they were
    // on, the crease chord is a crease, other diagonals are not.
    for (size_t i = 0; i < newTris.size(); i++)
    {
      ll fresh = (ll)(tris.size()/3);
      unsigned char flags[3];
      ll src = starTris[0];
      for (int j = 0; j < 3; j++)
      {
        ll p = newTris[i][(size_t)j], q = newTris[i][(size_t)(j+1)%3];
        flags[j] = 0;
        if (isRingEdge(p, q))
        {
          for (size_t s = 0; s < starTris.size(); s++)
          {
            if (Contains(starTris[s], p) && Contains(starTris[s], q))
            {
              flags[j] = FlagOf(starTris[s], p, q);
              src = starTris[s];
            }
          }
        }
        else if (creased && ((p == c1 && q == c2) || (p == c2 && q == c1)))
        {
          flags[j] = 1;
        }
      }
      for (int j = 0; j < 3; j++)
      {
        tris.push_back(newTris[i][(size_t)j]);
        creaseFlag.push_back(flags[j]);
      }
      alive.push_back(1);
      reshaped.push_back(1);
      source.push_back(source[(size_t)src]);
      extraTris.push_back(fresh);
      for (int j = 0; j < 3; j++) star[(size_t)newTris[i][(size_t)j]].push_back(fresh);
    }
    for (size_t s = 0; s < starTris.size(); s++) alive[(size_t)starTris[s]] = 0;
    star[(size_t)u].clear();
    return true;
  }

  // Flips edge j of triangle t when every check passes.
  bool Flip(ll t, int j)
  {
    ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3], c = tris[(size_t)3*t + (j+2)%3];
    if (IsCreaseEdge(a, b)) return false;
    ll onEdge[3];
    if (OnEdge(a, b, onEdge, 3) != 2) return false;
    ll m = (onEdge[0] == t) ? onEdge[1] : onEdge[0];
    if (m == t) return false;
    int ja = CornerOf(m, a), jb = CornerOf(m, b);
    if (ja < 0 || jb < 0 || (jb + 1)%3 != ja) return false;   // m has to run b -> a
    ll d = tris[(size_t)3*m + (ja+1)%3];
    if (d == c) return false;
    ll onNew[3];
    if (OnEdge(c, d, onNew, 3) != 0) return false;   // c and d already joined

    const double *ct[3], *cm[3];
    ll idsT[3], idsM[3];
    CornersAfter(t, -1, -1, ct, idsT);
    CornersAfter(m, -1, -1, cm, idsM);
    int s, l;
    double worstOld = std::max(TriangleAspect(ct, s, l), TriangleAspect(cm, s, l));
    double nt[3], nm[3];
    NormalOf(ct, nt);
    NormalOf(cm, nm);
    double lt = Norm(nt), lm = Norm(nm);
    if (!(lt > 0.0) || !(lm > 0.0)) return false;
    // Nearly coplanar: the flip is a change of diagonal in a flat quad.
    if (Dot(nt, nm) < 0.866*lt*lm) return false;
    double nRef[3] = {nt[0]/lt + nm[0]/lm, nt[1]/lt + nm[1]/lm, nt[2]/lt + nm[2]/lm};

    // New triangles (a, d, c) and (b, c, d).
    const double *c1[3] = {At(a), At(d), At(c)};
    const double *c2[3] = {At(b), At(c), At(d)};
    double worstNew = std::max(TriangleAspect(c1, s, l), TriangleAspect(c2, s, l));
    if (!(worstNew < worstOld)) return false;
    double n1[3], n2[3];
    NormalOf(c1, n1);
    NormalOf(c2, n2);
    if (!(Dot(n1, nRef) > 0.5*Norm(n1)*Norm(nRef)) || !(Dot(n2, nRef) > 0.5*Norm(n2)*Norm(nRef))) return false;

    // The surface may not move: the middle of the edge taken away has to lie
    // on a new triangle to within a twentieth of the quad's longest edge.
    {
      double mid[3] = {0.5*(At(a)[0] + At(b)[0]), 0.5*(At(a)[1] + At(b)[1]), 0.5*(At(a)[2] + At(b)[2])};
      double nearest = std::min(PointTriangleDistance(mid, c1[0], c1[1], c1[2]), PointTriangleDistance(mid, c2[0], c2[1], c2[2]));
      double longest = 0.0;
      const double *quad[4] = {At(a), At(c), At(b), At(d)};
      for (int k = 0; k < 4; k++)
      {
        longest = std::max(longest, Distance(quad[k], quad[(k+1)%4]));
      }
      if (!(nearest <= 0.05*longest)) return false;
    }

    // The four edges around the quad keep their outside neighbours; none
    // may fold against the new triangle on its side.
    struct Side { ll x, y; const double *n; ll self; };
    Side sides[4] = {{a, d, n1, m}, {c, a, n1, t}, {b, c, n2, t}, {d, b, n2, m}};
    for (int k = 0; k < 4; k++)
    {
      ll on[3];
      int n = OnEdge(sides[k].x, sides[k].y, on, 3);
      if (n > 2) return false;
      for (int i = 0; i < n; i++)
      {
        ll o = on[i];
        if (o == t || o == m) continue;
        if (IsCreaseEdge(sides[k].x, sides[k].y)) continue;
        const double *co[3];
        ll idso[3];
        CornersAfter(o, -1, -1, co, idso);
        double no[3];
        NormalOf(co, no);
        if (Folded(sides[k].n, no)) return false;
      }
    }

    // Neither new triangle may pass through the surface.
    {
      std::array<ll, 3> n1ids = {{a, d, c}}, n2ids = {{b, c, d}};
      std::vector<ll> changedIds;
      changedIds.push_back(t);
      changedIds.push_back(m);
      std::vector<ll> removed;
      std::vector<std::array<ll, 3> > other1(1, n2ids), other2(1, n1ids);
      if (WouldCross(n1ids.data(), changedIds, other1, removed)) return false;
      if (WouldCross(n2ids.data(), changedIds, other2, removed)) return false;
    }

    unsigned char fAD = FlagOf(m, a, d), fDB = FlagOf(m, d, b), fBC = FlagOf(t, b, c), fCA = FlagOf(t, c, a);
    tris[(size_t)3*t] = a; tris[(size_t)3*t + 1] = d; tris[(size_t)3*t + 2] = c;
    tris[(size_t)3*m] = b; tris[(size_t)3*m + 1] = c; tris[(size_t)3*m + 2] = d;
    creaseFlag[(size_t)3*t] = fAD; creaseFlag[(size_t)3*t + 1] = 0; creaseFlag[(size_t)3*t + 2] = fCA;
    creaseFlag[(size_t)3*m] = fBC; creaseFlag[(size_t)3*m + 1] = 0; creaseFlag[(size_t)3*m + 2] = fDB;
    reshaped[(size_t)t] = 1;
    reshaped[(size_t)m] = 1;
    star[(size_t)d].push_back(t);
    star[(size_t)c].push_back(m);
    return true;
  }
};

}

int CleanEnvelopeSlivers(Envelope &envelope, const std::vector<unsigned char> &fixedPoint,
    double aspectLimit, CleanReport &report)
{
  report = CleanReport();
  auto startTime = std::chrono::steady_clock::now();
  ll numPts = (ll)(envelope.points.size()/3);
  ll numTris = (ll)(envelope.triangles.size()/3);
  if (envelope.triangles.size() != (size_t)3*numTris || envelope.source.size() != (size_t)numTris ||
      envelope.whole.size() != (size_t)numTris || envelope.pointKind.size() != (size_t)numPts ||
      (!envelope.creaseEdge.empty() && envelope.creaseEdge.size() != (size_t)3*numTris) ||
      (!fixedPoint.empty() && fixedPoint.size() != (size_t)numPts) || !(aspectLimit > 1.0))
  {
    return 1;
  }
  for (size_t i = 0; i < envelope.triangles.size(); i++)
  {
    if (envelope.triangles[i] < 0 || envelope.triangles[i] >= numPts) return 1;
  }

  SliverMesh mesh(envelope.points, aspectLimit);
  mesh.tris = envelope.triangles;
  mesh.alive.assign((size_t)numTris, 1);
  mesh.reshaped.assign((size_t)numTris, 0);
  mesh.source = envelope.source;
  mesh.creaseFlag = envelope.creaseEdge;
  if (mesh.creaseFlag.empty())
  {
    // An envelope without the marks: fall back to reading a crease off the
    // points, an edge between two crease points whose pieces come from
    // different sheet triangles.
    mesh.creaseFlag.assign((size_t)3*numTris, 0);
    std::map<EdgeKey, std::vector<ll> > on;
    for (ll t = 0; t < numTris; t++)
      for (int j = 0; j < 3; j++)
        on[EdgeKey(mesh.tris[(size_t)3*t + j], mesh.tris[(size_t)3*t + (j+1)%3])].push_back(t);
    for (ll t = 0; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = mesh.tris[(size_t)3*t + j], b = mesh.tris[(size_t)3*t + (j+1)%3];
        if (envelope.pointKind[(size_t)a] == 0 || envelope.pointKind[(size_t)b] == 0) continue;
        const std::vector<ll> &o = on[EdgeKey(a, b)];
        bool differ = false;
        for (size_t k = 0; k < o.size(); k++) if (mesh.source[(size_t)o[k]] != mesh.source[(size_t)t]) differ = true;
        if (differ) mesh.creaseFlag[(size_t)3*t + j] = 1;
      }
    }
  }
  mesh.star.resize((size_t)numPts);
  mesh.fixed.assign((size_t)numPts, 0);
  mesh.crease.assign((size_t)numPts, 0);
  mesh.kind = &envelope.pointKind;
  for (ll t = 0; t < numTris; t++)
  {
    for (int j = 0; j < 3; j++)
    {
      mesh.star[(size_t)mesh.tris[(size_t)3*t + j]].push_back(t);
    }
  }
  mesh.BuildGrid();
  for (ll p = 0; p < numPts; p++)
  {
    if (envelope.pointKind[(size_t)p] != 0)
    {
      mesh.crease[(size_t)p] = 1;
    }
    if (!fixedPoint.empty() && fixedPoint[(size_t)p])
    {
      mesh.fixed[(size_t)p] = 1;
    }
  }
  // Boundary points stay: the boundary loops are the cap rims the caller
  // closes against. An edge on more than two pieces is not a surface this
  // can work on.
  {
    std::map<EdgeKey, int> edgeCount;
    for (ll t = 0; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        edgeCount[EdgeKey(mesh.tris[(size_t)3*t + j], mesh.tris[(size_t)3*t + (j+1)%3])]++;
      }
    }
    for (std::map<EdgeKey, int>::iterator it = edgeCount.begin(); it != edgeCount.end(); ++it)
    {
      if (it->second > 2) return 1;
      if (it->second == 1)
      {
        mesh.fixed[(size_t)it->first.a] = 1;
        mesh.fixed[(size_t)it->first.b] = 1;
      }
    }
  }

  auto survey = [&](ll &count, double &worst)
  {
    count = 0;
    worst = 0.0;
    ll numNow = (ll)(mesh.tris.size()/3);
    for (ll t = 0; t < numNow; t++)
    {
      if (!mesh.alive[(size_t)t]) continue;
      int s, l;
      double a = mesh.Aspect(t, s, l);
      if (a > aspectLimit) count++;
      if (a > worst) worst = a;
    }
  };
  survey(report.numSliversBefore, report.worstBefore);

  // The allowance a sliver gives its own removal. A piece a crease cut thin
  // is an artefact of the cutting, and the crumple of a fold around a crease
  // an artefact of the extrusion: such a piece may be flattened by its own
  // altitude, which is what the move needs and no more. A sliver of the
  // sheet itself, with no crease point on it, is geometry, and keeps the
  // twentieth-of-an-edge bound so that no ridge is chamfered.
  const int maxPasses = 30;
  for (int pass = 0; pass < maxPasses; pass++)
  {
    std::vector<std::pair<double, ll> > slivers;
    ll numNow = (ll)(mesh.tris.size()/3);
    for (ll t = 0; t < numNow; t++)
    {
      if (!mesh.alive[(size_t)t]) continue;
      int s, l;
      double a = mesh.Aspect(t, s, l);
      if (a > aspectLimit) slivers.push_back(std::make_pair(-a, t));
    }
    if (slivers.empty()) break;
    std::sort(slivers.begin(), slivers.end());
    report.numPasses = pass + 1;
    bool changed = false;
    for (size_t i = 0; i < slivers.size(); i++)
    {
      ll t = slivers[i].second;
      if (!mesh.alive[(size_t)t]) continue;
      int shortest, longest;
      if (!(mesh.Aspect(t, shortest, longest) > aspectLimit)) continue;
      // A snap moves the apex by the altitude; a collapse of the short edge
      // moves a point by that edge, which for a needle is about the altitude
      // and for a cap is half the long edge - a cap is snapped, not
      // collapsed, so a collapse is allowed the short edge only when that
      // edge is short against the long one or against the altitude.
      double allowanceSnap = 0.0, allowanceCollapse = 0.0;
      {
        const ll *tt = &mesh.tris[(size_t)3*t];
        bool creased = mesh.crease[(size_t)tt[0]] || mesh.crease[(size_t)tt[1]] || mesh.crease[(size_t)tt[2]];
        if (creased)
        {
          const double *c[3] = {mesh.At(tt[0]), mesh.At(tt[1]), mesh.At(tt[2])};
          double e1[3], e2[3], nn[3];
          Sub(c[1], c[0], e1);
          Sub(c[2], c[0], e2);
          Cross(e1, e2, nn);
          double L = Distance(c[longest], c[(longest+1)%3]);
          double s = Distance(c[shortest], c[(shortest+1)%3]);
          double h = (L > 0.0) ? Norm(nn)/L : 0.0;   // the altitude on the long edge
          allowanceSnap = h*(1.0 + 1.0e-6);
          allowanceCollapse = ((s <= 3.0*h || s <= 0.05*L) ? std::max(s, h) : h)*(1.0 + 1.0e-6);
        }
      }
      ll u = mesh.tris[(size_t)3*t + shortest], v = mesh.tris[(size_t)3*t + (shortest+1)%3];
      // Onto the end that must stay, else onto the crease end, else either.
      bool done = false;
      bool uStays = mesh.fixed[(size_t)u] || (mesh.crease[(size_t)u] && !mesh.crease[(size_t)v]);
      bool vStays = mesh.fixed[(size_t)v] || (mesh.crease[(size_t)v] && !mesh.crease[(size_t)u]);
      if (uStays || !vStays)
      {
        done = mesh.Collapse(v, u, allowanceCollapse);
      }
      if (!done)
      {
        done = mesh.Collapse(u, v, allowanceCollapse);
      }
      if (done)
      {
        report.numCollapsed++;
      }
      else if (mesh.SnapToEdge(t, longest, allowanceSnap))
      {
        done = true;
        report.numSnapped++;
      }
      else if (mesh.Flip(t, longest))
      {
        done = true;
        report.numFlipped++;
      }
      // The other two edges, when the shortest could not go.
      for (int j = 0; j < 3 && !done; j++)
      {
        if (j == shortest) continue;
        ll x = mesh.tris[(size_t)3*t + j], y = mesh.tris[(size_t)3*t + (j+1)%3];
        if (mesh.Collapse(x, y, allowanceCollapse) || mesh.Collapse(y, x, allowanceCollapse))
        {
          done = true;
          report.numCollapsed++;
        }
      }
      // Taking a point out altogether, when no joining was allowed: the
      // apex of the long edge first, then the ends of the short edge.
      if (!done)
      {
        ll order[3] = {mesh.tris[(size_t)3*t + (longest+2)%3], mesh.tris[(size_t)3*t + shortest], mesh.tris[(size_t)3*t + (shortest+1)%3]};
        for (int k = 0; k < 3 && !done; k++)
        {
          if (mesh.RemoveVertex(order[k], std::max(allowanceSnap, allowanceCollapse)))
          {
            done = true;
            report.numRemoved++;
          }
        }
      }
      if (done) changed = true;
    }
    if (!changed) break;
  }
  survey(report.numSliversAfter, report.worstAfter);
  report.numNoMoveAllowed = report.numSliversAfter;
  // Where the worst piece left is, and how many of those left touch a
  // crease, for the log.
  {
    ll numNow = (ll)(mesh.tris.size()/3);
    double worst = -1.0;
    for (ll t = 0; t < numNow; t++)
    {
      if (!mesh.alive[(size_t)t]) continue;
      int s, l;
      double a = mesh.Aspect(t, s, l);
      const ll *tt = &mesh.tris[(size_t)3*t];
      bool creased = mesh.crease[(size_t)tt[0]] || mesh.crease[(size_t)tt[1]] || mesh.crease[(size_t)tt[2]];
      if (a > aspectLimit && creased) report.numLeftCreased++;
      if (a > worst)
      {
        worst = a;
        for (int k = 0; k < 3; k++)
        {
          report.worstAfterAt[k] = (mesh.At(tt[0])[k] + mesh.At(tt[1])[k] + mesh.At(tt[2])[k])/3.0;
        }
      }
    }
  }

  // Write back, compacted.
  std::vector<ll> triangles, source;
  std::vector<unsigned char> whole, creaseEdge;
  ll numFinal = (ll)(mesh.tris.size()/3);
  for (ll t = 0; t < numFinal; t++)
  {
    if (!mesh.alive[(size_t)t]) continue;
    for (int j = 0; j < 3; j++)
    {
      triangles.push_back(mesh.tris[(size_t)3*t + j]);
      creaseEdge.push_back(mesh.creaseFlag[(size_t)3*t + j]);
    }
    source.push_back(mesh.source[(size_t)t]);
    whole.push_back((t >= numTris || mesh.reshaped[(size_t)t]) ? (unsigned char)2 : envelope.whole[(size_t)t]);
  }
  envelope.triangles.swap(triangles);
  envelope.source.swap(source);
  envelope.whole.swap(whole);
  envelope.creaseEdge.swap(creaseEdge);

  report.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
  return 0;
}

//---------------------
// CountCrossingTriangles
//---------------------

long long CountCrossingTriangles(const std::vector<double> &points,
    const std::vector<long long> &triangles, std::vector<unsigned char> &crossing,
    double firstAt[3])
{
  ll numTris = (ll)(triangles.size()/3);
  crossing.assign((size_t)numTris, 0);
  firstAt[0] = firstAt[1] = firstAt[2] = 0.0;
  if (numTris == 0)
  {
    return 0;
  }
  TriangleGeometry geometry;
  geometry.Build(points, triangles, numTris);
  double meanExtent = 0.0;
  for (ll t = 0; t < numTris; t++)
  {
    meanExtent += geometry.extent[(size_t)t];
  }
  meanExtent /= (double)numTris;
  Grid grid;
  grid.Build(geometry.box, numTris, meanExtent, 4000000);
  std::vector<ll> stamp((size_t)numTris, 0);
  ll stampValue = 0;
  const double baryTol = 1.0e-9;
  ll numCrossing = 0;
  ll first = -1;
  for (ll a = 0; a < numTris; a++)
  {
    if (geometry.degenerate[(size_t)a])
    {
      continue;
    }
    const ll *ta = &triangles[(size_t)3*a];
    const double *na = &geometry.normal[(size_t)3*a];
    const double *boxA = &geometry.box[(size_t)6*a];
    stampValue++;
    int i0 = grid.Bin(boxA[0], 0), i1 = grid.Bin(boxA[1], 0);
    int j0 = grid.Bin(boxA[2], 1), j1 = grid.Bin(boxA[3], 1);
    int k0 = grid.Bin(boxA[4], 2), k1 = grid.Bin(boxA[5], 2);
    for (int k = k0; k <= k1; k++)
    {
      for (int j = j0; j <= j1; j++)
      {
        for (int i = i0; i <= i1; i++)
        {
          ll bin = grid.Index(i, j, k);
          for (ll c = grid.start[(size_t)bin]; c < grid.start[(size_t)bin + 1]; c++)
          {
            ll b = grid.cells[(size_t)c];
            if (b <= a || stamp[(size_t)b] == stampValue)
            {
              continue;
            }
            stamp[(size_t)b] = stampValue;
            if (geometry.degenerate[(size_t)b] || !BoxesOverlap(boxA, &geometry.box[(size_t)6*b]))
            {
              continue;
            }
            CrossingEnd ends[7];
            int numShared = 0;
            int found = CrossingEnds(points, ta, na, &triangles[(size_t)3*b],
                &geometry.normal[(size_t)3*b], baryTol, ends, numShared);
            if (numShared >= 2 || found == 0)
            {
              continue;
            }
            // A crossing that a volume mesher would refuse is one of some
            // length: two triangles meeting at a shared crease point, or an
            // edge of one touching the other, are not.
            double minLength = 1.0e-7*std::min(geometry.extent[(size_t)a], geometry.extent[(size_t)b]);
            CrossingEnd pair[2];
            if (!PickCrossing(ends, found, numShared, minLength, pair))
            {
              continue;
            }
            if (!crossing[(size_t)a])
            {
              crossing[(size_t)a] = 1;
              numCrossing++;
              if (first < 0)
              {
                first = a;
              }
            }
            if (!crossing[(size_t)b])
            {
              crossing[(size_t)b] = 1;
              numCrossing++;
            }
          }
        }
      }
    }
  }
  if (first >= 0)
  {
    const ll *tp = &triangles[(size_t)3*first];
    for (int k = 0; k < 3; k++)
    {
      firstAt[k] = (points[(size_t)3*tp[0] + k] + points[(size_t)3*tp[1] + k] + points[(size_t)3*tp[2] + k])/3.0;
    }
  }
  return numCrossing;
}

}
