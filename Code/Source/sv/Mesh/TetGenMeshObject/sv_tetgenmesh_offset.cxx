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
 * @file sv_tetgenmesh_offset.cxx
 * @brief The offset outer surface of the wall, contoured from a distance
 * field sampled on a point cloud that follows the input. See the header.
 */

#include "sv_tetgenmesh_offset.h"
#include "sv_tetgenmesh_envelope.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <queue>
#include <unordered_map>

namespace svoffset
{

namespace
{

typedef long long ll;

//---------------------
// Vector helpers
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
  if (!(n > 0.0) || !std::isfinite(n))
  {
    return false;
  }
  a[0] /= n;
  a[1] /= n;
  a[2] /= n;
  return true;
}

// The closest point of triangle abc to p, with its barycentric weights;
// returns the distance. Ericson's region walk.
double ClosestOnTriangle(const double p[3], const double a[3], const double b[3],
    const double c[3], double q[3], double bary[3])
{
  double ab[3], ac[3], ap[3];
  Sub(b, a, ab);
  Sub(c, a, ac);
  Sub(p, a, ap);
  auto at = [&](double v, double w)
  {
    bary[0] = 1.0 - v - w;
    bary[1] = v;
    bary[2] = w;
    for (int k = 0; k < 3; k++)
    {
      q[k] = a[k] + v*ab[k] + w*ac[k];
    }
    return Distance(p, q);
  };
  double d1 = Dot(ab, ap), d2 = Dot(ac, ap);
  if (d1 <= 0.0 && d2 <= 0.0) return at(0.0, 0.0);
  double bp[3];
  Sub(p, b, bp);
  double d3 = Dot(ab, bp), d4 = Dot(ac, bp);
  if (d3 >= 0.0 && d4 <= d3) return at(1.0, 0.0);
  double vc = d1*d4 - d3*d2;
  if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
  {
    double v = (d1 - d3 != 0.0) ? d1/(d1 - d3) : 0.0;
    return at(v, 0.0);
  }
  double cp[3];
  Sub(p, c, cp);
  double d5 = Dot(ab, cp), d6 = Dot(ac, cp);
  if (d6 >= 0.0 && d5 <= d6) return at(0.0, 1.0);
  double vb = d5*d2 - d1*d6;
  if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
  {
    double w = (d2 - d6 != 0.0) ? d2/(d2 - d6) : 0.0;
    return at(0.0, w);
  }
  double va = d3*d6 - d5*d4;
  if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
  {
    double w = (d4 - d3)/((d4 - d3) + (d5 - d6));
    return at(1.0 - w, w);
  }
  double den = va + vb + vc;
  if (den == 0.0) return at(0.0, 0.0);
  return at(vb/den, vc/den);
}

//---------------------
// Triangle grid
//---------------------
// A uniform grid of bins over the field's surface, each bin listing the
// triangles whose boxes touch it and the largest thickness among their
// corners, so that a search for the least (distance - thickness) can skip a
// bin whose box is farther than the best found so far plus that thickness.
struct TriangleGrid
{
  const std::vector<double> *points = nullptr;
  const std::vector<ll> *triangles = nullptr;
  double lo[3] = {0.0, 0.0, 0.0};
  double cell = 1.0;
  int n[3] = {1, 1, 1};
  std::vector<ll> start, cells;
  std::vector<double> binMaxThickness;

  int Bin(double v, int k) const
  {
    double f = std::floor((v - lo[k])/cell);
    if (f < 0.0) return 0;
    if (f >= (double)n[k]) return n[k] - 1;
    return (int)f;
  }

  size_t Index(int i, int j, int k) const
  {
    return ((size_t)k*n[1] + j)*n[0] + i;
  }

  void Build(const std::vector<double> &pts, const std::vector<ll> &tris,
      const std::vector<double> &thickness, double cellSize, ll maxBins)
  {
    points = &pts;
    triangles = &tris;
    double hi[3];
    for (int k = 0; k < 3; k++)
    {
      lo[k] = std::numeric_limits<double>::max();
      hi[k] = -std::numeric_limits<double>::max();
    }
    for (size_t i = 0; i + 2 < pts.size(); i += 3)
    {
      for (int k = 0; k < 3; k++)
      {
        lo[k] = std::min(lo[k], pts[i + k]);
        hi[k] = std::max(hi[k], pts[i + k]);
      }
    }
    double volume = 1.0;
    for (int k = 0; k < 3; k++)
    {
      volume *= std::max(hi[k] - lo[k], 1.0e-12) + 2.0*cellSize;
    }
    cell = std::max(cellSize, std::cbrt(volume/(double)maxBins));
    for (int k = 0; k < 3; k++)
    {
      lo[k] -= cell;
      n[k] = std::max(1, (int)std::ceil((hi[k] - lo[k])/cell) + 2);
    }
    size_t numBins = (size_t)n[0]*n[1]*n[2];
    start.assign(numBins + 1, 0);
    auto range = [&](size_t t, int i0[3], int i1[3])
    {
      double blo[3], bhi[3];
      for (int k = 0; k < 3; k++)
      {
        blo[k] = std::numeric_limits<double>::max();
        bhi[k] = -std::numeric_limits<double>::max();
      }
      for (int j = 0; j < 3; j++)
      {
        const double *p = &pts[(size_t)3*tris[3*t + j]];
        for (int k = 0; k < 3; k++)
        {
          blo[k] = std::min(blo[k], p[k]);
          bhi[k] = std::max(bhi[k], p[k]);
        }
      }
      for (int k = 0; k < 3; k++)
      {
        i0[k] = Bin(blo[k], k);
        i1[k] = Bin(bhi[k], k);
      }
    };
    size_t numTris = tris.size()/3;
    for (size_t t = 0; t < numTris; t++)
    {
      int i0[3], i1[3];
      range(t, i0, i1);
      for (int k = i0[2]; k <= i1[2]; k++)
        for (int j = i0[1]; j <= i1[1]; j++)
          for (int i = i0[0]; i <= i1[0]; i++)
            start[Index(i, j, k) + 1]++;
    }
    for (size_t b = 0; b < numBins; b++)
    {
      start[b + 1] += start[b];
    }
    cells.assign((size_t)start[numBins], -1);
    std::vector<ll> cursor(start.begin(), start.end() - 1);
    for (size_t t = 0; t < numTris; t++)
    {
      int i0[3], i1[3];
      range(t, i0, i1);
      for (int k = i0[2]; k <= i1[2]; k++)
        for (int j = i0[1]; j <= i1[1]; j++)
          for (int i = i0[0]; i <= i1[0]; i++)
            cells[(size_t)cursor[Index(i, j, k)]++] = (ll)t;
    }
    binMaxThickness.assign(numBins, 0.0);
    for (size_t b = 0; b < numBins; b++)
    {
      for (ll m = start[b]; m < start[b + 1]; m++)
      {
        const ll *tt = &tris[(size_t)3*cells[(size_t)m]];
        for (int j = 0; j < 3; j++)
        {
          binMaxThickness[b] = std::max(binMaxThickness[b], thickness[(size_t)tt[j]]);
        }
      }
    }
  }

  // The distance from x to the box of bin (i, j, k).
  double BoxDistance(const double x[3], int i, int j, int k) const
  {
    int idx[3] = {i, j, k};
    double dd = 0.0;
    for (int a = 0; a < 3; a++)
    {
      double blo = lo[a] + idx[a]*cell, bhi = blo + cell;
      double e = (x[a] < blo) ? blo - x[a] : ((x[a] > bhi) ? x[a] - bhi : 0.0);
      dd += e*e;
    }
    return std::sqrt(dd);
  }
};

struct EdgeKey
{
  ll a, b;
  EdgeKey(ll p, ll q) : a(std::min(p, q)), b(std::max(p, q)) {}
  bool operator<(const EdgeKey &o) const
  {
    return (a != o.a) ? (a < o.a) : (b < o.b);
  }
};

void NoteFault(Report &report, const char *what, const double at[3])
{
  if (report.firstFault.empty())
  {
    report.firstFault = what;
    for (int k = 0; k < 3; k++)
    {
      report.firstFaultAt[k] = at[k];
    }
  }
}

}  // namespace

//---------------------
// CountEdges
//---------------------

// -------------------------------------
// TrimSurfaceAtCaps
// -------------------------------------
int TrimSurfaceAtCaps(Surface &surface, const std::vector<CapPlane> &planes, double snapFraction,
    TrimReport &report, std::string &error)
{
  std::vector<double> &pts = surface.points;
  std::vector<ll> &tris = surface.triangles;
  std::vector<ll> &owner = surface.pointRim;
  ll np = (ll)(pts.size()/3), nt = (ll)(tris.size()/3);
  if (owner.size() != (size_t)np)
  {
    error = "the surface carries no ownership of its points, so it cannot be trimmed by cap";
    return 1;
  }
  for (size_t i = 0; i < tris.size(); i++)
  {
    if (tris[i] < 0 || tris[i] >= np)
    {
      error = "a triangle of the surface refers to a point it does not have";
      return 1;
    }
  }
  report = TrimReport();
  report.numPointsBefore = np;
  report.numTrianglesBefore = nt;
  report.numCut.assign(planes.size(), 0);
  report.numSnapped.assign(planes.size(), 0);

  // The mean edge at each point, which sizes the band a point is snapped
  // onto the plane from.
  std::vector<double> meanEdge((size_t)np, 0.0);
  std::vector<int> degree((size_t)np, 0);
  for (ll t = 0; t < nt; t++)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3];
      double L = Distance(&pts[(size_t)3*a], &pts[(size_t)3*b]);
      meanEdge[(size_t)a] += L;
      degree[(size_t)a]++;
      meanEdge[(size_t)b] += L;
      degree[(size_t)b]++;
    }
  }
  for (ll v = 0; v < np; v++)
  {
    meanEdge[(size_t)v] = (degree[(size_t)v] > 0) ? meanEdge[(size_t)v]/(double)degree[(size_t)v] : 0.0;
  }

  for (size_t p = 0; p < planes.size(); p++)
  {
    const CapPlane &plane = planes[p];
    np = (ll)(pts.size()/3);
    nt = (ll)(tris.size()/3);
    // The signed height of each point over the plane, positive on the side
    // that is kept. A point of this plane's rim within the snap band is
    // moved onto the plane. A point that belongs to no plane or to another
    // is kept whatever side it lies on: its height is its own where that is
    // comfortably positive, and a mean edge otherwise, so that an edge from
    // it to a point coming off is cut somewhere along itself rather than at
    // infinity.
    std::vector<double> s((size_t)np, 0.0);
    for (ll v = 0; v < np; v++)
    {
      double off[3];
      Sub(plane.origin, &pts[(size_t)3*v], off);
      double sv = Dot(off, plane.outward);
      if (owner[(size_t)v] == plane.rim)
      {
        if (std::abs(sv) <= snapFraction*meanEdge[(size_t)v])
        {
          for (int k = 0; k < 3; k++)
          {
            pts[(size_t)3*v + k] += sv*plane.outward[k];
          }
          sv = 0.0;
          report.numSnapped[p]++;
        }
        if (sv < 0.0)
        {
          report.numCut[p]++;
        }
      }
      else
      {
        sv = std::max(sv, meanEdge[(size_t)v]);
      }
      s[(size_t)v] = sv;
    }

    std::vector<ll> kept;
    kept.reserve(tris.size());
    std::map<EdgeKey, ll> cutPoint;
    // The point where the plane crosses the edge v-w (heights of opposite
    // sign), made once per edge and put exactly on the plane.
    auto cutOn = [&](ll v, ll w) -> ll
    {
      EdgeKey key(v, w);
      std::map<EdgeKey, ll>::iterator it = cutPoint.find(key);
      if (it != cutPoint.end())
      {
        return it->second;
      }
      double f = s[(size_t)v]/(s[(size_t)v] - s[(size_t)w]);
      double q[3];
      for (int k = 0; k < 3; k++)
      {
        q[k] = pts[(size_t)3*v + k] + f*(pts[(size_t)3*w + k] - pts[(size_t)3*v + k]);
      }
      double off[3];
      Sub(plane.origin, q, off);
      double sq = Dot(off, plane.outward);
      for (int k = 0; k < 3; k++)
      {
        q[k] += sq*plane.outward[k];
      }
      ll id = (ll)(pts.size()/3);
      pts.insert(pts.end(), q, q + 3);
      owner.push_back(-1);
      meanEdge.push_back(0.5*(meanEdge[(size_t)v] + meanEdge[(size_t)w]));
      s.push_back(0.0);
      cutPoint[key] = id;
      return id;
    };

    for (ll t = 0; t < nt; t++)
    {
      ll v[3] = {tris[(size_t)3*t], tris[(size_t)3*t + 1], tris[(size_t)3*t + 2]};
      int numBelow = 0, numOn = 0;
      for (int j = 0; j < 3; j++)
      {
        if (s[(size_t)v[j]] < 0.0) numBelow++;
        else if (s[(size_t)v[j]] == 0.0) numOn++;
      }
      if (numBelow == 0)
      {
        if (numOn == 3)
        {
          // flat in the plane, where the annulus will be
          report.numDropped++;
          continue;
        }
        kept.insert(kept.end(), v, v + 3);
        continue;
      }
      if (numBelow == 3)
      {
        report.numRemoved++;
        continue;
      }
      // The part of the triangle on the kept side, walked in its own order:
      // a corner on or above the plane stays, an edge from above to below
      // (or back) gets its crossing point, and an edge from a corner on the
      // plane to one below gets nothing, since that corner is the crossing.
      ll poly[4];
      int n = 0;
      for (int j = 0; j < 3; j++)
      {
        ll a = v[j], b = v[(j+1)%3];
        double sa = s[(size_t)a], sb = s[(size_t)b];
        if (sa >= 0.0)
        {
          poly[n++] = a;
        }
        if ((sa > 0.0 && sb < 0.0) || (sa < 0.0 && sb > 0.0))
        {
          poly[n++] = cutOn(a, b);
        }
      }
      if (n < 3)
      {
        // only a point or an edge in the plane was on the kept side
        report.numDropped++;
        continue;
      }
      report.numSplit++;
      kept.push_back(poly[0]);
      kept.push_back(poly[1]);
      kept.push_back(poly[2]);
      if (n == 4)
      {
        kept.push_back(poly[0]);
        kept.push_back(poly[2]);
        kept.push_back(poly[3]);
      }
    }
    tris.swap(kept);
  }

  // Compact the points to those the triangles use, in their order.
  np = (ll)(pts.size()/3);
  {
    std::vector<ll> newId((size_t)np, -1);
    ll m = 0;
    for (size_t i = 0; i < tris.size(); i++)
    {
      if (newId[(size_t)tris[i]] < 0) newId[(size_t)tris[i]] = 0;
    }
    for (ll v = 0; v < np; v++)
    {
      if (newId[(size_t)v] == 0) newId[(size_t)v] = m++;
    }
    std::vector<double> cpts((size_t)3*m);
    std::vector<ll> cowner((size_t)m);
    std::vector<double> cedge((size_t)m);
    for (ll v = 0; v < np; v++)
    {
      ll u = newId[(size_t)v];
      if (u < 0) continue;
      for (int k = 0; k < 3; k++)
      {
        cpts[(size_t)3*u + k] = pts[(size_t)3*v + k];
      }
      cowner[(size_t)u] = owner[(size_t)v];
      cedge[(size_t)u] = meanEdge[(size_t)v];
    }
    for (size_t i = 0; i < tris.size(); i++)
    {
      tris[i] = newId[(size_t)tris[i]];
    }
    pts.swap(cpts);
    owner.swap(cowner);
    meanEdge.swap(cedge);
    np = m;
  }
  nt = (ll)(tris.size()/3);

  // Rim edges much shorter than the surface's own edges are collapsed: the
  // trim leaves them where the plane cut a triangle a hair from one of its
  // corners, or where a snapped point landed next to a cut point, and the
  // shell's annulus, which walks the rim by angle about the cap's centre,
  // folds a triangle over its neighbour on such a tooth (measured
  // 2026-09-23 on the user's 178k model with three layers: rim edges of 2
  // to 8 percent of the mean, and TetGen refused two coplanar annulus
  // facets on one of them). The collapse takes the edge's end with fewer
  // triangles onto the other, on the rim, by less than a quarter of an
  // edge, which is within what the decimation moves anyway; the link
  // condition of a boundary edge (the two ends share only the third corner
  // of their one triangle) keeps the surface a manifold.
  {
    const double shortFraction = 0.25;
    std::vector<unsigned char> deadTri((size_t)nt, 0), deadPt((size_t)np, 0);
    for (int pass = 0; pass < 4; pass++)
    {
      std::map<EdgeKey, int> edgeCount;
      std::vector<std::vector<ll> > incident((size_t)np);
      for (ll t = 0; t < nt; t++)
      {
        if (deadTri[(size_t)t]) continue;
        for (int j = 0; j < 3; j++)
        {
          edgeCount[EdgeKey(tris[(size_t)3*t + j], tris[(size_t)3*t + (j+1)%3])]++;
          incident[(size_t)tris[(size_t)3*t + j]].push_back(t);
        }
      }
      ll merged = 0;
      std::vector<unsigned char> touched((size_t)np, 0);
      for (ll t = 0; t < nt; t++)
      {
        if (deadTri[(size_t)t]) continue;
        for (int j = 0; j < 3; j++)
        {
          ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3], c = tris[(size_t)3*t + (j+2)%3];
          if (deadTri[(size_t)t] || touched[(size_t)a] || touched[(size_t)b] || touched[(size_t)c]) continue;
          if (edgeCount[EdgeKey(a, b)] != 1) continue;
          double mean = 0.5*(meanEdge[(size_t)a] + meanEdge[(size_t)b]);
          if (!(mean > 0.0) || Distance(&pts[(size_t)3*a], &pts[(size_t)3*b]) >= shortFraction*mean) continue;
          // the end with fewer triangles goes onto the other
          ll u = a, v = b;
          if (incident[(size_t)a].size() > incident[(size_t)b].size()) { u = b; v = a; }
          // the link condition: the neighbours of u and v share only c
          bool ok = true;
          for (size_t m = 0; m < incident[(size_t)u].size() && ok; m++)
          {
            ll tu = incident[(size_t)u][m];
            if (deadTri[(size_t)tu] || tu == t) continue;
            for (int q = 0; q < 3 && ok; q++)
            {
              ll w = tris[(size_t)3*tu + q];
              if (w == u || w == c) continue;
              for (size_t n = 0; n < incident[(size_t)v].size() && ok; n++)
              {
                ll tv = incident[(size_t)v][n];
                if (deadTri[(size_t)tv] || tv == t) continue;
                for (int r = 0; r < 3; r++)
                {
                  if (tris[(size_t)3*tv + r] == w) { ok = false; break; }
                }
              }
            }
          }
          if (!ok) continue;
          for (size_t m = 0; m < incident[(size_t)u].size(); m++)
          {
            ll tu = incident[(size_t)u][m];
            if (deadTri[(size_t)tu]) continue;
            if (tu == t) { deadTri[(size_t)tu] = 1; continue; }
            for (int q = 0; q < 3; q++)
            {
              if (tris[(size_t)3*tu + q] == u) tris[(size_t)3*tu + q] = v;
            }
            incident[(size_t)v].push_back(tu);
          }
          incident[(size_t)u].clear();
          deadPt[(size_t)u] = 1;
          touched[(size_t)u] = 1;
          touched[(size_t)v] = 1;
          touched[(size_t)c] = 1;
          merged++;
          break;   // this triangle is gone
        }
      }
      report.numRimEdgesMerged += merged;
      if (merged == 0) break;
    }
    if (report.numRimEdgesMerged > 0)
    {
      std::vector<ll> newId((size_t)np, -1);
      std::vector<double> cpts;
      std::vector<ll> cowner;
      std::vector<double> cedge;
      std::vector<ll> ctris;
      for (ll t = 0; t < nt; t++)
      {
        if (deadTri[(size_t)t]) continue;
        for (int j = 0; j < 3; j++)
        {
          ll v = tris[(size_t)3*t + j];
          if (newId[(size_t)v] < 0)
          {
            newId[(size_t)v] = (ll)(cpts.size()/3);
            cpts.insert(cpts.end(), &pts[(size_t)3*v], &pts[(size_t)3*v] + 3);
            cowner.push_back(owner[(size_t)v]);
            cedge.push_back(meanEdge[(size_t)v]);
          }
          ctris.push_back(newId[(size_t)v]);
        }
      }
      pts.swap(cpts);
      owner.swap(cowner);
      meanEdge.swap(cedge);
      tris.swap(ctris);
      np = (ll)(pts.size()/3);
      nt = (ll)(tris.size()/3);
    }
  }
  report.numPointsAfter = np;
  report.numTrianglesAfter = nt;

  // The rims: their shortest edge against the mean edge at its ends, and the
  // angle each rim triangle makes with its cap plane on the side the annulus
  // will lie, toward the plane's origin. A triangle leaning over the annulus
  // at a hair's angle is what TetGen refuses as two facets folded onto each
  // other.
  {
    std::map<EdgeKey, int> edgeCount;
    for (ll t = 0; t < nt; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        edgeCount[EdgeKey(tris[(size_t)3*t + j], tris[(size_t)3*t + (j+1)%3])]++;
      }
    }
    report.shortestRimEdgeRatio = std::numeric_limits<double>::max();
    report.smallestRimAngleDegrees = 180.0;
    for (ll t = 0; t < nt; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3], c = tris[(size_t)3*t + (j+2)%3];
        if (edgeCount[EdgeKey(a, b)] != 1) continue;
        report.numRimEdges++;
        const double *xa = &pts[(size_t)3*a], *xb = &pts[(size_t)3*b], *xc = &pts[(size_t)3*c];
        double L = Distance(xa, xb);
        double mean = 0.5*(meanEdge[(size_t)a] + meanEdge[(size_t)b]);
        if (mean > 0.0)
        {
          report.shortestRimEdgeRatio = std::min(report.shortestRimEdgeRatio, L/mean);
        }
        // the plane this edge lies on: both ends on it
        for (size_t p = 0; p < planes.size(); p++)
        {
          const CapPlane &plane = planes[p];
          double offA[3], offB[3], offC[3];
          Sub(plane.origin, xa, offA);
          Sub(plane.origin, xb, offB);
          Sub(plane.origin, xc, offC);
          double tol = 1.0e-6*std::max(mean, 1.0e-12);
          if (std::abs(Dot(offA, plane.outward)) > tol || std::abs(Dot(offB, plane.outward)) > tol) continue;
          double height = Dot(offC, plane.outward);   // positive on the kept side
          // inward in the plane: from the edge's midpoint toward the origin, perpendicular to the edge
          double mid[3], toOrigin[3], along[3];
          for (int k = 0; k < 3; k++)
          {
            mid[k] = 0.5*(xa[k] + xb[k]);
            along[k] = xb[k] - xa[k];
          }
          Sub(plane.origin, mid, toOrigin);
          double tn = Dot(toOrigin, plane.outward);
          for (int k = 0; k < 3; k++) toOrigin[k] -= tn*plane.outward[k];
          if (!Normalize(along)) break;
          double ta = Dot(toOrigin, along);
          for (int k = 0; k < 3; k++) toOrigin[k] -= ta*along[k];
          if (!Normalize(toOrigin)) break;
          double lean[3];
          Sub(xc, mid, lean);
          double inward = Dot(lean, toOrigin);
          double angle = std::atan2(std::abs(height), inward)*(180.0/3.14159265358979323846);
          report.smallestRimAngleDegrees = std::min(report.smallestRimAngleDegrees, angle);
          break;
        }
      }
    }
    if (report.numRimEdges == 0)
    {
      report.shortestRimEdgeRatio = 0.0;
      report.smallestRimAngleDegrees = 0.0;
    }
  }
  return 0;
}

void CountEdges(const std::vector<long long> &triangles, long long &numBoundary,
    long long &numNonManifold, long long &numMiswound)
{
  numBoundary = numNonManifold = numMiswound = 0;
  std::map<EdgeKey, std::pair<int, int> > use;   // (count, forward count)
  for (size_t i = 0; i + 2 < triangles.size(); i += 3)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = triangles[i + j], b = triangles[i + (j+1)%3];
      std::pair<int, int> &u = use[EdgeKey(a, b)];
      u.first++;
      if (a < b)
      {
        u.second++;
      }
    }
  }
  for (std::map<EdgeKey, std::pair<int, int> >::iterator it = use.begin(); it != use.end(); ++it)
  {
    if (it->second.first == 1)
    {
      numBoundary++;
    }
    else if (it->second.first > 2)
    {
      numNonManifold++;
    }
    else if (it->second.second != 1)
    {
      numMiswound++;
    }
  }
}

//---------------------
// OffsetField
//---------------------

struct OffsetField::Data
{
  std::vector<double> points, normals, thickness;
  std::vector<ll> triangles;
  std::vector<double> faceNormals;   // unit, per triangle
  std::vector<double> localSize;     // per triangle: the target edge length of the offset surface there
  std::vector<double> localThickness;   // per triangle: the mean thickness of its corners
  std::vector<unsigned char> boundaryPoint;   // per point: 1 on an open edge of the field's surface (the collar ends)
  std::vector<ll> triangleRim;       // per triangle: the cap rim it belongs to (its collar, or the interface within two collar lengths of the rim along the interface), -1 for none
  std::vector<std::vector<ll> > rims;
  std::vector<double> rimOutward;    // three per rim
  std::vector<double> collarLength;  // per rim
  ll numInterfaceTriangles = 0;
  ll numInterfacePoints = 0;
  double reach = 0.0;
  double largestThickness = 0.0;
  double smallestThickness = 0.0;
  double meanEdge = 0.0;
  TriangleGrid grid;
  mutable ll numEvaluations = 0;
};

OffsetField::OffsetField() : data_(new Data()) {}

OffsetField::~OffsetField()
{
  delete data_;
}

int OffsetField::Build(const Interface &input, Report &report, std::string &error, double chordTolerance)
{
  Data &d = *data_;
  d = Data();
  ll numPts = (ll)(input.points.size()/3);
  ll numTris = (ll)(input.triangles.size()/3);
  if (input.points.size() != (size_t)3*numPts || input.triangles.size() != (size_t)3*numTris ||
      input.normals.size() != (size_t)3*numPts || input.thickness.size() != (size_t)numPts)
  {
    error = "the interface arrays are not three per point, three per triangle and one thickness per point";
    return 1;
  }
  if (numTris == 0)
  {
    error = "the interface has no triangles";
    return 1;
  }
  for (size_t m = 0; m < input.triangles.size(); m++)
  {
    if (input.triangles[m] < 0 || input.triangles[m] >= numPts)
    {
      error = "a triangle refers to a point outside the point list";
      return 1;
    }
  }
  d.points = input.points;
  d.normals = input.normals;
  d.thickness = input.thickness;
  d.triangles = input.triangles;
  d.numInterfacePoints = numPts;
  d.numInterfaceTriangles = numTris;
  d.smallestThickness = std::numeric_limits<double>::max();
  for (ll i = 0; i < numPts; i++)
  {
    double *n = &d.normals[(size_t)3*i];
    if (!Normalize(n))
    {
      char what[160];
      snprintf(what, sizeof(what), "the normal at point %lld has no length", i);
      error = what;
      return 1;
    }
    double t = d.thickness[(size_t)i];
    if (!(t > 0.0) || !std::isfinite(t))
    {
      char what[160];
      snprintf(what, sizeof(what), "the thickness at point %lld is %g; every point needs a positive thickness", i, t);
      error = what;
      return 1;
    }
    d.largestThickness = std::max(d.largestThickness, t);
    d.smallestThickness = std::min(d.smallestThickness, t);
  }
  report.numInterfacePoints = numPts;
  report.numInterfaceTriangles = numTris;
  report.smallestThickness = d.smallestThickness;
  report.largestThickness = d.largestThickness;

  // The mean edge, for the grid and the collars.
  double edgeSum = 0.0;
  for (ll t = 0; t < numTris; t++)
  {
    for (int j = 0; j < 3; j++)
    {
      edgeSum += Distance(&d.points[(size_t)3*d.triangles[(size_t)3*t + j]],
          &d.points[(size_t)3*d.triangles[(size_t)3*t + (j+1)%3]]);
    }
  }
  d.meanEdge = edgeSum/(3.0*(double)numTris);

  // The cap rims: the boundary edges, each traversed as its triangle does,
  // chained into loops.
  {
    std::map<EdgeKey, int> edgeCount;
    for (ll t = 0; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        edgeCount[EdgeKey(d.triangles[(size_t)3*t + j], d.triangles[(size_t)3*t + (j+1)%3])]++;
      }
    }
    std::map<ll, ll> next;
    for (ll t = 0; t < numTris; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = d.triangles[(size_t)3*t + j], b = d.triangles[(size_t)3*t + (j+1)%3];
        std::map<EdgeKey, int>::iterator it = edgeCount.find(EdgeKey(a, b));
        if (it->second == 1)
        {
          next[a] = b;
        }
        else if (it->second > 2)
        {
          error = "the interface has an edge on more than two triangles";
          return 1;
        }
      }
    }
    std::map<ll, unsigned char> seen;
    for (std::map<ll, ll>::iterator it = next.begin(); it != next.end(); ++it)
    {
      if (seen.count(it->first))
      {
        continue;
      }
      std::vector<ll> loop;
      ll cur = it->first;
      while (!seen.count(cur))
      {
        seen[cur] = 1;
        loop.push_back(cur);
        std::map<ll, ll>::iterator nx = next.find(cur);
        if (nx == next.end())
        {
          break;
        }
        cur = nx->second;
      }
      if (loop.size() >= 3)
      {
        d.rims.push_back(loop);
      }
    }
  }
  report.numRims = (ll)d.rims.size();

  // A collar past each rim: the rim pushed along the outward direction of its
  // plane by the wall thickness there plus an edge, carrying the rim's normals
  // and thickness, so that the offset is a straight tube through the cap
  // plane rather than rounding off around the rim. The rim traversed in its
  // triangles' winding runs clockwise about the outward direction, so the
  // outward direction is the reverse of the loop's own normal; that is
  // checked against the interface itself, which lies on the inward side.
  // The interface's edges by point, for measuring along the interface how
  // far from a rim a triangle lies.
  std::vector<std::vector<std::pair<ll, double> > > adjacent((size_t)numPts);
  for (ll t = 0; t < numTris; t++)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = d.triangles[(size_t)3*t + j], b = d.triangles[(size_t)3*t + (j+1)%3];
      double L = Distance(&d.points[(size_t)3*a], &d.points[(size_t)3*b]);
      adjacent[(size_t)a].push_back(std::make_pair(b, L));
      adjacent[(size_t)b].push_back(std::make_pair(a, L));
    }
  }
  const double unreached = std::numeric_limits<double>::max();
  std::vector<double> alongInterface((size_t)numPts, unreached);
  d.rimOutward.assign(3*d.rims.size(), 0.0);
  d.collarLength.assign(d.rims.size(), 0.0);
  for (size_t r = 0; r < d.rims.size(); r++)
  {
    const std::vector<ll> &loop = d.rims[r];
    double centre[3] = {0.0, 0.0, 0.0};
    for (size_t m = 0; m < loop.size(); m++)
    {
      for (int k = 0; k < 3; k++)
      {
        centre[k] += d.points[(size_t)3*loop[m] + k]/(double)loop.size();
      }
    }
    double normal[3] = {0.0, 0.0, 0.0};
    double rimEdge = 0.0, capThickness = 0.0, radius = 0.0;
    for (size_t m = 0; m < loop.size(); m++)
    {
      const double *p = &d.points[(size_t)3*loop[m]], *q = &d.points[(size_t)3*loop[(m+1)%loop.size()]];
      normal[0] += (p[1] - q[1])*(p[2] + q[2]);
      normal[1] += (p[2] - q[2])*(p[0] + q[0]);
      normal[2] += (p[0] - q[0])*(p[1] + q[1]);
      rimEdge += Distance(p, q)/(double)loop.size();
      capThickness = std::max(capThickness, d.thickness[(size_t)loop[m]]);
      radius = std::max(radius, Distance(p, centre));
    }
    if (!Normalize(normal))
    {
      char what[160];
      snprintf(what, sizeof(what), "a cap rim of %zu points at (%.5g, %.5g, %.5g) encloses no area", loop.size(), centre[0], centre[1], centre[2]);
      error = what;
      return 1;
    }
    double outward[3] = {-normal[0], -normal[1], -normal[2]};
    // the interface near the rim should lie on the inward side
    double inward = 0.0, outwardSide = 0.0;
    for (ll i = 0; i < numPts; i++)
    {
      double off[3];
      Sub(&d.points[(size_t)3*i], centre, off);
      if (Norm(off) > 3.0*radius)
      {
        continue;
      }
      double s = Dot(off, outward);
      if (s > 0.05*radius) outwardSide += 1.0;
      else if (s < -0.05*radius) inward += 1.0;
    }
    if (outwardSide > inward)
    {
      for (int k = 0; k < 3; k++)
      {
        outward[k] = -outward[k];
      }
    }
    for (int k = 0; k < 3; k++)
    {
      d.rimOutward[3*r + k] = outward[k];
    }
    double length = capThickness + rimEdge;
    d.collarLength[r] = length;
    ll base = (ll)(d.points.size()/3);
    for (size_t m = 0; m < loop.size(); m++)
    {
      ll v = loop[m];
      for (int k = 0; k < 3; k++)
      {
        d.points.push_back(d.points[(size_t)3*v + k] + length*outward[k]);
      }
      for (int k = 0; k < 3; k++)
      {
        d.normals.push_back(d.normals[(size_t)3*v + k]);
      }
      d.thickness.push_back(d.thickness[(size_t)v]);
    }
    // The rim edge a->b is traversed a->b by its interface triangle, so the
    // collar triangles on it traverse it b->a.
    for (size_t m = 0; m < loop.size(); m++)
    {
      ll a = loop[m], b = loop[(m+1)%loop.size()];
      ll a2 = base + (ll)m, b2 = base + (ll)((m+1)%loop.size());
      d.triangles.push_back(b);
      d.triangles.push_back(a);
      d.triangles.push_back(a2);
      d.triangles.push_back(a2);
      d.triangles.push_back(b2);
      d.triangles.push_back(b);
      report.numCollarTriangles += 2;
    }
    // The interface within two collar lengths of the rim belongs to the rim
    // too, so that a point of the offset surface just inside the cap plane
    // is owned like one just past it and a trim along the plane cuts the
    // edges between them on the plane. The distance is measured along the
    // interface, not within a cylinder about the cap: another vessel passing
    // close by the cap end falls in such a cylinder, its wall is handed to
    // this cap, and the trim then cuts a hole in it past the cap plane (two
    // vessels 0.5 apart on the user's model did exactly that).
    if (d.triangleRim.size() < (size_t)numTris)
    {
      d.triangleRim.assign((size_t)numTris, -1);
    }
    {
      const double reachAlong = 2.0*length;
      std::vector<ll> touched;
      std::priority_queue<std::pair<double, ll>, std::vector<std::pair<double, ll> >, std::greater<std::pair<double, ll> > > queue;
      for (size_t m = 0; m < loop.size(); m++)
      {
        alongInterface[(size_t)loop[m]] = 0.0;
        touched.push_back(loop[m]);
        queue.push(std::make_pair(0.0, loop[m]));
      }
      while (!queue.empty())
      {
        std::pair<double, ll> top = queue.top();
        queue.pop();
        if (top.first > alongInterface[(size_t)top.second]) continue;
        const std::vector<std::pair<ll, double> > &around = adjacent[(size_t)top.second];
        for (size_t m = 0; m < around.size(); m++)
        {
          ll v = around[m].first;
          double cand = top.first + around[m].second;
          if (cand > reachAlong || cand >= alongInterface[(size_t)v]) continue;
          if (alongInterface[(size_t)v] == unreached) touched.push_back(v);
          alongInterface[(size_t)v] = cand;
          queue.push(std::make_pair(cand, v));
        }
      }
      for (ll t = 0; t < numTris; t++)
      {
        if (d.triangleRim[(size_t)t] >= 0) continue;
        for (int j = 0; j < 3; j++)
        {
          if (alongInterface[(size_t)d.triangles[(size_t)3*t + j]] != unreached)
          {
            d.triangleRim[(size_t)t] = (ll)r;
            break;
          }
        }
      }
      for (size_t m = 0; m < touched.size(); m++)
      {
        alongInterface[(size_t)touched[m]] = unreached;
      }
    }
    for (size_t m = 0; m < 2*loop.size(); m++)
    {
      d.triangleRim.push_back((ll)r);
    }
  }
  if (d.triangleRim.size() < d.triangles.size()/3)
  {
    d.triangleRim.resize(d.triangles.size()/3, -1);
  }

  // The open edges of the field's surface - the collar ends - and their
  // points. Beyond them there is no lumen side: a point whose closest point
  // lies on one is outside the wall by its distance, whichever way the
  // normal there happens to lean.
  ll numF = (ll)(d.triangles.size()/3);
  d.boundaryPoint.assign(d.points.size()/3, 0);
  {
    std::map<EdgeKey, int> edgeCount;
    for (ll t = 0; t < numF; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        edgeCount[EdgeKey(d.triangles[(size_t)3*t + j], d.triangles[(size_t)3*t + (j+1)%3])]++;
      }
    }
    for (std::map<EdgeKey, int>::iterator it = edgeCount.begin(); it != edgeCount.end(); ++it)
    {
      if (it->second == 1)
      {
        d.boundaryPoint[(size_t)it->first.a] = 1;
        d.boundaryPoint[(size_t)it->first.b] = 1;
      }
    }
  }

  // Face normals, the local size and the grid.
  d.faceNormals.assign((size_t)3*numF, 0.0);
  d.localSize.assign((size_t)numF, d.meanEdge);
  d.localThickness.assign((size_t)numF, d.largestThickness);
  for (ll t = 0; t < numF; t++)
  {
    const ll *tt = &d.triangles[(size_t)3*t];
    const double *a = &d.points[(size_t)3*tt[0]], *b = &d.points[(size_t)3*tt[1]], *c = &d.points[(size_t)3*tt[2]];
    double e1[3], e2[3];
    Sub(b, a, e1);
    Sub(c, a, e2);
    Cross(e1, e2, &d.faceNormals[(size_t)3*t]);
    Normalize(&d.faceNormals[(size_t)3*t]);
    // The target edge of the offset surface here: the interface's own edge,
    // and no longer than keeps the chord of a surface with the curvature of
    // the interface plus the wall within a twentieth of the wall. The curvature
    // is read off the turning of the normals along the edges; a collar's
    // edges along the outward direction carry no turning and drop out.
    double meanEdge = 0.0, radius = 1.0e3, thickness = 0.0;
    for (int j = 0; j < 3; j++)
    {
      ll u = tt[j], v = tt[(j+1)%3];
      double L = Distance(&d.points[(size_t)3*u], &d.points[(size_t)3*v]);
      meanEdge += L/3.0;
      thickness += d.thickness[(size_t)u]/3.0;
      double dn[3];
      Sub(&d.normals[(size_t)3*u], &d.normals[(size_t)3*v], dn);
      double turn = Norm(dn);
      if (turn > 1.0e-9)
      {
        radius = std::min(radius, L/turn);
      }
    }
    double chord = std::sqrt(8.0*(radius + thickness)*chordTolerance*thickness);
    d.localSize[(size_t)t] = std::max(std::min(meanEdge, chord), 1.0e-6*d.meanEdge);
    d.localThickness[(size_t)t] = thickness;
  }
  d.reach = d.largestThickness + 2.0*d.meanEdge;
  d.grid.Build(d.points, d.triangles, d.thickness, std::max(d.meanEdge, 0.25*d.largestThickness), 4000000);
  return 0;
}

double OffsetField::Evaluate(const double x[3], double thicknessScale) const
{
  const Data &d = *data_;
  d.numEvaluations++;
  const TriangleGrid &g = d.grid;
  int c[3];
  for (int k = 0; k < 3; k++)
  {
    c[k] = g.Bin(x[k], k);
  }
  double best = std::numeric_limits<double>::max();
  ll bestT = -1;
  double bestQ[3] = {0.0, 0.0, 0.0}, bestBary[3] = {0.0, 0.0, 0.0};
  double bestTerm = d.reach;
  int rings = (int)(d.reach/g.cell) + 2;
  for (int r = 0; r <= rings; r++)
  {
    // every bin of this ring is at least (r-1) cells away
    double bound = (r - 1)*g.cell;
    if (bound - d.largestThickness >= bestTerm && bound >= best)
    {
      break;
    }
    for (int k = c[2] - r; k <= c[2] + r; k++)
    {
      for (int j = c[1] - r; j <= c[1] + r; j++)
      {
        for (int i = c[0] - r; i <= c[0] + r; i++)
        {
          if (std::max(std::abs(k - c[2]), std::max(std::abs(j - c[1]), std::abs(i - c[0]))) != r) continue;
          if (i < 0 || j < 0 || k < 0 || i >= g.n[0] || j >= g.n[1] || k >= g.n[2]) continue;
          size_t b = g.Index(i, j, k);
          if (g.start[b + 1] == g.start[b]) continue;
          double dd = g.BoxDistance(x, i, j, k);
          if (dd - g.binMaxThickness[b] >= bestTerm && dd >= best) continue;
          for (ll m = g.start[b]; m < g.start[b + 1]; m++)
          {
            ll t = g.cells[(size_t)m];
            const ll *tt = &d.triangles[(size_t)3*t];
            double q[3], bary[3];
            double dist = ClosestOnTriangle(x, &d.points[(size_t)3*tt[0]], &d.points[(size_t)3*tt[1]], &d.points[(size_t)3*tt[2]], q, bary);
            if (dist < best)
            {
              best = dist;
              bestT = t;
              for (int a = 0; a < 3; a++)
              {
                bestQ[a] = q[a];
                bestBary[a] = bary[a];
              }
            }
            double tq = bary[0]*d.thickness[(size_t)tt[0]] + bary[1]*d.thickness[(size_t)tt[1]] + bary[2]*d.thickness[(size_t)tt[2]];
            bestTerm = std::min(bestTerm, dist - thicknessScale*tq);
          }
        }
      }
    }
  }
  if (bestT < 0)
  {
    return d.reach;
  }
  // Which side of the surface x is on: by the normal at the closest point of
  // the nearest triangle - the face normal inside it, the corners' normals
  // blended on an edge or at a corner.
  const ll *tt = &d.triangles[(size_t)3*bestT];
  double n[3];
  if (bestBary[0] > 1.0e-6 && bestBary[1] > 1.0e-6 && bestBary[2] > 1.0e-6)
  {
    for (int k = 0; k < 3; k++)
    {
      n[k] = d.faceNormals[(size_t)3*bestT + k];
    }
  }
  else
  {
    for (int k = 0; k < 3; k++)
    {
      n[k] = bestBary[0]*d.normals[(size_t)3*tt[0] + k] + bestBary[1]*d.normals[(size_t)3*tt[1] + k] + bestBary[2]*d.normals[(size_t)3*tt[2] + k];
    }
  }
  // A closest point on an open edge of the surface - at a boundary corner,
  // or on an edge between two boundary corners - has no lumen behind it.
  bool onOpenEdge = false;
  {
    const double tol = 1.0e-6;
    int numZero = 0, zeroAt = -1, oneAt = -1;
    for (int j = 0; j < 3; j++)
    {
      if (bestBary[j] <= tol) { numZero++; zeroAt = j; }
      if (bestBary[j] >= 1.0 - tol) oneAt = j;
    }
    if (oneAt >= 0)
    {
      onOpenEdge = d.boundaryPoint[(size_t)tt[oneAt]] != 0;
    }
    else if (numZero == 1)
    {
      int u = (zeroAt + 1)%3, v = (zeroAt + 2)%3;
      onOpenEdge = d.boundaryPoint[(size_t)tt[u]] != 0 && d.boundaryPoint[(size_t)tt[v]] != 0;
    }
  }
  double off[3];
  Sub(x, bestQ, off);
  if (!onOpenEdge && Dot(off, n) < 0.0)
  {
    // In the lumen: inside, by the wall there and the depth.
    double tq = bestBary[0]*d.thickness[(size_t)tt[0]] + bestBary[1]*d.thickness[(size_t)tt[1]] + bestBary[2]*d.thickness[(size_t)tt[2]];
    return -best - thicknessScale*tq;
  }
  return bestTerm;
}

void OffsetField::Local(const double x[3], double &size, double &thickness, long long &rim) const
{
  const Data &d = *data_;
  const TriangleGrid &g = d.grid;
  int c[3];
  for (int k = 0; k < 3; k++)
  {
    c[k] = g.Bin(x[k], k);
  }
  // Two triangles are found: the nearest by distance, which gives the size
  // and the wall the decimation works to (as it always has: taking them
  // from the other triangle let four collapses cross at a junction of the
  // user's model), and the one the field takes its value from at x, the
  // least distance less the wall there as Evaluate has it, which is the
  // piece of the interface x's offset stands on and so decides the cap. The
  // nearest by distance can be a thinner neighbour's, and a point would then
  // be handed to that neighbour's cap.
  double best = std::numeric_limits<double>::max(), bestTerm = std::numeric_limits<double>::max();
  ll bestT = -1, ownerT = -1;
  // Rings out until one lies farther than either could be beaten from; the
  // whole grid at most, since x may be far from the surface.
  int rings = std::max(g.n[0], std::max(g.n[1], g.n[2]));
  for (int r = 0; r <= rings; r++)
  {
    // every bin of this ring is at least (r-1) cells away
    double bound = (r - 1)*g.cell;
    if (bestT >= 0 && bound >= best && bound - d.largestThickness >= bestTerm)
    {
      break;
    }
    for (int k = c[2] - r; k <= c[2] + r; k++)
    {
      for (int j = c[1] - r; j <= c[1] + r; j++)
      {
        for (int i = c[0] - r; i <= c[0] + r; i++)
        {
          if (std::max(std::abs(k - c[2]), std::max(std::abs(j - c[1]), std::abs(i - c[0]))) != r) continue;
          if (i < 0 || j < 0 || k < 0 || i >= g.n[0] || j >= g.n[1] || k >= g.n[2]) continue;
          size_t b = g.Index(i, j, k);
          if (g.start[b + 1] == g.start[b]) continue;
          double dd = g.BoxDistance(x, i, j, k);
          if (dd >= best && dd - g.binMaxThickness[b] >= bestTerm) continue;
          for (ll m = g.start[b]; m < g.start[b + 1]; m++)
          {
            ll t = g.cells[(size_t)m];
            const ll *tt = &d.triangles[(size_t)3*t];
            double q[3], bary[3];
            double dist = ClosestOnTriangle(x, &d.points[(size_t)3*tt[0]], &d.points[(size_t)3*tt[1]], &d.points[(size_t)3*tt[2]], q, bary);
            if (dist < best)
            {
              best = dist;
              bestT = t;
            }
            double tq = bary[0]*d.thickness[(size_t)tt[0]] + bary[1]*d.thickness[(size_t)tt[1]] + bary[2]*d.thickness[(size_t)tt[2]];
            if (dist - tq < bestTerm)
            {
              bestTerm = dist - tq;
              ownerT = t;
            }
          }
        }
      }
    }
  }
  size = (bestT >= 0) ? d.localSize[(size_t)bestT] : d.meanEdge;
  thickness = (bestT >= 0) ? d.localThickness[(size_t)bestT] : d.largestThickness;
  rim = (ownerT >= 0) ? d.triangleRim[(size_t)ownerT] : -1;
}

double OffsetField::Reach() const
{
  return data_->reach;
}

const std::vector<std::vector<long long> > &OffsetField::Rims() const
{
  return data_->rims;
}

const std::vector<double> &OffsetField::Points() const { return data_->points; }
const std::vector<double> &OffsetField::Normals() const { return data_->normals; }
const std::vector<double> &OffsetField::Thickness() const { return data_->thickness; }
const std::vector<long long> &OffsetField::Triangles() const { return data_->triangles; }
long long OffsetField::NumEvaluations() const { return data_->numEvaluations; }

//---------------------
// BuildOffsetSurface
//---------------------

int BuildOffsetSurfaces(const Interface &input, const Options &options,
    const std::vector<double> &fractions, DelaunayFunction delaunay, void *context,
    std::vector<Surface> &surfaces, std::vector<Report> &reports,
    std::string &error, ProgressFunction progress)
{
  reports.assign(fractions.size(), Report());
  surfaces.assign(fractions.size(), Surface());
  if (fractions.empty())
  {
    error = "no fraction of the thickness was asked for";
    return 1;
  }
  for (size_t f = 0; f < fractions.size(); f++)
  {
    if (!(fractions[f] > 0.0 && fractions[f] <= 1.0))
    {
      error = "a fraction of the thickness must be in (0, 1]";
      return 1;
    }
  }
  // The field and the cloud are built once; what they report goes into
  // every level's report.
  Report base;
  Report &report = base;
  if (delaunay == nullptr)
  {
    error = "no tetrahedralization was supplied";
    return 1;
  }
  if (!(options.innerLayer > 0.0 && options.innerLayer < 1.0) || !(options.outerLayer > 1.0) ||
      !(options.farLayer > options.outerLayer) || !(options.farSpacing > 0.0) ||
      !(options.layerJitter >= 0.0 && options.layerJitter < 0.5))
  {
    error = "the cloud layers must be 0 < inner < 1 < outer < far, with a positive far spacing and a sideways offset under half";
    return 1;
  }
  auto say = [&](const char *stage)
  {
    if (progress != nullptr)
    {
      progress(stage, context);
    }
  };
  auto t0 = std::chrono::steady_clock::now();
  say("the distance field: the interface, its collars and their search grid");
  OffsetField field;
  if (field.Build(input, report, error, options.chordTolerance) != 0)
  {
    return 1;
  }

  // The cloud: every field-surface point, its inner, outer and far offsets
  // along its normal, and past each collar's end the outer and far offsets
  // along the outward direction, so that the dome the field closes around
  // the collar's end is bounded and contoured (it comes off with the trim).
  // The value at a surface point is -t, and at the inner layer at most
  // -(1 - inner)t: both inside, so only the outer layers are evaluated.
  say("the point cloud and the field on its outer layer");
  const std::vector<double> &fp = field.Points();
  const std::vector<double> &fnrm = field.Normals();
  const std::vector<double> &ft = field.Thickness();
  const std::vector<ll> &ftris = field.Triangles();
  ll numFP = (ll)(fp.size()/3);
  // The edge length at each point, for the far layer: the mean of its edges.
  std::vector<double> edgeAt((size_t)numFP, 0.0);
  std::vector<int> edgeCount((size_t)numFP, 0);
  for (size_t m = 0; m + 2 < ftris.size(); m += 3)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = ftris[m + j], b = ftris[m + (j+1)%3];
      double L = Distance(&fp[(size_t)3*a], &fp[(size_t)3*b]);
      edgeAt[(size_t)a] += L;
      edgeAt[(size_t)b] += L;
      edgeCount[(size_t)a]++;
      edgeCount[(size_t)b]++;
    }
  }
  for (ll i = 0; i < numFP; i++)
  {
    edgeAt[(size_t)i] = (edgeCount[(size_t)i] > 0) ? edgeAt[(size_t)i]/edgeCount[(size_t)i] : 0.0;
  }
  // The far distance at a point: the far layer of the thickness, and no less
  // than farSpacing edges (see Options).
  auto farDistance = [&](ll i)
  {
    return std::max(options.farLayer*ft[(size_t)i], options.farSpacing*edgeAt[(size_t)i]);
  };
  // The value of a cloud point at another fraction of the thickness is
  // recomputed per level below; a surface point's is -fraction*t without
  // evaluation, every other point's is evaluated.
  std::vector<double> cloud, value;
  std::vector<signed char> cloudIsSurface;
  std::vector<ll> cloudBase;
  auto addPoint = [&](const double p[3], double v, bool onSurface, ll basePoint)
  {
    cloud.insert(cloud.end(), p, p + 3);
    value.push_back(v);
    cloudIsSurface.push_back(onSurface ? 1 : 0);
    cloudBase.push_back(basePoint);
  };
  // A layer point leaves its ray sideways by the option's fraction of its
  // distance out, in one of three directions 120 degrees apart around the
  // ray, the first chosen by a hash of the ray's point: the four points of
  // a ray are then the corners of a proper tetrahedron and not a line (see
  // Options::layerJitter for what a line did to the Delaunay).
  auto jitter = [&](double x[3], const double n[3], ll i, int third, double distance)
  {
    if (!(options.layerJitter > 0.0)) return;
    int axis = (std::fabs(n[0]) <= std::fabs(n[1]) && std::fabs(n[0]) <= std::fabs(n[2])) ? 0 : ((std::fabs(n[1]) <= std::fabs(n[2])) ? 1 : 2);
    double e[3] = {0.0, 0.0, 0.0};
    e[axis] = 1.0;
    double u[3], w[3];
    Cross(n, e, u);
    if (!Normalize(u)) return;
    Cross(n, u, w);
    unsigned long long h = (unsigned long long)(i + 1)*11400714819323198485ULL;
    h ^= h >> 29;
    double theta = 6.283185307179586*((double)(h & 0xFFFFFFFFULL)/4294967296.0 + third/3.0);
    double d = options.layerJitter*distance;
    for (int k = 0; k < 3; k++)
    {
      x[k] += d*(std::cos(theta)*u[k] + std::sin(theta)*w[k]);
    }
  };
  for (ll i = 0; i < numFP; i++)
  {
    const double *p = &fp[(size_t)3*i], *n = &fnrm[(size_t)3*i];
    double t = ft[(size_t)i];
    addPoint(p, -t, true, i);
    double a[3], b[3], c[3];
    double far = farDistance(i);
    for (int k = 0; k < 3; k++)
    {
      a[k] = p[k] + options.innerLayer*t*n[k];
      b[k] = p[k] + options.outerLayer*t*n[k];
      c[k] = p[k] + far*n[k];
    }
    jitter(a, n, i, 0, options.innerLayer*t);
    jitter(b, n, i, 1, options.outerLayer*t);
    jitter(c, n, i, 2, far);
    // The inner point's value is put to the field like the others': on the
    // ray it would be -(1 - innerLayer) t exactly, off it not quite.
    addPoint(a, field.Evaluate(a), false, i);
    addPoint(b, field.Evaluate(b), false, i);
    addPoint(c, field.Evaluate(c), false, i);
  }
  {
    const std::vector<std::vector<ll> > &rims = field.Rims();
    ll numInterfacePts = report.numInterfacePoints;
    ll base = numInterfacePts;
    for (size_t r = 0; r < rims.size(); r++)
    {
      const std::vector<ll> &loop = rims[r];
      double outward[3];
      // the collar points follow the interface points in rim order
      for (size_t m = 0; m < loop.size(); m++)
      {
        ll v = base + (ll)m;
        const double *p = &fp[(size_t)3*v];
        // the outward direction: from the rim point to its collar point
        Sub(p, &fp[(size_t)3*loop[m]], outward);
        if (!Normalize(outward))
        {
          continue;
        }
        double b[3], c[3];
        double t = ft[(size_t)v];
        double far = farDistance(v);
        for (int k = 0; k < 3; k++)
        {
          b[k] = p[k] + options.outerLayer*t*outward[k];
          c[k] = p[k] + far*outward[k];
        }
        jitter(b, outward, v, 1, options.outerLayer*t);
        jitter(c, outward, v, 2, far);
        addPoint(b, field.Evaluate(b), false, v);
        addPoint(c, field.Evaluate(c), false, v);
      }
      base += (ll)loop.size();
    }
  }
  report.numCloudPoints = (ll)(cloud.size()/3);
  auto t1 = std::chrono::steady_clock::now();
  report.secondsField = std::chrono::duration<double>(t1 - t0).count();

  for (size_t level = 0; level < fractions.size(); level++)
  {
  const double fraction = fractions[level];
  Report &report = reports[level];
  report = base;
  Surface &surface = surfaces[level];
  char stage[160];
  auto sayLevel = [&](const char *what)
  {
    if (fractions.size() > 1)
    {
      snprintf(stage, sizeof(stage), "%s (surface %zu of %zu, at %.3g of the thickness)", what, level + 1, fractions.size(), fraction);
      say(stage);
    }
    else
    {
      say(what);
    }
  };
  // The field on the cloud at this fraction of the thickness.
  std::vector<double> levelValue;
  if (fraction == 1.0)
  {
    levelValue = value;
  }
  else
  {
    sayLevel("the field on the cloud at this fraction of the thickness");
    levelValue.resize(value.size());
    for (size_t i = 0; i < value.size(); i++)
    {
      levelValue[i] = cloudIsSurface[i] ? -fraction*ft[(size_t)cloudBase[i]] : field.Evaluate(&cloud[3*i], fraction);
    }
  }
  {
  std::vector<double> &value = levelValue;
  // A cloud point can lie on the zero level itself, within rounding: a
  // branch's far layer landing on its parent's offset, say. Marching
  // tetrahedra keyed by edge would then put one contour point per edge on
  // that one spot, and the surface would carry triangles that touch without
  // sharing a point, which the checks read as crossings and the mesher as
  // folded facets. Such a value is made exactly zero here, and the marching
  // below gives the point one contour point of its own.
  ll numZero = 0;
  {
    double tiny = 1.0e-9*std::max(report.largestThickness, 0.0);
    for (size_t i = 0; i < value.size(); i++)
    {
      if (value[i] != 0.0 && std::fabs(value[i]) <= tiny)
      {
        value[i] = 0.0;
      }
      if (value[i] == 0.0)
      {
        numZero++;
      }
    }
    report.numZeroCloudPoints = numZero;
  }

  sayLevel("the Delaunay tetrahedralization of the cloud");
  auto t1 = std::chrono::steady_clock::now();
  std::vector<ll> tets;
  if (!delaunay(cloud, tets, context, error))
  {
    if (error.empty())
    {
      error = "the tetrahedralization of the cloud failed";
    }
    return 1;
  }
  if (tets.size() % 4 != 0)
  {
    error = "the tetrahedralization did not return four points per tetrahedron";
    return 1;
  }
  ll numCloud = report.numCloudPoints;
  for (size_t m = 0; m < tets.size(); m++)
  {
    if (tets[m] < 0 || tets[m] >= numCloud)
    {
      error = "a tetrahedron refers to a point outside the cloud";
      return 1;
    }
  }
  report.numTetrahedra = (ll)(tets.size()/4);
  auto t2 = std::chrono::steady_clock::now();
  report.secondsDelaunay = std::chrono::duration<double>(t2 - t1).count();

  sayLevel("marching tetrahedra over the cloud");
  // Marching tetrahedra: one contour point per tetrahedron edge whose ends
  // have values of opposite sign (a value of zero counts as outside), placed
  // by linear interpolation and then put on the zero level by secant steps
  // along the edge, bisection when the secant leaves the bracket, until the
  // field there is within the tolerance or the steps run out; a tetrahedron
  // with one corner inside gives one triangle, with two, a quad as two
  // triangles. The steps are counted out one by one rather than taken once:
  // where two walls' bands come close - a thin branch's root on a thick
  // parent - an edge crosses the ridge of the field between them and is
  // far from linear, and one step from the linear guess left points off
  // the level by up to one and a half thicknesses (measured 2026-09-22 on
  // the user's model: the decimation then folded triangles over one another
  // at two such roots, 12 crossings, gone with the steps run to convergence). A cloud point whose value
  // is zero is on the level already: every cut edge ending there gets the
  // one contour point placed on it (not one per edge), a triangle that then
  // names a point twice has no area and is left out, and a triangle emitted
  // twice (a face with all three corners on the level and the inside on both
  // sides of it) is a sheet of no thickness and is left out with its twin. The winding comes from the
  // tetrahedron itself, not from the triangle's geometry (which a sliver
  // cannot be trusted for): with the corners ordered so that the
  // tetrahedron is positively oriented, the face opposite corner i wound
  // outward is (1,2,3), (0,3,2), (0,1,3), (0,2,1) for i = 0..3, and the cut
  // around an inside corner is wound like the face opposite it; the quad
  // between an inside pair (i,j) and an outside pair (k,l) runs
  // e(i,k), e(i,l), e(j,l), e(j,k) when (i,j,k,l) is an even permutation
  // of (0,1,2,3). Both faces of the cut on a face shared by two tetrahedra
  // then agree, so every edge of the surface is traversed once each way.
  std::vector<double> &pts = surface.points;
  std::vector<ll> &tris = surface.triangles;
  std::unordered_map<unsigned long long, ll> onEdge;
  onEdge.reserve((size_t)numCloud*2);
  std::vector<ll> zeroPoint;
  if (numZero > 0)
  {
    zeroPoint.assign((size_t)numCloud, -1);
  }
  auto pointAt = [&](ll v) -> ll
  {
    if (zeroPoint[(size_t)v] < 0)
    {
      zeroPoint[(size_t)v] = (ll)(pts.size()/3);
      pts.insert(pts.end(), &cloud[(size_t)3*v], &cloud[(size_t)3*v] + 3);
    }
    return zeroPoint[(size_t)v];
  };
  auto edgePoint = [&](ll a, ll b) -> ll
  {
    if (a > b) std::swap(a, b);
    // a cut edge has one end inside, so at most one on the level
    if (value[(size_t)a] == 0.0)
    {
      return pointAt(a);
    }
    if (value[(size_t)b] == 0.0)
    {
      return pointAt(b);
    }
    unsigned long long key = (unsigned long long)a*4294967311ULL + (unsigned long long)b;
    std::unordered_map<unsigned long long, ll>::iterator it = onEdge.find(key);
    if (it != onEdge.end())
    {
      return it->second;
    }
    double va = value[(size_t)a], vb = value[(size_t)b];
    const double *pa = &cloud[(size_t)3*a], *pb = &cloud[(size_t)3*b];
    double s = va/(va - vb);
    double lo = 0.0, hi = 1.0, flo = va, fhi = vb;
    const double closeEnough = options.snapTolerance*std::max(std::fabs(va), std::fabs(vb));
    bool onLevel = options.snapIterations <= 0;
    for (int it2 = 0; it2 < options.snapIterations; it2++)
    {
      double pm[3];
      for (int k = 0; k < 3; k++)
      {
        pm[k] = pa[k] + s*(pb[k] - pa[k]);
      }
      double fm = field.Evaluate(pm, fraction);
      if (std::fabs(fm) <= closeEnough)
      {
        onLevel = true;
        break;
      }
      if ((fm < 0.0) == (flo < 0.0))
      {
        lo = s;
        flo = fm;
      }
      else
      {
        hi = s;
        fhi = fm;
      }
      double next = (flo*hi - fhi*lo)/(flo - fhi);
      s = (next > lo && next < hi) ? next : 0.5*(lo + hi);
    }
    if (!onLevel)
    {
      report.numContourPointsOffLevel++;
    }
    ll id = (ll)(pts.size()/3);
    for (int k = 0; k < 3; k++)
    {
      pts.push_back(pa[k] + s*(pb[k] - pa[k]));
    }
    onEdge[key] = id;
    return id;
  };
  auto emit = [&](ll a, ll b, ll d)
  {
    if (a == b || b == d || d == a)
    {
      report.numDegenerateContourTriangles++;
      return;
    }
    tris.push_back(a); tris.push_back(b); tris.push_back(d);
  };
  const int oppositeFace[4][3] = {{1, 2, 3}, {0, 3, 2}, {0, 1, 3}, {0, 2, 1}};
  for (size_t q = 0; q + 3 < tets.size(); q += 4)
  {
    ll c[4] = {tets[q], tets[q + 1], tets[q + 2], tets[q + 3]};
    bool in[4];
    int numIn = 0;
    for (int m = 0; m < 4; m++)
    {
      in[m] = value[(size_t)c[m]] < 0.0;
      if (in[m]) numIn++;
    }
    if (numIn == 0 || numIn == 4)
    {
      continue;
    }
    report.numTetrahedraCut++;
    // positively oriented
    {
      double e1[3], e2[3], e3[3], cr[3];
      Sub(&cloud[(size_t)3*c[1]], &cloud[(size_t)3*c[0]], e1);
      Sub(&cloud[(size_t)3*c[2]], &cloud[(size_t)3*c[0]], e2);
      Sub(&cloud[(size_t)3*c[3]], &cloud[(size_t)3*c[0]], e3);
      Cross(e2, e3, cr);
      if (Dot(e1, cr) < 0.0)
      {
        std::swap(c[2], c[3]);
        std::swap(in[2], in[3]);
      }
    }
    if (numIn == 1 || numIn == 3)
    {
      int apex = 0;
      for (int m = 0; m < 4; m++)
      {
        if (in[m] == (numIn == 1)) apex = m;
      }
      const int *f = oppositeFace[apex];
      ll a = edgePoint(c[apex], c[f[0]]), b = edgePoint(c[apex], c[f[1]]), d = edgePoint(c[apex], c[f[2]]);
      if (numIn == 1)
      {
        emit(a, b, d);
      }
      else
      {
        emit(a, d, b);
      }
    }
    else
    {
      int i = -1, j = -1, k = -1, l = -1;
      for (int m = 0; m < 4; m++)
      {
        if (in[m]) { if (i < 0) i = m; else j = m; }
        else { if (k < 0) k = m; else l = m; }
      }
      // the parity of (i, j, k, l) as a permutation of (0, 1, 2, 3)
      int perm[4] = {i, j, k, l};
      int inversions = 0;
      for (int x = 0; x < 4; x++) for (int y = x + 1; y < 4; y++) if (perm[x] > perm[y]) inversions++;
      if (inversions % 2 == 1)
      {
        std::swap(k, l);
      }
      ll p00 = edgePoint(c[i], c[k]), p01 = edgePoint(c[i], c[l]), p11 = edgePoint(c[j], c[l]), p10 = edgePoint(c[j], c[k]);
      emit(p00, p01, p11);
      emit(p00, p11, p10);
    }
  }
  if (numZero > 0)
  {
    // A triangle emitted twice: only possible around points on the level,
    // as the cut of a face lies in the face and is shared by two tetrahedra
    // exactly. Both copies go, and their sheet of no thickness with them.
    std::unordered_map<unsigned long long, ll> seen;
    std::vector<unsigned char> twin(tris.size()/3, 0);
    for (size_t m = 0; m + 2 < tris.size(); m += 3)
    {
      ll u[3] = {tris[m], tris[m + 1], tris[m + 2]};
      std::sort(u, u + 3);
      unsigned long long key = ((unsigned long long)u[0]*4294967311ULL + (unsigned long long)u[1])*4294967291ULL + (unsigned long long)u[2];
      std::unordered_map<unsigned long long, ll>::iterator it = seen.find(key);
      if (it == seen.end())
      {
        seen[key] = (ll)(m/3);
      }
      else
      {
        ll other = it->second;
        const ll *w = &tris[(size_t)3*other];
        ll v[3] = {w[0], w[1], w[2]};
        std::sort(v, v + 3);
        if (v[0] == u[0] && v[1] == u[1] && v[2] == u[2])
        {
          twin[m/3] = 1;
          twin[(size_t)other] = 1;
        }
      }
    }
    size_t kept = 0;
    for (size_t m = 0; m < twin.size(); m++)
    {
      if (twin[m])
      {
        report.numDegenerateContourTriangles++;
        continue;
      }
      tris[3*kept] = tris[3*m];
      tris[3*kept + 1] = tris[3*m + 1];
      tris[3*kept + 2] = tris[3*m + 2];
      kept++;
    }
    tris.resize(3*kept);
  }
  report.numContourPoints = (ll)(pts.size()/3);
  report.numContourTriangles = (ll)(tris.size()/3);
  // The cloud and its tetrahedra are not needed past here, and the surface
  // that follows is large too.
  std::vector<ll>().swap(tets);
  std::vector<double>().swap(levelValue);
  std::unordered_map<unsigned long long, ll>().swap(onEdge);
  auto t3 = std::chrono::steady_clock::now();
  report.secondsContour = std::chrono::duration<double>(t3 - t2).count();
  sayLevel("the decimation of the contour to the interface's size");

  // Decimation: the contour has a point wherever the zero level crosses a
  // tetrahedron edge, several per interface point. Edges shorter than a
  // fraction of the local size are collapsed onto one of their ends, shortest
  // first, so every point that remains stays where it was, on the zero level,
  // and the surface stays a chord of the true offset. A collapse must keep
  // the surface a manifold (the link condition), make no edge longer than
  // the local size allows, and turn no triangle by more than the option's
  // cosine.
  if (options.collapseRatio > 0.0)
  {
    ll nv = (ll)(pts.size()/3), nt = (ll)(tris.size()/3);
    std::vector<double> target((size_t)nv, field.Reach()), thick((size_t)nv, fraction*report.largestThickness);
    for (ll v = 0; v < nv; v++)
    {
      ll rim;
      field.Local(&pts[(size_t)3*v], target[(size_t)v], thick[(size_t)v], rim);
      thick[(size_t)v] *= fraction;
    }
    std::vector<std::vector<ll> > incident((size_t)nv);
    std::vector<unsigned char> dead((size_t)nt, 0);
    for (ll t = 0; t < nt; t++)
    {
      for (int j = 0; j < 3; j++)
      {
        incident[(size_t)tris[(size_t)3*t + j]].push_back(t);
      }
    }
    std::vector<ll> version((size_t)nv, 0);
    struct Entry
    {
      double length;
      ll u, v, vu, vv;
      bool operator<(const Entry &o) const { return length > o.length; }
    };
    std::priority_queue<Entry> queue;
    const double growFactor = 1.25;
    auto pushEdges = [&](ll v)
    {
      for (size_t m = 0; m < incident[(size_t)v].size(); m++)
      {
        ll t = incident[(size_t)v][m];
        if (dead[(size_t)t]) continue;
        for (int j = 0; j < 3; j++)
        {
          ll a = tris[(size_t)3*t + j], b = tris[(size_t)3*t + (j+1)%3];
          if (a != v && b != v) continue;
          ll u = (a == v) ? b : a;
          double L = Distance(&pts[(size_t)3*u], &pts[(size_t)3*v]);
          if (L < options.collapseRatio*std::min(target[(size_t)u], target[(size_t)v]))
          {
            queue.push(Entry{L, u, v, version[(size_t)u], version[(size_t)v]});
          }
        }
      }
    };
    for (ll v = 0; v < nv; v++)
    {
      pushEdges(v);
    }
    auto normalOf = [&](ll a, ll b, ll c, double n[3])
    {
      double e1[3], e2[3];
      Sub(&pts[(size_t)3*b], &pts[(size_t)3*a], e1);
      Sub(&pts[(size_t)3*c], &pts[(size_t)3*a], e2);
      Cross(e1, e2, n);
    };
    // The crossing guard: every operation's new triangles are tested against
    // the live triangles near them, and one that passes through another is
    // not made. A uniform grid over the contour holds each live triangle in
    // the cells its box covers; a collapse re-indexes triangles (points never
    // move), so a triangle leaves its cells before and enters them after.
    struct GuardGrid
    {
      double origin[3];
      double cell;
      int n[3];
      std::unordered_map<ll, std::vector<ll> > cells;
      ll Index(int i, int j, int k) const { return ((ll)k*n[1] + j)*n[0] + i; }
      int Bin(double x, int k) const { int b = (int)std::floor((x - origin[k])/cell); return std::min(std::max(b, 0), n[k] - 1); }
    };
    GuardGrid guard;
    const bool guardOn = options.guardCrossings;
    auto triangleBox = [&](const ll *T, int lo[3], int hi[3])
    {
      for (int k = 0; k < 3; k++)
      {
        double a = pts[(size_t)3*T[0] + k], b = pts[(size_t)3*T[1] + k], c = pts[(size_t)3*T[2] + k];
        lo[k] = guard.Bin(std::min(a, std::min(b, c)), k);
        hi[k] = guard.Bin(std::max(a, std::max(b, c)), k);
      }
    };
    auto guardAdd = [&](ll t)
    {
      int lo[3], hi[3];
      triangleBox(&tris[(size_t)3*t], lo, hi);
      for (int k = lo[2]; k <= hi[2]; k++) for (int j = lo[1]; j <= hi[1]; j++) for (int i = lo[0]; i <= hi[0]; i++) guard.cells[guard.Index(i, j, k)].push_back(t);
    };
    auto guardRemove = [&](ll t)
    {
      int lo[3], hi[3];
      triangleBox(&tris[(size_t)3*t], lo, hi);
      for (int k = lo[2]; k <= hi[2]; k++) for (int j = lo[1]; j <= hi[1]; j++) for (int i = lo[0]; i <= hi[0]; i++)
      {
        std::unordered_map<ll, std::vector<ll> >::iterator it = guard.cells.find(guard.Index(i, j, k));
        if (it == guard.cells.end()) continue;
        std::vector<ll> &list = it->second;
        for (size_t m = 0; m < list.size(); m++) { if (list[m] == t) { list[m] = list.back(); list.pop_back(); break; } }
      }
    };
    std::vector<ll> guardStamp;
    ll guardStampValue = 0;
    if (guardOn)
    {
      double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
      for (ll v = 0; v < nv; v++) for (int k = 0; k < 3; k++) { lo[k] = std::min(lo[k], pts[(size_t)3*v + k]); hi[k] = std::max(hi[k], pts[(size_t)3*v + k]); }
      double meanTarget = 0.0;
      for (ll v = 0; v < nv; v++) meanTarget += target[(size_t)v];
      meanTarget = (nv > 0) ? meanTarget/(double)nv : 1.0;
      guard.cell = std::max(2.0*meanTarget, 1.0e-12);
      for (;;)
      {
        double count = 1.0;
        for (int k = 0; k < 3; k++) { guard.n[k] = std::max(1, (int)std::ceil((hi[k] - lo[k])/guard.cell) + 1); count *= guard.n[k]; }
        if (count < 4.0e8) break;
        guard.cell *= 2.0;
      }
      for (int k = 0; k < 3; k++) guard.origin[k] = lo[k];
      guard.cells.reserve((size_t)nt/4 + 1);
      for (ll t = 0; t < nt; t++) guardAdd(t);
      guardStamp.assign((size_t)nt, 0);
    }
    // Whether the triangle T (its corners' current positions) passes through
    // a live triangle near it other than those stamped as leaving.
    auto guardCrosses = [&](const ll *T) -> bool
    {
      int lo[3], hi[3];
      triangleBox(T, lo, hi);
      for (int k = lo[2]; k <= hi[2]; k++) for (int j = lo[1]; j <= hi[1]; j++) for (int i = lo[0]; i <= hi[0]; i++)
      {
        std::unordered_map<ll, std::vector<ll> >::iterator it = guard.cells.find(guard.Index(i, j, k));
        if (it == guard.cells.end()) continue;
        const std::vector<ll> &list = it->second;
        for (size_t m = 0; m < list.size(); m++)
        {
          ll c = list[m];
          if (dead[(size_t)c] || guardStamp[(size_t)c] == guardStampValue) continue;
          guardStamp[(size_t)c] = guardStampValue;
          if (svenvelope::TrianglesCross(pts, T, &tris[(size_t)3*c])) return true;
        }
      }
      return false;
    };
    const int ringMax = 64;
    ll ringU[64], ringV[64], shared[2];
    auto tryCollapse = [&](ll u, ll v) -> bool
    {
      // u goes onto v
      double moved = Distance(&pts[(size_t)3*u], &pts[(size_t)3*v]);
      int nu = 0, nvv = 0, ns = 0;
      for (size_t m = 0; m < incident[(size_t)u].size(); m++)
      {
        ll t = incident[(size_t)u][m];
        if (dead[(size_t)t]) continue;
        bool hasV = false;
        for (int j = 0; j < 3; j++)
        {
          ll w = tris[(size_t)3*t + j];
          if (w == v) hasV = true;
          else if (w != u)
          {
            bool seen = false;
            for (int q = 0; q < nu; q++) if (ringU[q] == w) seen = true;
            if (!seen) { if (nu >= ringMax) return false; ringU[nu++] = w; }
          }
        }
        if (hasV)
        {
          if (ns < 2) shared[ns] = t;
          ns++;
        }
      }
      if (ns != 2) return false;
      for (size_t m = 0; m < incident[(size_t)v].size(); m++)
      {
        ll t = incident[(size_t)v][m];
        if (dead[(size_t)t]) continue;
        for (int j = 0; j < 3; j++)
        {
          ll w = tris[(size_t)3*t + j];
          if (w != v && w != u)
          {
            bool seen = false;
            for (int q = 0; q < nvv; q++) if (ringV[q] == w) seen = true;
            if (!seen) { if (nvv >= ringMax) return false; ringV[nvv++] = w; }
          }
        }
      }
      int common = 0;
      for (int q = 0; q < nu; q++) for (int r = 0; r < nvv; r++) if (ringU[q] == ringV[r]) common++;
      if (common != 2) return false;
      // The edges that appear, from v to every neighbour w of u, may be no
      // longer than the local size allows; one that is already longer than
      // that (u-w, which v-w replaces) may grow by a tenth at most, so that a
      // flat cap whose long edges are what they are can still lose its short
      // one without the long edges creeping.
      for (int q = 0; q < nu; q++)
      {
        ll w = ringU[q];
        double L = Distance(&pts[(size_t)3*w], &pts[(size_t)3*v]);
        double allowed = growFactor*std::min(target[(size_t)w], target[(size_t)v]);
        if (L > allowed && L > 1.1*Distance(&pts[(size_t)3*w], &pts[(size_t)3*u])) return false;
      }
      for (size_t m = 0; m < incident[(size_t)u].size(); m++)
      {
        ll t = incident[(size_t)u][m];
        if (dead[(size_t)t] || t == shared[0] || t == shared[1]) continue;
        ll a = tris[(size_t)3*t], b = tris[(size_t)3*t + 1], c = tris[(size_t)3*t + 2];
        ll a2 = (a == u) ? v : a, b2 = (b == u) ? v : b, c2 = (c == u) ? v : c;
        double nOld[3], nNew[3];
        normalOf(a, b, c, nOld);
        normalOf(a2, b2, c2, nNew);
        double lo = Norm(nOld), ln = Norm(nNew);
        if (!(ln > 1.0e-12*lo) || Dot(nOld, nNew) < options.collapseTurnCosine*lo*ln) return false;
        // The surface moves by no more than the collapsed edge, so a short
        // collapse cannot take a triangle off the zero level by more than
        // the tolerance and is not put to the field.
        if (options.collapseFieldTolerance > 0.0 && moved > options.collapseFieldTolerance*thick[(size_t)v])
        {
          double centre[3];
          for (int k = 0; k < 3; k++)
          {
            centre[k] = (pts[(size_t)3*a2 + k] + pts[(size_t)3*b2 + k] + pts[(size_t)3*c2 + k])/3.0;
          }
          double thinnest = std::min(thick[(size_t)a2], std::min(thick[(size_t)b2], thick[(size_t)c2]));
          if (std::abs(field.Evaluate(centre, fraction)) > options.collapseFieldTolerance*thinnest) return false;
        }
      }
      if (guardOn)
      {
        // the triangles on u leave (the two on the edge die, the rest move)
        guardStampValue++;
        for (size_t m = 0; m < incident[(size_t)u].size(); m++) guardStamp[(size_t)incident[(size_t)u][m]] = guardStampValue;
        ll moved[ringMax][3];
        int numMoved = 0;
        for (size_t m = 0; m < incident[(size_t)u].size(); m++)
        {
          ll t = incident[(size_t)u][m];
          if (dead[(size_t)t] || t == shared[0] || t == shared[1]) continue;
          if (numMoved >= ringMax) return false;
          for (int j = 0; j < 3; j++) moved[numMoved][j] = (tris[(size_t)3*t + j] == u) ? v : tris[(size_t)3*t + j];
          numMoved++;
        }
        for (int i = 0; i < numMoved; i++)
        {
          if (guardCrosses(moved[i])) { report.numRefusedForCrossing++; return false; }
          for (int j = i + 1; j < numMoved; j++) if (svenvelope::TrianglesCross(pts, moved[i], moved[j])) { report.numRefusedForCrossing++; return false; }
        }
      }
      dead[(size_t)shared[0]] = 1;
      dead[(size_t)shared[1]] = 1;
      if (guardOn) { guardRemove(shared[0]); guardRemove(shared[1]); }
      for (size_t m = 0; m < incident[(size_t)u].size(); m++)
      {
        ll t = incident[(size_t)u][m];
        if (dead[(size_t)t]) continue;
        if (guardOn) guardRemove(t);
        for (int j = 0; j < 3; j++)
        {
          if (tris[(size_t)3*t + j] == u) tris[(size_t)3*t + j] = v;
        }
        if (guardOn) guardAdd(t);
        incident[(size_t)v].push_back(t);
      }
      incident[(size_t)u].clear();
      version[(size_t)u]++;
      version[(size_t)v]++;
      report.numCollapsed++;
      return true;
    };
    auto runCollapses = [&]()
    {
      while (!queue.empty())
      {
        Entry e = queue.top();
        queue.pop();
        if (e.vu != version[(size_t)e.u] || e.vv != version[(size_t)e.v]) continue;
        if (incident[(size_t)e.u].empty() || incident[(size_t)e.v].empty()) continue;
        bool done = tryCollapse(e.u, e.v);
        if (!done) done = tryCollapse(e.v, e.u);
        if (done) pushEdges(incident[(size_t)e.u].empty() ? e.v : e.u);
      }
    };
    // The smallest angle of a triangle, as its sine.
    auto smallestAngleSine = [&](ll a, ll b, ll c)
    {
      double la = Distance(&pts[(size_t)3*b], &pts[(size_t)3*c]), lb = Distance(&pts[(size_t)3*c], &pts[(size_t)3*a]), lc = Distance(&pts[(size_t)3*a], &pts[(size_t)3*b]);
      double n[3];
      normalOf(a, b, c, n);
      double twiceArea = Norm(n);
      double longest = std::max(la, std::max(lb, lc));
      double product = la*lb*lc;
      if (!(product > 0.0)) return 0.0;
      // sin(smallest angle) = 2A / (product of the two longer sides) = 2A * shortest / (la lb lc)
      double shortest = std::min(la, std::min(lb, lc));
      (void)longest;
      return twiceArea*shortest/product;
    };
    auto removeFromIncident = [&](ll v, ll t)
    {
      std::vector<ll> &list = incident[(size_t)v];
      for (size_t m = 0; m < list.size(); m++)
      {
        if (list[m] == t)
        {
          list[m] = list.back();
          list.pop_back();
          return;
        }
      }
    };
    // Flips: the edge shared by two nearly coplanar triangles is flipped to
    // the other diagonal of their quad when that raises the smaller of their
    // smallest angles and the quad is convex (both new triangles face the way
    // the old ones did). Nothing moves, so nothing can leave the zero level.
    std::vector<ll> live;
    auto runFlips = [&]() -> ll
    {
      ll numFlipped = 0;
      for (size_t li = 0; li < live.size(); li++)
      {
        ll t1 = live[li];
        if (dead[(size_t)t1]) continue;
        for (int j = 0; j < 3; j++)
        {
          if (dead[(size_t)t1]) break;
          ll a = tris[(size_t)3*t1 + j], b = tris[(size_t)3*t1 + (j+1)%3], c = tris[(size_t)3*t1 + (j+2)%3];
          // the other triangle on edge a-b
          ll t2 = -1;
          for (size_t m = 0; m < incident[(size_t)a].size(); m++)
          {
            ll t = incident[(size_t)a][m];
            if (t == t1 || dead[(size_t)t]) continue;
            for (int q = 0; q < 3; q++)
            {
              if (tris[(size_t)3*t + q] == b && tris[(size_t)3*t + (q+1)%3] == a)
              {
                t2 = t;
              }
            }
          }
          if (t2 < 0) continue;
          ll d = -1;
          for (int q = 0; q < 3; q++)
          {
            ll w = tris[(size_t)3*t2 + q];
            if (w != a && w != b) d = w;
          }
          if (d < 0 || d == c) continue;
          // c-d must not be an edge already
          bool exists = false;
          for (size_t m = 0; m < incident[(size_t)c].size() && !exists; m++)
          {
            ll t = incident[(size_t)c][m];
            if (dead[(size_t)t]) continue;
            for (int q = 0; q < 3; q++) if (tris[(size_t)3*t + q] == d) exists = true;
          }
          if (exists) continue;
          double n1[3], n2[3], m1[3], m2[3];
          normalOf(a, b, c, n1);
          normalOf(b, a, d, n2);
          normalOf(c, a, d, m1);
          normalOf(d, b, c, m2);
          double l1 = Norm(n1), l2 = Norm(n2), k1 = Norm(m1), k2 = Norm(m2);
          if (!(l1 > 0.0 && l2 > 0.0 && k1 > 0.0 && k2 > 0.0)) continue;
          // The new pair must face the way the old pair does, taken together
          // with their areas (a sliver's own normal means little) and within
          // 25 degrees, which also keeps the flip off a crease and out of a
          // quad that is not convex.
          double avg[3] = {n1[0] + n2[0], n1[1] + n2[1], n1[2] + n2[2]};
          double la = Norm(avg);
          if (!(la > 0.0) || Dot(m1, avg) < 0.9*k1*la || Dot(m2, avg) < 0.9*k2*la) continue;
          // Two well-shaped triangles that meet at a crease stay as they
          // are; a sliver's normal is noise and does not get that say.
          auto aspectOf = [&](ll x, ll y, ll z, double twiceArea)
          {
            double lx = Distance(&pts[(size_t)3*y], &pts[(size_t)3*z]), ly = Distance(&pts[(size_t)3*z], &pts[(size_t)3*x]), lz = Distance(&pts[(size_t)3*x], &pts[(size_t)3*y]);
            double perimeter = lx + ly + lz;
            return (twiceArea > 0.0 && perimeter > 0.0) ? std::max(lx, std::max(ly, lz))*perimeter/(2.0*twiceArea) : 1.0e300;
          };
          if (aspectOf(a, b, c, l1) < 20.0 && aspectOf(b, a, d, l2) < 20.0 && Dot(n1, n2) < 0.9*l1*l2) continue;
          // and the new diagonal no longer than the old edge or the local size
          // allows - or, when one of the pair is a flat cap (its apex on the
          // shared edge), no longer than the other's own edges: the flip then
          // only splits that other triangle along a segment from its edge to
          // its apex, which is never longer than its sides (9 caps of the
          // user's model, apex valence 4, were kept by the size rule alone).
          double newDiagonal = Distance(&pts[(size_t)3*c], &pts[(size_t)3*d]);
          double allowedDiagonal = std::max(Distance(&pts[(size_t)3*a], &pts[(size_t)3*b]), growFactor*std::min(target[(size_t)c], target[(size_t)d]));
          if (aspectOf(a, b, c, l1) > 40.0 || aspectOf(b, a, d, l2) > 40.0)
          {
            double sides = std::max(std::max(Distance(&pts[(size_t)3*a], &pts[(size_t)3*d]), Distance(&pts[(size_t)3*b], &pts[(size_t)3*d])),
                std::max(Distance(&pts[(size_t)3*a], &pts[(size_t)3*c]), Distance(&pts[(size_t)3*b], &pts[(size_t)3*c])));
            allowedDiagonal = std::max(allowedDiagonal, sides);
          }
          if (newDiagonal > allowedDiagonal) continue;
          // Worth it when the smaller of the smallest angles grows - or when
          // the new diagonal is short enough to be collapsed in the next
          // round: two flat caps on one long edge form a thin quad that
          // neither a collapse nor a better-angle flip can take apart, but
          // its short diagonal, once an edge, can.
          double before = std::min(smallestAngleSine(a, b, c), smallestAngleSine(b, a, d));
          double after = std::min(smallestAngleSine(c, a, d), smallestAngleSine(d, b, c));
          bool shortDiagonal = newDiagonal < options.collapseRatio*std::min(target[(size_t)c], target[(size_t)d]) && before < 0.05;
          if (!(after > 1.05*before) && !shortDiagonal) continue;
          if (guardOn)
          {
            guardStampValue++;
            guardStamp[(size_t)t1] = guardStampValue;
            guardStamp[(size_t)t2] = guardStampValue;
            ll new1[3] = {c, a, d}, new2[3] = {d, b, c};
            if (guardCrosses(new1) || guardCrosses(new2)) { report.numRefusedForCrossing++; continue; }
            guardRemove(t1);
            guardRemove(t2);
          }
          // do it: t1 = (c, a, d), t2 = (d, b, c)
          tris[(size_t)3*t1] = c; tris[(size_t)3*t1 + 1] = a; tris[(size_t)3*t1 + 2] = d;
          tris[(size_t)3*t2] = d; tris[(size_t)3*t2 + 1] = b; tris[(size_t)3*t2 + 2] = c;
          if (guardOn) { guardAdd(t1); guardAdd(t2); }
          removeFromIncident(b, t1);
          incident[(size_t)d].push_back(t1);
          removeFromIncident(a, t2);
          incident[(size_t)c].push_back(t2);
          version[(size_t)a]++; version[(size_t)b]++; version[(size_t)c]++; version[(size_t)d]++;
          numFlipped++;
          break;   // t1 changed; its other edges are seen in the next pass
        }
      }
      return numFlipped;
    };
    // A point on three triangles whose edges are all short is taken out and
    // its three triangles made one, when that one faces the way they did.
    auto removeValenceThree = [&]() -> ll
    {
      ll numRemoved = 0;
      for (ll u = 0; u < nv; u++)
      {
        ll live[3];
        int n = 0;
        bool tooMany = false;
        for (size_t m = 0; m < incident[(size_t)u].size() && !tooMany; m++)
        {
          ll t = incident[(size_t)u][m];
          if (dead[(size_t)t]) continue;
          if (n < 3) live[n] = t;
          n++;
          if (n > 3) tooMany = true;
        }
        if (tooMany || n != 3) continue;
        // the ring a -> b -> c, from the triangles (u, a, b), (u, b, c), (u, c, a)
        ll ring[3];
        int nr = 0;
        ll t0 = live[0];
        for (int q = 0; q < 3; q++)
        {
          if (tris[(size_t)3*t0 + q] == u)
          {
            ring[0] = tris[(size_t)3*t0 + (q+1)%3];
            ring[1] = tris[(size_t)3*t0 + (q+2)%3];
            nr = 2;
          }
        }
        if (nr != 2) continue;
        for (int i = 1; i < 3 && nr == 2; i++)
        {
          ll t = live[i];
          for (int q = 0; q < 3; q++)
          {
            ll w = tris[(size_t)3*t + q];
            if (w != u && w != ring[0] && w != ring[1])
            {
              ring[2] = w;
              nr = 3;
            }
          }
        }
        if (nr != 3) continue;
        double nNew[3];
        normalOf(ring[0], ring[1], ring[2], nNew);
        double ln = Norm(nNew);
        if (!(ln > 0.0)) continue;
        // A point lying (all but) on one of its ring's edges - the apex of a
        // flat cap, which is what the collapses leave behind when they snap
        // a point onto the zero level in a flat spot - has three triangles
        // whose union is the ring triangle to within the cap's height:
        // taking it out is not sized by the edge lengths (nothing grows),
        // and the flat triangle's own normal, which is noise, gets no say.
        // Measured on the user's model: the 199 triangles left above aspect
        // 100 were all such caps, 190 of them in pairs on one long edge with
        // apexes of valence 3, which no flip can take apart (the flipped
        // pair is as flat) and which this removal refused for the noise
        // normal.
        bool onRingEdge = false;
        for (int i = 0; i < 3; i++)
        {
          const ll *tt = &tris[(size_t)3*live[i]];
          double nOld[3];
          normalOf(tt[0], tt[1], tt[2], nOld);
          double lo = Norm(nOld);
          double e1 = Distance(&pts[(size_t)3*tt[0]], &pts[(size_t)3*tt[1]]), e2 = Distance(&pts[(size_t)3*tt[1]], &pts[(size_t)3*tt[2]]);
          if (lo <= 0.05*e1*e2) onRingEdge = true;   // sine of the angle at the middle corner under 0.05: an aspect ratio past about 40
        }
        if (!onRingEdge)
        {
          bool shortEdges = true;
          for (int i = 0; i < 3; i++)
          {
            if (Distance(&pts[(size_t)3*u], &pts[(size_t)3*ring[i]]) >= options.collapseRatio*std::min(target[(size_t)u], target[(size_t)ring[i]])) shortEdges = false;
          }
          if (!shortEdges) continue;
        }
        bool agree = true;
        for (int i = 0; i < 3 && agree; i++)
        {
          const ll *tt = &tris[(size_t)3*live[i]];
          double nOld[3];
          normalOf(tt[0], tt[1], tt[2], nOld);
          double lo = Norm(nOld);
          double e1 = Distance(&pts[(size_t)3*tt[0]], &pts[(size_t)3*tt[1]]), e2 = Distance(&pts[(size_t)3*tt[1]], &pts[(size_t)3*tt[2]]);
          if (lo <= 0.05*e1*e2) continue;   // a flat triangle has no direction to agree with
          if (Dot(nOld, nNew) < options.collapseTurnCosine*lo*ln) agree = false;
        }
        if (!agree) continue;
        double centre[3];
        for (int k = 0; k < 3; k++)
        {
          centre[k] = (pts[(size_t)3*ring[0] + k] + pts[(size_t)3*ring[1] + k] + pts[(size_t)3*ring[2] + k])/3.0;
        }
        if (options.collapseFieldTolerance > 0.0 && std::abs(field.Evaluate(centre, fraction)) > options.collapseFieldTolerance*thick[(size_t)u]) continue;
        if (guardOn)
        {
          guardStampValue++;
          for (int i = 0; i < 3; i++) guardStamp[(size_t)live[i]] = guardStampValue;
          if (guardCrosses(ring)) { report.numRefusedForCrossing++; continue; }
          for (int i = 0; i < 3; i++) guardRemove(live[i]);
        }
        // live[0] becomes the ring triangle; the other two die
        tris[(size_t)3*t0] = ring[0]; tris[(size_t)3*t0 + 1] = ring[1]; tris[(size_t)3*t0 + 2] = ring[2];
        if (guardOn) guardAdd(t0);
        dead[(size_t)live[1]] = 1;
        dead[(size_t)live[2]] = 1;
        incident[(size_t)ring[2]].push_back(t0);
        incident[(size_t)u].clear();
        for (int i = 0; i < 3; i++) version[(size_t)ring[i]]++;
        version[(size_t)u]++;
        numRemoved++;
        report.numCollapsed++;
      }
      return numRemoved;
    };
    // A flat cap whose apex has more than three triangles: the apex is taken
    // out and its ring - a polygon with the cap's long edge as one side and
    // the apex lying on that side - is fanned from one end of that edge,
    // whichever gives the better smallest angle. The union is what it was to
    // within the cap's height; no new edge is longer than the ring's own
    // edges or the local size allow, none doubles an edge that exists, and
    // every new triangle faces the way the old ones did. Measured on the
    // user's model after the valence-three pass learned to take flat apexes:
    // the 101 triangles left above aspect 100 were all such caps, in pairs
    // on one long edge with both apexes of valence 4, which no flip can help
    // (a flat pair flipped is a flat pair) and no collapse took (the edges
    // it would make grew past the size rule).
    auto removeFlatApexes = [&]() -> ll
    {
      ll numRemoved = 0;
      ll slots[ringMax];
      for (ll t0 = 0; t0 < nt; t0++)
      {
        if (dead[(size_t)t0]) continue;
        // the apex: the corner whose angle is within three degrees of straight
        int apexCorner = -1;
        for (int q = 0; q < 3 && apexCorner < 0; q++)
        {
          ll c = tris[(size_t)3*t0 + q], a = tris[(size_t)3*t0 + (q+1)%3], b = tris[(size_t)3*t0 + (q+2)%3];
          double ca[3], cb[3];
          Sub(&pts[(size_t)3*a], &pts[(size_t)3*c], ca);
          Sub(&pts[(size_t)3*b], &pts[(size_t)3*c], cb);
          double la = Norm(ca), lb = Norm(cb);
          if (!(la > 0.0 && lb > 0.0)) continue;
          if (Dot(ca, cb) < -0.9986*la*lb) apexCorner = q;
        }
        if (apexCorner < 0) continue;
        ll c = tris[(size_t)3*t0 + apexCorner], a = tris[(size_t)3*t0 + (apexCorner+1)%3], b = tris[(size_t)3*t0 + (apexCorner+2)%3];
        // the ring around c, walked from (c, a, b) through the triangles (c, r_i, r_i+1)
        int k = 0;
        ringU[k++] = a;
        ringU[k++] = b;
        slots[0] = t0;
        int numSlots = 1;
        ll prevT = t0;
        bool closed = false, bad = false;
        while (!bad && !closed)
        {
          ll last = ringU[k-1];
          ll found = -1, w = -1;
          for (size_t m = 0; m < incident[(size_t)c].size(); m++)
          {
            ll t = incident[(size_t)c][m];
            if (t == prevT || dead[(size_t)t]) continue;
            for (int q = 0; q < 3; q++)
            {
              if (tris[(size_t)3*t + q] == c && tris[(size_t)3*t + (q+1)%3] == last)
              {
                found = t;
                w = tris[(size_t)3*t + (q+2)%3];
              }
            }
          }
          if (found < 0 || numSlots >= ringMax) { bad = true; break; }
          slots[numSlots++] = found;
          prevT = found;
          if (w == a) { closed = true; break; }
          if (k >= ringMax) { bad = true; break; }
          ringU[k++] = w;
        }
        if (bad || !closed || k < 4 || numSlots != k) continue;
        // every live triangle on c must be in the ring
        int numLive = 0;
        for (size_t m = 0; m < incident[(size_t)c].size(); m++)
        {
          if (!dead[(size_t)incident[(size_t)c][m]]) numLive++;
        }
        if (numLive != k) continue;
        // the way the old triangles face, with their areas
        double avg[3] = {0.0, 0.0, 0.0};
        double longestRingEdge = 0.0;
        for (int i = 0; i < k; i++)
        {
          double n[3];
          normalOf(c, ringU[i], ringU[(i+1)%k], n);
          for (int d = 0; d < 3; d++) avg[d] += n[d];
          longestRingEdge = std::max(longestRingEdge, Distance(&pts[(size_t)3*ringU[i]], &pts[(size_t)3*ringU[(i+1)%k]]));
        }
        double la = Norm(avg);
        if (!(la > 0.0)) continue;
        // the two fans: from ringU[0] over (r_i, r_i+1), i = 1..k-2, and from
        // ringU[1] over (r_i, r_i+1), i = 2..k-1 with r_k = r_0
        int bestPivot = -1;
        double bestScore = 0.0;
        for (int pivot = 0; pivot < 2 && k - 2 <= ringMax; pivot++)
        {
          ll p = ringU[pivot];
          bool ok = true;
          double score = 1.0e300;
          for (int i = 0; i < k - 2 && ok; i++)
          {
            ll r1 = ringU[(pivot + 1 + i)%k], r2 = ringU[(pivot + 2 + i)%k];
            double n[3];
            normalOf(p, r1, r2, n);
            double ln = Norm(n);
            if (!(ln > 0.0) || Dot(n, avg) < options.collapseTurnCosine*ln*la) { ok = false; break; }
            score = std::min(score, smallestAngleSine(p, r1, r2));
            // the new edge p-r2 (r1 is p's ring neighbour on the first step): not already an edge, not too long
            if (i < k - 3)
            {
              bool exists = false;
              for (size_t m = 0; m < incident[(size_t)p].size() && !exists; m++)
              {
                ll t = incident[(size_t)p][m];
                if (dead[(size_t)t]) continue;
                for (int q = 0; q < 3; q++) if (tris[(size_t)3*t + q] == r2) exists = true;
              }
              if (exists) { ok = false; break; }
              double L = Distance(&pts[(size_t)3*p], &pts[(size_t)3*r2]);
              if (L > std::max(growFactor*std::min(target[(size_t)p], target[(size_t)r2]), 1.25*longestRingEdge)) { ok = false; break; }
            }
            if (options.collapseFieldTolerance > 0.0)
            {
              double centre[3];
              for (int d = 0; d < 3; d++)
              {
                centre[d] = (pts[(size_t)3*p + d] + pts[(size_t)3*r1 + d] + pts[(size_t)3*r2 + d])/3.0;
              }
              if (std::abs(field.Evaluate(centre, fraction)) > options.collapseFieldTolerance*thick[(size_t)c]) { ok = false; break; }
            }
          }
          if (ok && score > bestScore)
          {
            bestScore = score;
            bestPivot = pivot;
          }
        }
        if (bestPivot < 0) continue;
        if (guardOn)
        {
          guardStampValue++;
          for (int i = 0; i < k; i++) guardStamp[(size_t)slots[i]] = guardStampValue;
          ll pv = ringU[bestPivot];
          bool crosses = false;
          for (int i = 0; i < k - 2 && !crosses; i++)
          {
            ll fan[3] = {pv, ringU[(bestPivot + 1 + i)%k], ringU[(bestPivot + 2 + i)%k]};
            if (guardCrosses(fan)) crosses = true;
            for (int j = i + 1; j < k - 2 && !crosses; j++)
            {
              ll other[3] = {pv, ringU[(bestPivot + 1 + j)%k], ringU[(bestPivot + 2 + j)%k]};
              if (svenvelope::TrianglesCross(pts, fan, other)) crosses = true;
            }
          }
          if (crosses) { report.numRefusedForCrossing++; continue; }
          for (int i = 0; i < k; i++) guardRemove(slots[i]);
        }
        // apply: the old triangles leave their corners' lists, the fan takes k-2 of their slots
        for (int i = 0; i < k; i++)
        {
          ll t = slots[i];
          for (int q = 0; q < 3; q++) removeFromIncident(tris[(size_t)3*t + q], t);
        }
        ll p = ringU[bestPivot];
        for (int i = 0; i < k - 2; i++)
        {
          ll t = slots[i];
          ll r1 = ringU[(bestPivot + 1 + i)%k], r2 = ringU[(bestPivot + 2 + i)%k];
          tris[(size_t)3*t] = p; tris[(size_t)3*t + 1] = r1; tris[(size_t)3*t + 2] = r2;
          if (guardOn) guardAdd(t);
          incident[(size_t)p].push_back(t);
          incident[(size_t)r1].push_back(t);
          incident[(size_t)r2].push_back(t);
        }
        dead[(size_t)slots[k-2]] = 1;
        dead[(size_t)slots[k-1]] = 1;
        incident[(size_t)c].clear();
        for (int i = 0; i < k; i++) version[(size_t)ringU[i]]++;
        version[(size_t)c]++;
        numRemoved++;
        report.numCollapsed++;
      }
      return numRemoved;
    };
    for (int round = 0; round < 3; round++)
    {
      runCollapses();
      live.clear();
      for (ll t = 0; t < nt; t++)
      {
        if (!dead[(size_t)t]) live.push_back(t);
      }
      for (int pass = 0; pass < 8; pass++)
      {
        if (runFlips() == 0) break;
      }
      removeValenceThree();
      removeFlatApexes();
      for (ll v = 0; v < nv; v++)
      {
        if (!incident[(size_t)v].empty()) pushEdges(v);
      }
    }
    runCollapses();
    // The last collapses leave flat caps of their own (a point snapped onto
    // the zero level in a flat spot lands on the edge across from it), and
    // nothing came after them to take those apart; now the flips and the
    // valence-three removal do, as in every round before.
    live.clear();
    for (ll t = 0; t < nt; t++)
    {
      if (!dead[(size_t)t]) live.push_back(t);
    }
    for (int pass = 0; pass < 8; pass++)
    {
      if (runFlips() == 0) break;
    }
    removeValenceThree();
    removeFlatApexes();
    std::vector<ll> newTris;
    std::vector<ll> newId((size_t)nv, -1);
    std::vector<double> newPts;
    for (ll t = 0; t < nt; t++)
    {
      if (dead[(size_t)t]) continue;
      for (int j = 0; j < 3; j++)
      {
        ll v = tris[(size_t)3*t + j];
        if (newId[(size_t)v] < 0)
        {
          newId[(size_t)v] = (ll)(newPts.size()/3);
          newPts.insert(newPts.end(), &pts[(size_t)3*v], &pts[(size_t)3*v] + 3);
        }
        newTris.push_back(newId[(size_t)v]);
      }
    }
    pts = newPts;
    tris = newTris;
  }
  // Which vessel end each point belongs to (see Surface::pointRim).
  surface.rims = field.Rims();
  surface.pointRim.assign(pts.size()/3, -1);
  for (size_t v = 0; v < pts.size()/3; v++)
  {
    double size, thickness;
    field.Local(&pts[3*v], size, thickness, surface.pointRim[v]);
  }
  auto t4 = std::chrono::steady_clock::now();
  report.secondsDecimate = std::chrono::duration<double>(t4 - t3).count();
  report.numPoints = (ll)(pts.size()/3);
  report.numTriangles = (ll)(tris.size()/3);
  report.numFieldEvaluations = field.NumEvaluations();
  CountEdges(tris, report.numBoundaryEdges, report.numNonManifoldEdges, report.numMiswoundEdges);
  if (report.numNonManifoldEdges > 0 || report.numMiswoundEdges > 0)
  {
    // find the first such edge for the report
    std::map<EdgeKey, std::pair<int, int> > use;
    for (size_t i = 0; i + 2 < tris.size(); i += 3)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = tris[i + j], b = tris[i + (j+1)%3];
        std::pair<int, int> &u = use[EdgeKey(a, b)];
        u.first++;
        if (a < b) u.second++;
      }
    }
    for (std::map<EdgeKey, std::pair<int, int> >::iterator it = use.begin(); it != use.end(); ++it)
    {
      if (it->second.first > 2 || (it->second.first == 2 && it->second.second != 1))
      {
        double at[3];
        for (int k = 0; k < 3; k++)
        {
          at[k] = 0.5*(pts[(size_t)3*it->first.a + k] + pts[(size_t)3*it->first.b + k]);
        }
        NoteFault(report, it->second.first > 2 ? "an edge of the offset surface on more than two triangles" : "an edge of the offset surface traversed the same way by both its triangles", at);
        break;
      }
    }
  }
  }  // the level's value
  }  // each level
  std::vector<double>().swap(cloud);
  std::vector<double>().swap(value);
  return 0;
}

int BuildOffsetSurface(const Interface &input, const Options &options,
    DelaunayFunction delaunay, void *context, Surface &surface, Report &report,
    std::string &error, ProgressFunction progress)
{
  std::vector<double> fractions(1, 1.0);
  std::vector<Surface> surfaces;
  std::vector<Report> reports;
  int rc = BuildOffsetSurfaces(input, options, fractions, delaunay, context, surfaces, reports, error, progress);
  surface = surfaces.empty() ? Surface() : surfaces[0];
  report = reports.empty() ? Report() : reports[0];
  return rc;
}

}  // namespace svoffset
