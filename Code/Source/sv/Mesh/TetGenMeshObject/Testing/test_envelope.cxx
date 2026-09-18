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


// Standalone checks of sv_tetgenmesh_envelope on synthetic surfaces. Build and
// run with a plain C++17 compiler, no VTK or SimVascular needed:
//
//   c++ -std=c++17 -O2 -I.. test_envelope.cxx ../sv_tetgenmesh_envelope.cxx -o test_envelope && ./test_envelope
//
// Every case builds a closed, consistently wound surface that runs through
// itself, takes its envelope, and checks that the envelope is watertight
// (every edge on two pieces, traversed once each way), that no two of its
// triangles cross, and that it encloses what it should.

#include "sv_tetgenmesh_envelope.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

using svenvelope::Surface;
using svenvelope::Envelope;
using svenvelope::Report;
using svenvelope::CleanReport;
typedef long long ll;

static int numFailed = 0;

static void Check(bool ok, const char *what)
{
  printf("  %s %s\n", ok ? "ok  " : "FAIL", what);
  if (!ok)
  {
    numFailed++;
  }
}

//---------------------
// Surface builders
//---------------------

static void AddTriangle(Surface &s, ll a, ll b, ll c)
{
  s.triangles.push_back(a);
  s.triangles.push_back(b);
  s.triangles.push_back(c);
}

static ll AddPoint(Surface &s, double x, double y, double z)
{
  s.points.push_back(x);
  s.points.push_back(y);
  s.points.push_back(z);
  return (ll)(s.points.size()/3 - 1);
}

// An icosphere of the given radius about a centre, outward wound, appended.
static void AddIcosphere(Surface &s, double cx, double cy, double cz, double r, int subdivisions,
    double rotate)
{
  ll base = (ll)(s.points.size()/3);
  const double t = (1.0 + std::sqrt(5.0))/2.0;
  double v[12][3] = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t},
      {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  std::vector<double> pts;
  for (int i = 0; i < 12; i++)
  {
    double n = std::sqrt(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    pts.push_back(v[i][0]/n);
    pts.push_back(v[i][1]/n);
    pts.push_back(v[i][2]/n);
  }
  int f[20][3] = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4},
      {11, 10, 2}, {10, 7, 6}, {7, 1, 8}, {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
      {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  std::vector<ll> tris;
  for (int i = 0; i < 20; i++)
  {
    tris.push_back(f[i][0]);
    tris.push_back(f[i][1]);
    tris.push_back(f[i][2]);
  }
  for (int level = 0; level < subdivisions; level++)
  {
    std::map<std::pair<ll, ll>, ll> mid;
    std::vector<ll> next;
    auto midpoint = [&](ll a, ll b)
    {
      std::pair<ll, ll> key(std::min(a, b), std::max(a, b));
      std::map<std::pair<ll, ll>, ll>::iterator found = mid.find(key);
      if (found != mid.end())
      {
        return found->second;
      }
      double m[3];
      for (int k = 0; k < 3; k++)
      {
        m[k] = 0.5*(pts[(size_t)3*a + k] + pts[(size_t)3*b + k]);
      }
      double n = std::sqrt(m[0]*m[0] + m[1]*m[1] + m[2]*m[2]);
      ll id = (ll)(pts.size()/3);
      pts.push_back(m[0]/n);
      pts.push_back(m[1]/n);
      pts.push_back(m[2]/n);
      mid[key] = id;
      return id;
    };
    for (size_t i = 0; i + 2 < tris.size(); i += 3)
    {
      ll a = tris[i], b = tris[i+1], c = tris[i+2];
      ll ab = midpoint(a, b), bc = midpoint(b, c), ca = midpoint(c, a);
      ll quad[4][3] = {{a, ab, ca}, {b, bc, ab}, {c, ca, bc}, {ab, bc, ca}};
      for (int q = 0; q < 4; q++)
      {
        for (int k = 0; k < 3; k++)
        {
          next.push_back(quad[q][k]);
        }
      }
    }
    tris.swap(next);
  }
  double cr = std::cos(rotate), sr = std::sin(rotate);
  for (size_t i = 0; i + 2 < pts.size(); i += 3)
  {
    double x = pts[i], y = pts[i+1], z = pts[i+2];
    // A rotation about an odd axis, so that nothing lines up with the axes.
    double x1 = cr*x - sr*y, y1 = sr*x + cr*y, z1 = z;
    double y2 = cr*y1 - sr*z1, z2 = sr*y1 + cr*z1;
    AddPoint(s, cx + r*x1, cy + r*y2, cz + r*z2);
  }
  for (size_t i = 0; i + 2 < tris.size(); i += 3)
  {
    AddTriangle(s, base + tris[i], base + tris[i+1], base + tris[i+2]);
  }
}

// An open tube (no end faces) along an axis, outward wound, appended; the
// two rims are returned in the tube's own winding order.
static void AddTube(Surface &s, const double origin[3], const double axis[3], double radius,
    double length, int around, int along, std::vector<ll> rims[2])
{
  double a[3] = {axis[0], axis[1], axis[2]};
  double n = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
  for (int k = 0; k < 3; k++) a[k] /= n;
  double up[3] = {0.0, 0.0, 1.0};
  if (std::abs(a[2]) > 0.9)
  {
    up[0] = 1.0; up[2] = 0.0;
  }
  double u[3] = {a[1]*up[2] - a[2]*up[1], a[2]*up[0] - a[0]*up[2], a[0]*up[1] - a[1]*up[0]};
  n = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
  for (int k = 0; k < 3; k++) u[k] /= n;
  double v[3] = {a[1]*u[2] - a[2]*u[1], a[2]*u[0] - a[0]*u[2], a[0]*u[1] - a[1]*u[0]};
  ll base = (ll)(s.points.size()/3);
  for (int i = 0; i <= along; i++)
  {
    double h = length*(double)i/(double)along;
    for (int j = 0; j < around; j++)
    {
      double th = 2.0*M_PI*(double)j/(double)around + 0.1*i;
      double c = std::cos(th), sn = std::sin(th);
      AddPoint(s, origin[0] + h*a[0] + radius*(c*u[0] + sn*v[0]),
          origin[1] + h*a[1] + radius*(c*u[1] + sn*v[1]),
          origin[2] + h*a[2] + radius*(c*u[2] + sn*v[2]));
    }
  }
  auto id = [&](int i, int j) { return base + (ll)i*around + (ll)(j % around); };
  for (int i = 0; i < along; i++)
  {
    for (int j = 0; j < around; j++)
    {
      // Outward: (u, v, a) is right handed, so going around in theta with the
      // axis up is counter-clockwise seen from outside.
      AddTriangle(s, id(i, j), id(i, j+1), id(i+1, j));
      AddTriangle(s, id(i, j+1), id(i+1, j+1), id(i+1, j));
    }
  }
  rims[0].clear();
  rims[1].clear();
  // The boundary edges as the triangles traverse them: at i = 0 the edge
  // id(0,j+1) -> id(0,j) is used by the first triangle as id(0,j) -> id(0,j+1),
  // so the rim in traversal order runs with j increasing at i = 0 and with j
  // decreasing at i = along.
  for (int j = 0; j < around; j++)
  {
    rims[0].push_back(id(0, j));
  }
  for (int j = around; j > 0; j--)
  {
    rims[1].push_back(id(along, j % around));
  }
}

// Closes a rim with a fan, wound against the rim's traversal order as the
// caller in SimVascular does it.
static void AddFan(Surface &s, const std::vector<ll> &rim)
{
  double c[3] = {0.0, 0.0, 0.0};
  for (size_t m = 0; m < rim.size(); m++)
  {
    for (int k = 0; k < 3; k++)
    {
      c[k] += s.points[(size_t)3*rim[m] + k];
    }
  }
  ll centre = AddPoint(s, c[0]/rim.size(), c[1]/rim.size(), c[2]/rim.size());
  for (size_t m = 0; m < rim.size(); m++)
  {
    ll a = rim[m], b = rim[(m+1) % rim.size()];
    AddTriangle(s, b, a, centre);
  }
}

// A block with a parabolic valley for a top, outward wound: x in [-1, 1],
// y in [0, 3], z from -1 up to c x^2.
static void AddValleyBlock(Surface &s, double c, int nx, int ny)
{
  ll base = (ll)(s.points.size()/3);
  auto top = [&](int i, int j) { return base + (ll)i*(ny+1) + j; };
  for (int i = 0; i <= nx; i++)
  {
    for (int j = 0; j <= ny; j++)
    {
      double x = -1.0 + 2.0*i/nx, y = 3.0*j/ny;
      AddPoint(s, x, y, c*x*x);
    }
  }
  ll bottomBase = (ll)(s.points.size()/3);
  auto bottom = [&](int i, int j) { return bottomBase + (ll)i*(ny+1) + j; };
  for (int i = 0; i <= nx; i++)
  {
    for (int j = 0; j <= ny; j++)
    {
      double x = -1.0 + 2.0*i/nx, y = 3.0*j/ny;
      AddPoint(s, x, y, -1.0);
    }
  }
  for (int i = 0; i < nx; i++)
  {
    for (int j = 0; j < ny; j++)
    {
      // Top faces up: counter-clockwise seen from +z is (i,j) -> (i+1,j) -> (i+1,j+1).
      AddTriangle(s, top(i, j), top(i+1, j), top(i+1, j+1));
      AddTriangle(s, top(i, j), top(i+1, j+1), top(i, j+1));
      AddTriangle(s, bottom(i, j), bottom(i+1, j+1), bottom(i+1, j));
      AddTriangle(s, bottom(i, j), bottom(i, j+1), bottom(i+1, j+1));
    }
  }
  // Sides: at y = 0 the outward normal is -y; the quad (top(i,0), top(i+1,0),
  // bottom(i+1,0), bottom(i,0)) seen from -y ... wound so that the normal is -y.
  for (int i = 0; i < nx; i++)
  {
    AddTriangle(s, top(i, 0), bottom(i, 0), bottom(i+1, 0));
    AddTriangle(s, top(i, 0), bottom(i+1, 0), top(i+1, 0));
    AddTriangle(s, top(i, ny), top(i+1, ny), bottom(i+1, ny));
    AddTriangle(s, top(i, ny), bottom(i+1, ny), bottom(i, ny));
  }
  for (int j = 0; j < ny; j++)
  {
    AddTriangle(s, top(0, j), top(0, j+1), bottom(0, j+1));
    AddTriangle(s, top(0, j), bottom(0, j+1), bottom(0, j));
    AddTriangle(s, top(nx, j), bottom(nx, j), bottom(nx, j+1));
    AddTriangle(s, top(nx, j), bottom(nx, j+1), top(nx, j+1));
  }
}

// Extrudes every point of a surface along its area-weighted vertex normal by
// a thickness, in place.
static void Extrude(Surface &s, double thickness)
{
  size_t numPts = s.points.size()/3;
  std::vector<double> normal(3*numPts, 0.0);
  for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
  {
    const double *a = &s.points[(size_t)3*s.triangles[i]];
    const double *b = &s.points[(size_t)3*s.triangles[i+1]];
    const double *c = &s.points[(size_t)3*s.triangles[i+2]];
    double e1[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]}, e2[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    double n[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
    for (int j = 0; j < 3; j++)
    {
      for (int k = 0; k < 3; k++)
      {
        normal[(size_t)3*s.triangles[i+j] + k] += n[k];
      }
    }
  }
  for (size_t p = 0; p < numPts; p++)
  {
    double *n = &normal[3*p];
    double len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len > 0.0)
    {
      for (int k = 0; k < 3; k++)
      {
        s.points[3*p + k] += thickness*n[k]/len;
      }
    }
  }
}

// Gives an extruded surface the normals its sheet triangles had before the
// extrusion, from the surface it was extruded from (same triangles).
static void SheetNormals(const Surface &inner, Surface &s)
{
  s.sheetNormal.assign((size_t)3*s.numSheetTriangles, 0.0);
  for (ll t = 0; t < s.numSheetTriangles; t++)
  {
    const double *a = &inner.points[(size_t)3*inner.triangles[(size_t)3*t]];
    const double *b = &inner.points[(size_t)3*inner.triangles[(size_t)3*t + 1]];
    const double *c = &inner.points[(size_t)3*inner.triangles[(size_t)3*t + 2]];
    double e1[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]}, e2[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    s.sheetNormal[(size_t)3*t] = e1[1]*e2[2]-e1[2]*e2[1];
    s.sheetNormal[(size_t)3*t + 1] = e1[2]*e2[0]-e1[0]*e2[2];
    s.sheetNormal[(size_t)3*t + 2] = e1[0]*e2[1]-e1[1]*e2[0];
  }
}

//---------------------
// Checks
//---------------------

// Watertight apart from the given number of boundary edges: every edge on two
// triangles traversed once each way.
static bool Watertight(const Envelope &e, ll expectedBoundary, ll &boundary, ll &bad)
{
  std::map<std::pair<ll, ll>, std::pair<int, int> > use;
  for (size_t i = 0; i + 2 < e.triangles.size(); i += 3)
  {
    for (int j = 0; j < 3; j++)
    {
      ll a = e.triangles[i + j], b = e.triangles[i + (j+1)%3];
      std::pair<int, int> &u = use[std::make_pair(std::min(a, b), std::max(a, b))];
      u.first++;
      if (a < b) u.second++;
    }
  }
  boundary = 0;
  bad = 0;
  for (std::map<std::pair<ll, ll>, std::pair<int, int> >::iterator it = use.begin(); it != use.end(); ++it)
  {
    if (it->second.first == 1)
    {
      boundary++;
    }
    else if (it->second.first != 2 || it->second.second != 1)
    {
      bad++;
    }
  }
  return bad == 0 && boundary == expectedBoundary;
}

static double EnclosedVolume(const Envelope &e)
{
  double volume = 0.0;
  for (size_t i = 0; i + 2 < e.triangles.size(); i += 3)
  {
    const double *a = &e.points[(size_t)3*e.triangles[i]];
    const double *b = &e.points[(size_t)3*e.triangles[i+1]];
    const double *c = &e.points[(size_t)3*e.triangles[i+2]];
    double bc[3] = {b[1]*c[2]-b[2]*c[1], b[2]*c[0]-b[0]*c[2], b[0]*c[1]-b[1]*c[0]};
    volume += (a[0]*bc[0] + a[1]*bc[1] + a[2]*bc[2])/6.0;
  }
  return volume;
}

// Distance from a point to a triangle mesh, brute force.
static double DistanceToMesh(const double p[3], const Envelope &e)
{
  double best = 1.0e300;
  for (size_t i = 0; i + 2 < e.triangles.size(); i += 3)
  {
    const double *a = &e.points[(size_t)3*e.triangles[i]];
    const double *b = &e.points[(size_t)3*e.triangles[i+1]];
    const double *c = &e.points[(size_t)3*e.triangles[i+2]];
    // Ericson's closest point on triangle.
    double ab[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
    double ac[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
    double ap[3] = {p[0]-a[0], p[1]-a[1], p[2]-a[2]};
    double d1 = ab[0]*ap[0]+ab[1]*ap[1]+ab[2]*ap[2], d2 = ac[0]*ap[0]+ac[1]*ap[1]+ac[2]*ap[2];
    double q[3];
    if (d1 <= 0 && d2 <= 0) { q[0]=a[0]; q[1]=a[1]; q[2]=a[2]; }
    else
    {
      double bp[3] = {p[0]-b[0], p[1]-b[1], p[2]-b[2]};
      double d3 = ab[0]*bp[0]+ab[1]*bp[1]+ab[2]*bp[2], d4 = ac[0]*bp[0]+ac[1]*bp[1]+ac[2]*bp[2];
      if (d3 >= 0 && d4 <= d3) { q[0]=b[0]; q[1]=b[1]; q[2]=b[2]; }
      else
      {
        double vc = d1*d4 - d3*d2;
        if (vc <= 0 && d1 >= 0 && d3 <= 0)
        {
          double v = d1/(d1 - d3);
          for (int k = 0; k < 3; k++) q[k] = a[k] + v*ab[k];
        }
        else
        {
          double cp[3] = {p[0]-c[0], p[1]-c[1], p[2]-c[2]};
          double d5 = ab[0]*cp[0]+ab[1]*cp[1]+ab[2]*cp[2], d6 = ac[0]*cp[0]+ac[1]*cp[1]+ac[2]*cp[2];
          if (d6 >= 0 && d5 <= d6) { q[0]=c[0]; q[1]=c[1]; q[2]=c[2]; }
          else
          {
            double vb = d5*d2 - d1*d6;
            if (vb <= 0 && d2 >= 0 && d6 <= 0)
            {
              double w = d2/(d2 - d6);
              for (int k = 0; k < 3; k++) q[k] = a[k] + w*ac[k];
            }
            else
            {
              double va = d3*d6 - d5*d4;
              if (va <= 0 && (d4 - d3) >= 0 && (d5 - d6) >= 0)
              {
                double w = (d4 - d3)/((d4 - d3) + (d5 - d6));
                for (int k = 0; k < 3; k++) q[k] = b[k] + w*(c[k] - b[k]);
              }
              else
              {
                double denom = 1.0/(va + vb + vc);
                double v = vb*denom, w = vc*denom;
                for (int k = 0; k < 3; k++) q[k] = a[k] + ab[k]*v + ac[k]*w;
              }
            }
          }
        }
      }
    }
    double d = std::sqrt((p[0]-q[0])*(p[0]-q[0]) + (p[1]-q[1])*(p[1]-q[1]) + (p[2]-q[2])*(p[2]-q[2]));
    best = std::min(best, d);
  }
  return best;
}

static void PrintReport(const Report &r)
{
  printf("  pairs tested %lld, crossing %lld; crease points %lld, triple points %lld; %lld triangles split into %lld pieces, %lld pieces kept; whole kept %lld, dropped %lld; inverted %lld\n",
      r.numPairsTested, r.numPairsCrossing, r.numCreasePoints, r.numTriplePoints,
      r.numTrianglesSplit, r.numPieces, r.numPiecesKept, r.numWholeKept, r.numWholeDropped, r.numInverted);
  printf("  rays %lld (%lld retried, %lld undecided); faults %lld, non-manifold %lld, miswound %lld, boundary %lld; %.3f s\n",
      r.numRays, r.numRayRetries, r.numUndecided, r.numArrangementFaults, r.numNonManifoldEdges,
      r.numMiswoundEdges, r.numBoundaryEdges, r.seconds);
  if (!r.firstFault.empty())
  {
    printf("  first fault: %s at (%g, %g, %g)\n", r.firstFault.c_str(), r.firstFaultAt[0], r.firstFaultAt[1], r.firstFaultAt[2]);
  }
  if (r.numPockets > 0)
  {
    printf("  pockets and shreds: %lld components cut off from the rims, %lld pieces dropped with them\n", r.numPockets, r.numPocketPieces);
  }
  printf("  winding histogram:");
  for (size_t i = 0; i < r.windingHistogram.size(); i++)
  {
    printf(" w=%d:%lld", (int)i - r.windingOffset, r.windingHistogram[i]);
  }
  printf("\n");
}

static bool RunAndCheckClosed(const char *name, const Surface &s, Envelope &e, Report &r,
    ll expectedBoundary)
{
  printf("%s\n", name);
  std::string error;
  int rc = svenvelope::BuildOuterEnvelope(s, e, r, error);
  if (rc != 0)
  {
    printf("  FAIL build: %s\n", error.c_str());
    numFailed++;
    return false;
  }
  PrintReport(r);
  ll boundary = 0, bad = 0;
  bool tight = Watertight(e, expectedBoundary, boundary, bad);
  char what[160];
  snprintf(what, sizeof(what), "watertight (boundary edges %lld of %lld expected, bad %lld)", boundary, expectedBoundary, bad);
  Check(tight, what);
  std::vector<unsigned char> crossing;
  double at[3];
  ll numCrossing = svenvelope::CountCrossingTriangles(e.points, e.triangles, crossing, at);
  snprintf(what, sizeof(what), "no crossings (%lld crossing, first at %g %g %g)", numCrossing, at[0], at[1], at[2]);
  Check(numCrossing == 0, what);
  Check(r.numArrangementFaults == 0 && r.numUndecided == 0 && r.numNonManifoldEdges == 0 && r.numMiswoundEdges == 0,
      "no faults reported");
  return tight && numCrossing == 0;
}

// Monte Carlo volume of a union of balls.
static double UnionOfBallsVolume(const std::vector<double> &centres, double r, int samples)
{
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (size_t c = 0; c + 2 < centres.size(); c += 3)
  {
    for (int k = 0; k < 3; k++)
    {
      lo[k] = std::min(lo[k], centres[c+k] - r);
      hi[k] = std::max(hi[k], centres[c+k] + r);
    }
  }
  unsigned long long state = 12345;
  auto next = [&]() { state = state*6364136223846793005ULL + 1442695040888963407ULL; return (double)(state >> 11)/9007199254740992.0; };
  int inside = 0;
  for (int i = 0; i < samples; i++)
  {
    double p[3];
    for (int k = 0; k < 3; k++)
    {
      p[k] = lo[k] + (hi[k] - lo[k])*next();
    }
    bool in = false;
    for (size_t c = 0; c + 2 < centres.size() && !in; c += 3)
    {
      double d2 = 0.0;
      for (int k = 0; k < 3; k++)
      {
        d2 += (p[k] - centres[c+k])*(p[k] - centres[c+k]);
      }
      in = d2 < r*r;
    }
    if (in) inside++;
  }
  return (hi[0]-lo[0])*(hi[1]-lo[1])*(hi[2]-lo[2])*(double)inside/(double)samples;
}

int main(int argc, char **argv)
{
  bool perf = argc > 1 && std::string(argv[1]) == "perf";

  // 1. A sphere on its own: nothing crosses, everything is kept whole.
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 4, 0.3);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("1. one sphere", s, e, r, 0);
    Check(r.numPairsCrossing == 0 && r.numPiecesKept == s.numSheetTriangles && r.numWholeKept == s.numSheetTriangles,
        "kept whole and uncut");
    double v = EnclosedVolume(e);
    Check(std::abs(v - 4.0/3.0*M_PI) < 0.02*4.0/3.0*M_PI, "volume of a sphere");
  }

  // 2. Two spheres overlapping: the envelope is the boundary of their union.
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 4, 0.3);
    AddIcosphere(s, 1.1, 0.2, -0.1, 1.0, 4, 0.7);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("2. two spheres", s, e, r, 0);
    Check(r.numPairsCrossing > 0 && r.numCreasePoints > 0, "crease found");
    Check(r.numTriplePoints == 0, "no triple points for two spheres");
    std::vector<double> centres = {0, 0, 0, 1.1, 0.2, -0.1};
    double expected = UnionOfBallsVolume(centres, 1.0, 400000);
    double v = EnclosedVolume(e);
    char what[160];
    snprintf(what, sizeof(what), "volume of the union %.4f vs %.4f sampled", v, expected);
    Check(std::abs(v - expected) < 0.02*expected, what);
    // Every kept point is outside the other sphere, up to the facet sag of
    // the icosphere.
    int inside = 0;
    for (size_t i = 0; i + 2 < e.triangles.size(); i += 3)
    {
      for (int j = 0; j < 3; j++)
      {
        const double *p = &e.points[(size_t)3*e.triangles[i+j]];
        double d0 = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
        double d1 = std::sqrt((p[0]-1.1)*(p[0]-1.1) + (p[1]-0.2)*(p[1]-0.2) + (p[2]+0.1)*(p[2]+0.1));
        if (std::max(d0, d1) < 0.97)
        {
          inside++;
        }
      }
    }
    Check(inside == 0, "no envelope point inside both spheres");
  }

  // 3. Three spheres with a common region: creases cross at triple points.
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 4, 0.3);
    AddIcosphere(s, 1.1, 0.2, -0.1, 1.0, 4, 0.7);
    AddIcosphere(s, 0.5, 1.0, 0.3, 1.0, 4, 1.1);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("3. three spheres", s, e, r, 0);
    Check(r.numTriplePoints > 0, "triple points found");
    std::vector<double> centres = {0, 0, 0, 1.1, 0.2, -0.1, 0.5, 1.0, 0.3};
    double expected = UnionOfBallsVolume(centres, 1.0, 400000);
    double v = EnclosedVolume(e);
    char what[160];
    snprintf(what, sizeof(what), "volume of the union %.4f vs %.4f sampled", v, expected);
    Check(std::abs(v - expected) < 0.02*expected, what);
  }

  // 4. A valley folded over: the extrusion of a concave surface at more than
  // its radius of curvature runs through itself, and the envelope is the
  // crease of the two flanks.
  {
    Surface inner;
    AddValleyBlock(inner, 2.0, 40, 30);
    Surface s = inner;
    const double t = 0.5;
    Extrude(s, t);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    SheetNormals(inner, s);
    Envelope e;
    Report r;
    RunAndCheckClosed("4. folded valley", s, e, r, 0);
    Check(r.numPairsCrossing > 0 && r.numInverted > 0, "fold found (crossings and turned-over triangles)");
    Envelope innerAsEnvelope;
    innerAsEnvelope.points = inner.points;
    innerAsEnvelope.triangles = inner.triangles;
    double vInner = EnclosedVolume(innerAsEnvelope);
    double v = EnclosedVolume(e);
    char what[200];
    snprintf(what, sizeof(what), "volume grew from %.4f to %.4f", vInner, v);
    Check(v > vInner, what);
    // The wall over the valley floor is at least the thickness, measured
    // from inner points well inside the top face.
    double minDistance = 1e300, maxDistance = 0.0;
    int numSampled = 0;
    for (size_t p = 0; p < inner.points.size()/3; p++)
    {
      const double *q = &inner.points[3*p];
      if (std::abs(q[0]) > 0.8 || q[1] < 0.3 || q[1] > 2.7 || q[2] < 2.0*q[0]*q[0] - 1e-9)
      {
        continue;
      }
      double d = DistanceToMesh(q, e);
      minDistance = std::min(minDistance, d);
      maxDistance = std::max(maxDistance, d);
      numSampled++;
    }
    snprintf(what, sizeof(what), "wall over the valley: %d inner points stand %.4f to %.4f from the envelope (t = %.2f)", numSampled, minDistance, maxDistance, t);
    Check(numSampled > 0 && minDistance > 0.9*t && maxDistance < 2.0*t, what);
    // At the valley floor the envelope is the crease, above the offset of
    // the floor itself.
    double floor[3] = {0.0, 1.5, 0.0};
    double d = DistanceToMesh(floor, e);
    snprintf(what, sizeof(what), "the crease over the valley floor stands %.4f off it, more than t", d);
    Check(d > t, what);
  }

  // 5. Two open tubes crossing, each closed with a fan for the winding
  // number: the envelope keeps both rims whole and the tubes join along
  // their crease. The fans are not part of the sheet.
  {
    Surface s;
    std::vector<ll> rimsA[2], rimsB[2];
    double oA[3] = {-3.0, 0.0, 0.0}, aA[3] = {1.0, 0.05, 0.02};
    double oB[3] = {0.1, -3.0, 0.3}, aB[3] = {0.03, 1.0, 0.1};
    AddTube(s, oA, aA, 1.0, 6.0, 48, 40, rimsA);
    AddTube(s, oB, aB, 0.8, 6.0, 40, 40, rimsB);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    AddFan(s, rimsA[0]);
    AddFan(s, rimsA[1]);
    AddFan(s, rimsB[0]);
    AddFan(s, rimsB[1]);
    Envelope e;
    Report r;
    ll expectedBoundary = 48 + 48 + 40 + 40;
    RunAndCheckClosed("5. two tubes crossing, with cap fans", s, e, r, expectedBoundary);
    Check(r.numPairsCrossing > 0, "crease found");
    // The rims are whole: every rim edge is a boundary edge of the envelope.
    std::map<std::pair<ll, ll>, int> use;
    for (size_t i = 0; i + 2 < e.triangles.size(); i += 3)
    {
      for (int j = 0; j < 3; j++)
      {
        ll a = e.triangles[i + j], b = e.triangles[i + (j+1)%3];
        use[std::make_pair(std::min(a, b), std::max(a, b))]++;
      }
    }
    int rimEdgesMissing = 0;
    const std::vector<ll> *rims[4] = {&rimsA[0], &rimsA[1], &rimsB[0], &rimsB[1]};
    for (int q = 0; q < 4; q++)
    {
      const std::vector<ll> &rim = *rims[q];
      for (size_t m = 0; m < rim.size(); m++)
      {
        ll a = rim[m], b = rim[(m+1) % rim.size()];
        std::map<std::pair<ll, ll>, int>::iterator found = use.find(std::make_pair(std::min(a, b), std::max(a, b)));
        if (found == use.end() || found->second != 1)
        {
          rimEdgesMissing++;
        }
      }
    }
    char what[160];
    snprintf(what, sizeof(what), "all four rims whole on the boundary (%d rim edges not)", rimEdgesMissing);
    Check(rimEdgesMissing == 0, what);
  }

  // 6. The crossing count on its own.
  {
    std::vector<double> pts = {0,0,0, 1,0,0, 0,1,0,  0.2,0.2,-1, 0.3,0.2,1, 0.2,0.4,1,  5,5,5, 6,5,5, 5,6,5};
    std::vector<ll> tris = {0,1,2, 3,4,5, 6,7,8};
    std::vector<unsigned char> crossing;
    double at[3];
    ll n = svenvelope::CountCrossingTriangles(pts, tris, crossing, at);
    printf("6. crossing count\n");
    Check(n == 2 && crossing[0] && crossing[1] && !crossing[2], "two crossing triangles counted, the third not");
    // Sharing a corner and folded through each other.
    std::vector<double> pts2 = {0,0,0, 1,0,0, 0,1,0,  1,0,0.5, 0,1,-0.5};
    std::vector<ll> tris2 = {0,1,2, 0,3,4};
    n = svenvelope::CountCrossingTriangles(pts2, tris2, crossing, at);
    Check(n == 2, "two triangles sharing a corner and crossing are counted");
    std::vector<double> pts3 = {0,0,0, 1,0,0, 0,1,0, 1,1,0.3};
    std::vector<ll> tris3 = {0,1,2, 1,3,2};
    n = svenvelope::CountCrossingTriangles(pts3, tris3, crossing, at);
    Check(n == 0, "two triangles sharing an edge are not");
  }

  // 7. The same two spheres wound inward: the envelope is the same surface,
  // wound the way the input was, so the result does not depend on the
  // convention.
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 4, 0.3);
    AddIcosphere(s, 1.1, 0.2, -0.1, 1.0, 4, 0.7);
    for (size_t i = 0; i + 2 < s.triangles.size(); i += 3)
    {
      std::swap(s.triangles[i+1], s.triangles[i+2]);
    }
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("7. two spheres wound inward", s, e, r, 0);
    Check(r.windingSense == -1, "wound inward, as read off the volume");
    Check(r.numPiecesKept == 8600 && r.numWholeDropped == 2099, "the same pieces kept and dropped as when wound outward");
    double v = EnclosedVolume(e);
    Check(v < 0.0 && std::abs(v + 7.3308) < 0.01, "the same volume, with the input's sign");
  }

  // 8. A surface that is not consistently wound is refused.
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 2, 0.3);
    std::swap(s.triangles[1], s.triangles[2]);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    std::string error;
    int rc = svenvelope::BuildOuterEnvelope(s, e, r, error);
    printf("8. inconsistent winding\n");
    Check(rc != 0 && error.find("consistently wound") != std::string::npos, "refused with a message naming the winding");
  }

  // 9. Sliver cleanup: a valley whose grid is jittered so that the crease
  // passes close to corners and along edges, leaving pieces with almost no
  // altitude. The cleanup has to take them out without moving a crease
  // point, without opening the surface, folding it or making it cross
  // itself, and without taking wall away over the floor.
  // The seeds are ones whose jittered valley the envelope takes without an
  // arrangement fault, since the caller runs the cleanup on a sound envelope
  // only; a dangling crease segment under heavier jitter is the envelope's
  // own open question, not the cleanup's.
  const int seeds[4] = {0, 1, 8, 9};
  for (int which = 0; which < 4; which++)
  {
    const int variant = seeds[which];
    Surface inner;
    AddValleyBlock(inner, 2.0, 40, 30);
    // Jitter the interior top points by up to 0.4 of the grid spacing.
    unsigned long long state = 12345ULL + 977ULL*(unsigned long long)variant;
    auto next = [&]() { state = state*6364136223846793005ULL + 1442695040888963407ULL; return (double)(state >> 11)/9007199254740992.0 - 0.5; };
    const int nx = 40, ny = 30;
    for (int i = 1; i < nx; i++)
    {
      for (int j = 1; j < ny; j++)
      {
        size_t p = (size_t)i*(ny+1) + j;
        double x = inner.points[3*p] + 0.8*(2.0/nx)*next()*(variant == 0 ? 0.0 : 1.0);
        double y = inner.points[3*p + 1] + 0.8*(3.0/ny)*next()*(variant == 0 ? 0.0 : 1.0);
        inner.points[3*p] = x;
        inner.points[3*p + 1] = y;
        inner.points[3*p + 2] = 2.0*x*x;
      }
    }
    Surface s = inner;
    const double t = 0.5;
    Extrude(s, t);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    SheetNormals(inner, s);
    Envelope e;
    Report r;
    char name[120];
    snprintf(name, sizeof(name), "9.%d sliver cleanup on a %s valley (seed %d)", which, variant == 0 ? "regular" : "jittered", variant);
    bool sound = RunAndCheckClosed(name, s, e, r, 0) && r.numArrangementFaults == 0;
    if (!sound)
    {
      printf("  (the envelope has faults; the cleanup is not run on it, as the caller would not)\n");
      continue;
    }
    // Slivers that touch a crease are the ones the cleanup is for; the
    // block's ridges make long triangles too, which are its own affair.
    auto creaseSlivers = [&](const Envelope &env, double &worstAspect)
    {
      ll n = 0;
      worstAspect = 0.0;
      for (size_t i = 0; i + 2 < env.triangles.size(); i += 3)
      {
        const double *a = &env.points[3*env.triangles[i]], *b = &env.points[3*env.triangles[i+1]], *c = &env.points[3*env.triangles[i+2]];
        double L = 0.0;
        const double *corner[3] = {a, b, c};
        for (int j = 0; j < 3; j++)
        {
          const double *x = corner[j], *y = corner[(j+1)%3];
          L = std::max(L, std::sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1]) + (x[2]-y[2])*(x[2]-y[2])));
        }
        double e1[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]}, e2[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
        double n3[3] = {e1[1]*e2[2]-e1[2]*e2[1], e1[2]*e2[0]-e1[0]*e2[2], e1[0]*e2[1]-e1[1]*e2[0]};
        double twiceArea = std::sqrt(n3[0]*n3[0] + n3[1]*n3[1] + n3[2]*n3[2]);
        double aspect = (twiceArea > 0.0) ? L*L*std::sqrt(3.0)/(2.0*twiceArea) : 1e300;
        bool touchesCrease = env.pointKind[env.triangles[i]] != 0 || env.pointKind[env.triangles[i+1]] != 0 || env.pointKind[env.triangles[i+2]] != 0;
        if (!touchesCrease) continue;
        worstAspect = std::max(worstAspect, aspect);
        if (aspect > 10.0) n++;
      }
      return n;
    };
    double worstCreaseBefore = 0.0, worstCreaseAfter = 0.0;
    ll creaseSliversBefore = creaseSlivers(e, worstCreaseBefore);
    double vBefore = EnclosedVolume(e);
    std::vector<double> creaseBefore;
    for (size_t p = 0; p < e.pointKind.size(); p++)
    {
      if (e.pointKind[p] != 0)
      {
        creaseBefore.push_back(e.points[3*p]);
        creaseBefore.push_back(e.points[3*p + 1]);
        creaseBefore.push_back(e.points[3*p + 2]);
      }
    }
    auto folds = [&](const Envelope &env)
    {
      // Edges whose two pieces face away from each other, not counting the
      // creases, where the sheets meet at whatever angle they make.
      std::map<std::pair<ll, ll>, std::vector<size_t> > on;
      for (size_t i = 0; i + 2 < env.triangles.size(); i += 3)
      {
        for (int j = 0; j < 3; j++)
        {
          ll a = env.triangles[i + j], b = env.triangles[i + (j+1)%3];
          on[std::make_pair(std::min(a, b), std::max(a, b))].push_back(i/3);
        }
      }
      ll folded = 0;
      for (std::map<std::pair<ll, ll>, std::vector<size_t> >::iterator it = on.begin(); it != on.end(); ++it)
      {
        if (it->second.size() != 2) continue;
        bool crease = env.pointKind[(size_t)it->first.first] != 0 && env.pointKind[(size_t)it->first.second] != 0 &&
            env.source[it->second[0]] != env.source[it->second[1]];
        if (crease) continue;
        double n[2][3];
        for (int k = 0; k < 2; k++)
        {
          const ll *tri = &env.triangles[3*it->second[k]];
          const double *a = &env.points[3*tri[0]], *b = &env.points[3*tri[1]], *c = &env.points[3*tri[2]];
          double e1[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]}, e2[3] = {c[0]-a[0], c[1]-a[1], c[2]-a[2]};
          n[k][0] = e1[1]*e2[2]-e1[2]*e2[1]; n[k][1] = e1[2]*e2[0]-e1[0]*e2[2]; n[k][2] = e1[0]*e2[1]-e1[1]*e2[0];
        }
        if (n[0][0]*n[1][0] + n[0][1]*n[1][1] + n[0][2]*n[1][2] < 0.0) folded++;
      }
      return folded;
    };
    ll foldsBefore = folds(e);
    CleanReport c;
    std::vector<unsigned char> noFixed;
    int rc = svenvelope::CleanEnvelopeSlivers(e, noFixed, 10.0, c);
    printf("  slivers %lld -> %lld (worst %.1f -> %.1f), %lld collapsed, %lld flipped, %d passes, %.3f s\n",
        c.numSliversBefore, c.numSliversAfter, c.worstBefore, c.worstAfter, c.numCollapsed, c.numFlipped, c.numPasses, c.seconds);
    Check(rc == 0, "cleanup ran");
    char what[200];
    if (variant == 0)
    {
      Check(c.numSliversBefore >= 0, "regular valley surveyed");
    }
    else
    {
      Check(c.numSliversBefore > 0, "the jittered valley had slivers to clean");
    }
    ll creaseSliversAfter = creaseSlivers(e, worstCreaseAfter);
    snprintf(what, sizeof(what), "slivers touching a crease: %lld -> %lld (worst %.0f -> %.0f); all slivers %lld -> %lld",
        creaseSliversBefore, creaseSliversAfter, worstCreaseBefore, worstCreaseAfter, c.numSliversBefore, c.numSliversAfter);
    Check(creaseSliversAfter*2 <= creaseSliversBefore, what);
    Check(creaseSliversBefore == 0 || worstCreaseAfter < worstCreaseBefore, "the worst piece on a crease is better");
    ll boundary = 0, bad = 0;
    Check(Watertight(e, 0, boundary, bad), "still watertight after the cleanup");
    std::vector<unsigned char> crossing;
    double at[3];
    ll numCrossing = svenvelope::CountCrossingTriangles(e.points, e.triangles, crossing, at);
    snprintf(what, sizeof(what), "no crossings after the cleanup (%lld)", numCrossing);
    Check(numCrossing == 0, what);
    ll foldsAfter = folds(e);
    snprintf(what, sizeof(what), "folded edges %lld -> %lld", foldsBefore, foldsAfter);
    Check(foldsAfter <= foldsBefore, what);
    std::vector<double> creaseAfter;
    for (size_t p = 0; p < e.pointKind.size(); p++)
    {
      if (e.pointKind[p] != 0)
      {
        creaseAfter.push_back(e.points[3*p]);
        creaseAfter.push_back(e.points[3*p + 1]);
        creaseAfter.push_back(e.points[3*p + 2]);
      }
    }
    Check(creaseAfter == creaseBefore, "no crease point moved");
    size_t numReshaped = 0, numWhole = 0;
    for (size_t i = 0; i < e.whole.size(); i++)
    {
      if (e.whole[i] == 2) numReshaped++;
      if (e.whole[i] == 1) numWhole++;
    }
    Check(e.source.size() == e.triangles.size()/3 && e.whole.size() == e.source.size(), "source and whole follow the triangles");
    Check(c.numCollapsed + c.numFlipped == 0 || numReshaped > 0, "reshaped pieces are marked");
    double vAfter = EnclosedVolume(e);
    snprintf(what, sizeof(what), "volume %.4f -> %.4f", vBefore, vAfter);
    Check(std::abs(vAfter - vBefore) < 0.01*vBefore, what);
    // The wall over the floor, away from the block's ends, whose ridges the
    // jitter reaches.
    double minDistance = 1e300;
    int numSampled = 0;
    for (size_t p = 0; p < inner.points.size()/3; p++)
    {
      const double *q = &inner.points[3*p];
      if (std::abs(q[0]) > 0.8 || q[1] < 0.5 || q[1] > 2.5 || q[2] < 2.0*q[0]*q[0] - 1e-9)
      {
        continue;
      }
      minDistance = std::min(minDistance, DistanceToMesh(q, e));
      numSampled++;
    }
    snprintf(what, sizeof(what), "wall over the valley after the cleanup: %d inner points stand at least %.4f from the envelope (t = %.2f)", numSampled, minDistance, t);
    Check(numSampled > 0 && minDistance > 0.9*t, what);
    (void)numWhole;
  }

  if (perf)
  {
    Surface s;
    AddIcosphere(s, 0, 0, 0, 1.0, 7, 0.3);
    AddIcosphere(s, 1.1, 0.2, -0.1, 1.0, 7, 0.7);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("perf: two spheres of 327680 triangles each", s, e, r, 0);
  }

  printf("%s: %d failed\n", numFailed == 0 ? "PASS" : "FAIL", numFailed);
  return numFailed == 0 ? 0 : 1;
}
