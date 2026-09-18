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
#include <array>
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
// The slab [-5, 5] x [-5, 5] x [-1, 0], twelve triangles wound outward; its
// top is the two triangles (t0, t1, t2), the half x > y, and (t0, t2, t3).
static void AddSlab(Surface &s)
{
  ll b0 = AddPoint(s, -5, -5, -1), b1 = AddPoint(s, 5, -5, -1), b2 = AddPoint(s, 5, 5, -1), b3 = AddPoint(s, -5, 5, -1);
  ll t0 = AddPoint(s, -5, -5, 0), t1 = AddPoint(s, 5, -5, 0), t2 = AddPoint(s, 5, 5, 0), t3 = AddPoint(s, -5, 5, 0);
  AddTriangle(s, t0, t1, t2); AddTriangle(s, t0, t2, t3);          // top, up
  AddTriangle(s, b0, b2, b1); AddTriangle(s, b0, b3, b2);          // bottom, down
  AddTriangle(s, b0, b1, t1); AddTriangle(s, b0, t1, t0);          // -y
  AddTriangle(s, b1, b2, t2); AddTriangle(s, b1, t2, t1);          // +x
  AddTriangle(s, b2, b3, t3); AddTriangle(s, b2, t3, t2);          // +y
  AddTriangle(s, b3, b0, t0); AddTriangle(s, b3, t0, t3);          // -x
}

// A tetrahedron wound outward: a base triangle at height zb, counter-clockwise
// seen from above, and an apex.
static void AddTetrahedron(Surface &s, const double b0[2], const double b1[2], const double b2[2], double zb, const double apex[3])
{
  ll a0 = AddPoint(s, b0[0], b0[1], zb), a1 = AddPoint(s, b1[0], b1[1], zb), a2 = AddPoint(s, b2[0], b2[1], zb);
  ll ap = AddPoint(s, apex[0], apex[1], apex[2]);
  AddTriangle(s, a0, a2, a1);
  AddTriangle(s, a0, a1, ap); AddTriangle(s, a1, a2, ap); AddTriangle(s, a2, a0, ap);
}

// The valley block of AddValleyBlock with its inner grid points jittered by
// up to 0.4 of a cell in x and y (and z following the valley), by a seed;
// seed 0 is the regular grid.
static void AddJitteredValley(Surface &inner, int nx, int ny, int seed)
{
  AddValleyBlock(inner, 2.0, nx, ny);
  unsigned long long state = 12345ULL + 977ULL*(unsigned long long)seed;
  auto next = [&]() { state = state*6364136223846793005ULL + 1442695040888963407ULL; return (double)(state >> 11)/9007199254740992.0 - 0.5; };
  for (int i = 1; i < nx; i++)
  {
    for (int j = 1; j < ny; j++)
    {
      size_t p = (size_t)i*(ny+1) + j;
      double x = inner.points[3*p] + 0.8*(2.0/nx)*next()*(seed == 0 ? 0.0 : 1.0);
      double y = inner.points[3*p + 1] + 0.8*(3.0/ny)*next()*(seed == 0 ? 0.0 : 1.0);
      inner.points[3*p] = x;
      inner.points[3*p + 1] = y;
      inner.points[3*p + 2] = 2.0*x*x;
    }
  }
}

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
  if (r.numPockets > 0 || r.numShreds > 0)
  {
    printf("  pockets and shreds: %lld components cut off from the rims, %lld pieces dropped with them; %lld open shreds, %lld pieces\n", r.numPockets, r.numPocketPieces, r.numShreds, r.numShredPieces);
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
  // All twelve seeds of test 12: the cleanup runs on every sound envelope
  // the caller would hand it.
  for (int variant = 0; variant < 12; variant++)
  {
    Surface inner;
    AddJitteredValley(inner, 40, 30, variant);
    Surface s = inner;
    const double t = 0.5;
    Extrude(s, t);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    SheetNormals(inner, s);
    Envelope e;
    Report r;
    char name[120];
    snprintf(name, sizeof(name), "9.%d sliver cleanup on a %s valley (seed %d)", variant, variant == 0 ? "regular" : "jittered", variant);
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
    std::vector<double> creaseBefore(e.points);
    std::vector<int> kindBefore(e.pointKind);
    // The crease edges going in, as segments: both ends crease points and the
    // two pieces on the edge from different sheet triangles.
    std::vector<std::pair<ll, ll> > creaseIdsBefore;
    auto creaseEdges = [&](const Envelope &env, std::vector<std::array<double, 6> > &segs, std::vector<std::pair<ll, ll> > *ids)
    {
      if (ids) ids->clear();
      std::map<std::pair<ll, ll>, std::vector<size_t> > on;
      for (size_t i = 0; i + 2 < env.triangles.size(); i += 3)
        for (int j = 0; j < 3; j++)
        {
          ll a = env.triangles[i + j], b = env.triangles[i + (j+1)%3];
          on[std::make_pair(std::min(a, b), std::max(a, b))].push_back(i/3);
        }
      segs.clear();
      for (auto &kv : on)
      {
        ll a = kv.first.first, b = kv.first.second;
        // The envelope marks its crease edges; take the edge if any piece on
        // it says so.
        bool isCrease = false;
        for (size_t k = 0; k < kv.second.size(); k++)
        {
          size_t tri = kv.second[k];
          for (int j = 0; j < 3; j++)
          {
            ll p = env.triangles[3*tri + j], q = env.triangles[3*tri + (j+1)%3];
            if (((p == a && q == b) || (p == b && q == a)) && env.creaseEdge[3*tri + j]) isCrease = true;
          }
        }
        if (!isCrease) continue;
        std::array<double, 6> s = {{env.points[3*a], env.points[3*a+1], env.points[3*a+2], env.points[3*b], env.points[3*b+1], env.points[3*b+2]}};
        segs.push_back(s);
        if (ids) ids->push_back(std::make_pair(a, b));
      }
    };
    std::vector<std::array<double, 6> > creaseSegsBefore;
    creaseEdges(e, creaseSegsBefore, &creaseIdsBefore);
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
        bool crease = false;
        for (int k = 0; k < 2; k++)
          for (int j = 0; j < 3; j++)
          {
            ll p = env.triangles[3*it->second[k] + j], q = env.triangles[3*it->second[k] + (j+1)%3];
            if (((p == it->first.first && q == it->first.second) || (p == it->first.second && q == it->first.first)) && env.creaseEdge[3*it->second[k] + j]) crease = true;
          }
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
    // The wall over the whole top, ends included, where the rims and the
    // block's ridges are: the cleanup may not thin it anywhere.
    auto wallOverTop = [&](const Envelope &env)
    {
      double minDistance = 1e300;
      for (size_t p = 0; p < inner.points.size()/3; p++)
      {
        const double *q = &inner.points[3*p];
        if (std::abs(q[0]) > 0.8 || q[2] < 2.0*q[0]*q[0] - 1e-9)
        {
          continue;
        }
        minDistance = std::min(minDistance, DistanceToMesh(q, env));
      }
      return minDistance;
    };
    // Edges on two pieces whose crease marks disagree.
    auto flagMismatches = [&](const Envelope &env)
    {
      std::map<std::pair<ll, ll>, std::vector<std::pair<size_t, int> > > on;
      for (size_t i = 0; i + 2 < env.triangles.size(); i += 3)
        for (int j = 0; j < 3; j++)
        {
          ll a = env.triangles[i + j], b = env.triangles[i + (j+1)%3];
          on[std::make_pair(std::min(a, b), std::max(a, b))].push_back(std::make_pair(i/3, j));
        }
      ll mismatches = 0;
      for (auto &kv : on)
      {
        if (kv.second.size() != 2) continue;
        if (env.creaseEdge[3*kv.second[0].first + kv.second[0].second] != env.creaseEdge[3*kv.second[1].first + kv.second[1].second]) mismatches++;
      }
      return mismatches;
    };
    double wallBefore = wallOverTop(e);
    ll mismatchesBefore = flagMismatches(e);
    // The wall over the floor, away from the block's ends, whose ridges the
    // jitter reaches.
    auto wallOverFloor = [&](const Envelope &env, int &numSampled)
    {
      double minDistance = 1e300;
      numSampled = 0;
      for (size_t p = 0; p < inner.points.size()/3; p++)
      {
        const double *q = &inner.points[3*p];
        if (std::abs(q[0]) > 0.8 || q[1] < 0.5 || q[1] > 2.5 || q[2] < 2.0*q[0]*q[0] - 1e-9)
        {
          continue;
        }
        minDistance = std::min(minDistance, DistanceToMesh(q, env));
        numSampled++;
      }
      return minDistance;
    };
    int numSampledBefore = 0;
    double floorBefore = wallOverFloor(e, numSampledBefore);
    std::vector<ll> trianglesBefore(e.triangles);
    CleanReport c;
    std::vector<unsigned char> noFixed;
    int rc = svenvelope::CleanEnvelopeSlivers(e, noFixed, 10.0, c);
    printf("  slivers %lld -> %lld (worst %.1f -> %.1f), %lld collapsed, %lld snapped, %lld flipped, %lld points removed, %d passes, %.3f s\n",
        c.numSliversBefore, c.numSliversAfter, c.worstBefore, c.worstAfter, c.numCollapsed, c.numSnapped, c.numFlipped, c.numRemoved, c.numPasses, c.seconds);
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
    // The creases stay where they were, to within a tenth of an edge: every
    // crease point going in that is still used has not moved by more than
    // that, and the middle of every crease edge going in lies that close to
    // a crease edge coming out.
    {
      std::vector<unsigned char> used(e.points.size()/3, 0);
      for (size_t i = 0; i < e.triangles.size(); i++) used[(size_t)e.triangles[i]] = 1;
      // The scale a move is judged against is the surface's own edge length,
      // not the crease edge's: the crumple leaves crease edges far shorter
      // than the mesh, and merging those is the point. It is the mean edge
      // of the surface going in, over every edge of every triangle, which
      // is the scale the cleanup's own budget is set by.
      double meanEdge = 0.0;
      {
        ll numEdges = 0;
        for (size_t i = 0; i + 2 < trianglesBefore.size(); i += 3)
          for (int j = 0; j < 3; j++)
          {
            const double *pa = &creaseBefore[3*trianglesBefore[i + j]], *pb = &creaseBefore[3*trianglesBefore[i + (j+1)%3]];
            meanEdge += std::sqrt((pa[0]-pb[0])*(pa[0]-pb[0]) + (pa[1]-pb[1])*(pa[1]-pb[1]) + (pa[2]-pb[2])*(pa[2]-pb[2]));
            numEdges++;
          }
        meanEdge /= std::max<ll>(numEdges, 1);
      }
      double meanCrease = meanEdge;
      double maxMove = 0.0;
      for (size_t p = 0; p < kindBefore.size(); p++)
      {
        if (kindBefore[p] == 0 || !used[p]) continue;
        double d = std::sqrt((e.points[3*p]-creaseBefore[3*p])*(e.points[3*p]-creaseBefore[3*p]) +
            (e.points[3*p+1]-creaseBefore[3*p+1])*(e.points[3*p+1]-creaseBefore[3*p+1]) +
            (e.points[3*p+2]-creaseBefore[3*p+2])*(e.points[3*p+2]-creaseBefore[3*p+2]));
        maxMove = std::max(maxMove, d);
      }
      snprintf(what, sizeof(what), "crease points moved at most %.4g (a tenth of the mean edge is %.4g)", maxMove, 0.1*meanCrease);
      Check(maxMove <= 0.1*meanCrease, what);
      std::vector<std::array<double, 6> > creaseSegsAfter;
      creaseEdges(e, creaseSegsAfter, nullptr);
      double worstOff = 0.0;
      int numOff = 0;
      // A crease edge one of whose ends the cleanup took away was collapsed
      // within its allowance, which the point check above bounds; the edges
      // whose ends both remain have to still be creases, where they were.
      int numSkipped = 0;
      for (size_t i = 0; i < creaseSegsBefore.size(); i++)
      {
        if (!used[(size_t)creaseIdsBefore[i].first] || !used[(size_t)creaseIdsBefore[i].second])
        {
          numSkipped++;
          continue;
        }
        const std::array<double, 6> &s = creaseSegsBefore[i];
        double mid[3] = {0.5*(s[0]+s[3]), 0.5*(s[1]+s[4]), 0.5*(s[2]+s[5])};
        double len = std::sqrt((s[3]-s[0])*(s[3]-s[0]) + (s[4]-s[1])*(s[4]-s[1]) + (s[5]-s[2])*(s[5]-s[2]));
        double best = 1e300;
        for (size_t k = 0; k < creaseSegsAfter.size(); k++)
        {
          const std::array<double, 6> &q = creaseSegsAfter[k];
          double ab[3] = {q[3]-q[0], q[4]-q[1], q[5]-q[2]}, ap[3] = {mid[0]-q[0], mid[1]-q[1], mid[2]-q[2]};
          double L2 = ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2];
          double u = (L2 > 0.0) ? std::max(0.0, std::min(1.0, (ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2])/L2)) : 0.0;
          double d2 = 0.0;
          for (int m = 0; m < 3; m++) { double dd = mid[m] - (q[m] + u*ab[m]); d2 += dd*dd; }
          best = std::min(best, std::sqrt(d2));
        }
        (void)len;
        worstOff = std::max(worstOff, best);
        if (best > 0.1*meanEdge) numOff++;
      }
      snprintf(what, sizeof(what), "crease edges stay in place: %d of %zu midpoints lie more than a tenth of the mean edge (%.4g) off the new creases (worst %.4g; %d edges collapsed away)", numOff, creaseSegsBefore.size(), 0.1*meanEdge, worstOff, numSkipped);
      Check(numOff == 0, what);
    }
    size_t numReshaped = 0, numWhole = 0;
    for (size_t i = 0; i < e.whole.size(); i++)
    {
      if (e.whole[i] == 2) numReshaped++;
      if (e.whole[i] == 1) numWhole++;
    }
    Check(e.source.size() == e.triangles.size()/3 && e.whole.size() == e.source.size() && e.creaseEdge.size() == e.triangles.size(), "source, whole and creaseEdge follow the triangles");
    Check(c.numCollapsed + c.numFlipped == 0 || numReshaped > 0, "reshaped pieces are marked");
    double vAfter = EnclosedVolume(e);
    snprintf(what, sizeof(what), "volume %.4f -> %.4f", vBefore, vAfter);
    Check(std::abs(vAfter - vBefore) < 0.01*vBefore, what);
    // The wall over the floor: the envelope's own, near the thickness, and
    // not thinned by the cleanup.
    int numSampled = 0;
    double minDistance = wallOverFloor(e, numSampled);
    snprintf(what, sizeof(what), "wall over the valley: %d inner points stand at least %.4f from the envelope before the cleanup and %.4f after (t = %.2f)", numSampled, floorBefore, minDistance, t);
    Check(numSampled > 0 && numSampled == numSampledBefore && floorBefore > 0.85*t && minDistance >= 0.95*floorBefore, what);
    double wallAfter = wallOverTop(e);
    snprintf(what, sizeof(what), "wall over the whole top %.4f -> %.4f (not thinned by more than a twentieth)", wallBefore, wallAfter);
    Check(wallAfter >= 0.95*wallBefore, what);
    ll mismatchesAfter = flagMismatches(e);
    snprintf(what, sizeof(what), "crease marks agree on both pieces of every edge: %lld -> %lld mismatches", mismatchesBefore, mismatchesAfter);
    Check(mismatchesBefore == 0 && mismatchesAfter == 0, what);
    (void)numWhole;
  }

  // 10. A closed surface of its own, however small or coarse, is envelope
  // in full: a unit tetrahedron at the origin and one far from it, whose
  // pieces would have read as a slit under a volume-over-area test (a
  // tetrahedron spreads to 0.07 of its edge) or as a solid under a
  // signed volume taken about the origin.
  for (int where = 0; where < 2; where++)
  {
    Surface s;
    double o[3] = {0.0, 0.0, 0.0};
    if (where == 1) { o[0] = 100.0; o[1] = 50.0; o[2] = -30.0; }
    ll a = AddPoint(s, o[0], o[1], o[2]);
    ll b = AddPoint(s, o[0] + 1.0, o[1], o[2]);
    ll c = AddPoint(s, o[0] + 0.5, o[1] + std::sqrt(3.0)/2.0, o[2]);
    ll d = AddPoint(s, o[0] + 0.5, o[1] + std::sqrt(3.0)/6.0, o[2] + std::sqrt(2.0/3.0));
    // Outward wound: seen from outside, each face runs counter-clockwise.
    AddTriangle(s, a, c, b);
    AddTriangle(s, a, b, d);
    AddTriangle(s, b, c, d);
    AddTriangle(s, c, a, d);
    s.numSheetTriangles = 4;
    Envelope e;
    Report r;
    RunAndCheckClosed(where == 0 ? "10.0 a unit tetrahedron at the origin" : "10.1 a unit tetrahedron far from the origin", s, e, r, 0);
    Check(r.numPockets == 0 && r.numPiecesKept == 4 && r.numWholeKept == 4, "all four faces kept, none taken for a pocket");
    double v = EnclosedVolume(e);
    char what[120];
    snprintf(what, sizeof(what), "volume %.5f of %.5f", v, 1.0/(6.0*std::sqrt(2.0)));
    Check(std::abs(v - 1.0/(6.0*std::sqrt(2.0))) < 1e-9, what);
  }

  // 11. A crease loop closed inside one triangle: a small tetrahedron whose
  // apex pokes up through the top of a large slab, so that the loop where
  // the two cross lies inside a single triangle of the slab's top and is a
  // hole in the piece around it. The union's boundary has the slab's top
  // with the loop cut out and the tetrahedron's tip standing in it.
  {
    Surface s;
    // The slab: [-5, 5] x [-5, 5] x [-1, 0], twelve triangles wound outward.
    ll b0 = AddPoint(s, -5, -5, -1), b1 = AddPoint(s, 5, -5, -1), b2 = AddPoint(s, 5, 5, -1), b3 = AddPoint(s, -5, 5, -1);
    ll t0 = AddPoint(s, -5, -5, 0), t1 = AddPoint(s, 5, -5, 0), t2 = AddPoint(s, 5, 5, 0), t3 = AddPoint(s, -5, 5, 0);
    AddTriangle(s, t0, t1, t2); AddTriangle(s, t0, t2, t3);          // top, up
    AddTriangle(s, b0, b2, b1); AddTriangle(s, b0, b3, b2);          // bottom, down
    AddTriangle(s, b0, b1, t1); AddTriangle(s, b0, t1, t0);          // -y
    AddTriangle(s, b1, b2, t2); AddTriangle(s, b1, t2, t1);          // +x
    AddTriangle(s, b2, b3, t3); AddTriangle(s, b2, t3, t2);          // +y
    AddTriangle(s, b3, b0, t0); AddTriangle(s, b3, t0, t3);          // -x
    // The tetrahedron: base at z = -0.5 inside the slab, apex at z = 0.5
    // above it, placed in the top triangle (t0, t1, t2), which is the half
    // x > y of the top.
    ll a0 = AddPoint(s, 1.5, -2.5, -0.5), a1 = AddPoint(s, 2.5, -2.5, -0.5), a2 = AddPoint(s, 2.0, -1.5, -0.5);
    ll ap = AddPoint(s, 2.0, -2.0, 0.5);
    AddTriangle(s, a0, a2, a1);                                        // base, down
    AddTriangle(s, a0, a1, ap); AddTriangle(s, a1, a2, ap); AddTriangle(s, a2, a0, ap);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed("11. a tetrahedron poking through one triangle of a slab", s, e, r, 0);
    Check(r.numPairsCrossing == 3, "the three faces of the tip cross the top");
    // The tip above z = 0 is the tetrahedron scaled by a half about its
    // apex: an eighth of its volume, which is (1/3)(1/2)(1) = 1/6.
    double v = EnclosedVolume(e);
    char what[160];
    snprintf(what, sizeof(what), "volume %.6f of the slab's 100 plus the tip's %.6f", v, 1.0/48.0);
    Check(std::abs(v - (100.0 + 1.0/48.0)) < 1.0e-6, what);
  }

  // 12. The jittered valley over a dozen seeds: the envelope has to come out
  // sound (no arrangement fault, watertight, no crossing) whatever the
  // jitter, which puts crease loops inside triangles, creases through
  // corners and along edges, and triple points everywhere.
  for (int variant = 0; variant < 12; variant++)
  {
    Surface inner;
    AddJitteredValley(inner, 40, 30, variant);
    Surface s = inner;
    Extrude(s, 0.5);
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    SheetNormals(inner, s);
    Envelope e;
    Report r;
    std::string error;
    int rc = svenvelope::BuildOuterEnvelope(s, e, r, error);
    ll boundary = 0, bad = 0;
    std::vector<unsigned char> crossing;
    double at[3];
    ll numCrossing = (rc == 0) ? svenvelope::CountCrossingTriangles(e.points, e.triangles, crossing, at) : -1;
    char what[200];
    snprintf(what, sizeof(what), "12.%d jittered valley seed %d: built %d, faults %lld, watertight %d, crossings %lld, pockets %lld", variant, variant, rc == 0, r.numArrangementFaults, rc == 0 && Watertight(e, 0, boundary, bad), numCrossing, r.numPockets);
    Check(rc == 0 && r.numArrangementFaults == 0 && Watertight(e, 0, boundary, bad) && numCrossing == 0, what);
  }

  // 13. Crease loops closed inside one triangle, in the arrangements the
  // hole bridging has to get right: two islands in one face, where the
  // right one, a tall thin loop, screens the left one's rightmost vertex
  // from every outer vertex to its right, so that the bridge of the left
  // island has to go to a vertex of the right island already joined in;
  // an island inside another's loop (the small tip pokes through inside
  // the big loop, whose inside is a piece of its own with a hole); and an
  // island whose loop touches the triangle's edge at a vertex, exactly and
  // within a hair. Either order of the islands in the input.
  for (int variant = 0; variant < 6; variant++)
  {
    Surface s;
    AddSlab(s);
    const char *name = "";
    double expected = 100.0;
    if (variant < 2)
    {
      name = (variant == 0) ? "13.0 two islands in one face, left first" : "13.1 two islands in one face, right first";
      const double l0[2] = {1.5, -3.5}, l1[2] = {2.5, -3.5}, l2[2] = {2.0, -2.5}, lap[3] = {2.0, -3.0, 0.5};
      const double r0[2] = {3.4, -4.89}, r1[2] = {3.6, -4.89}, r2[2] = {3.5, 3.975}, rap[3] = {3.5, -0.75, 3.0};
      if (variant == 0) { AddTetrahedron(s, l0, l1, l2, -0.5, lap); AddTetrahedron(s, r0, r1, r2, -0.5, rap); }
      else { AddTetrahedron(s, r0, r1, r2, -0.5, rap); AddTetrahedron(s, l0, l1, l2, -0.5, lap); }
      // The tips above z = 0: the left one as in test 11, the right one the
      // tetrahedron scaled by 3/3.5 about its apex.
      double baseArea = 0.5*0.2*(3.975 + 4.89);
      expected += 1.0/48.0 + (1.0/3.0)*baseArea*3.5*std::pow(3.0/3.5, 3);
    }
    else if (variant < 4)
    {
      name = (variant == 2) ? "13.2 an island inside another's loop, big first" : "13.3 an island inside another's loop, small first";
      const double A0[2] = {0.0, -4.0}, A1[2] = {4.0, -4.0}, A2[2] = {2.0, -1.0}, Aap[3] = {2.0, -3.0, 1.0};
      const double B0[2] = {1.7, -3.3}, B1[2] = {2.3, -3.3}, B2[2] = {2.0, -2.8}, Bap[3] = {2.0, -3.1, 0.5};
      if (variant == 2) { AddTetrahedron(s, A0, A1, A2, -0.5, Aap); AddTetrahedron(s, B0, B1, B2, -0.5, Bap); }
      else { AddTetrahedron(s, B0, B1, B2, -0.5, Bap); AddTetrahedron(s, A0, A1, A2, -0.5, Aap); }
      // The big tip above z = 0 is the big tetrahedron scaled by 2/3 about
      // its apex; the small one lies wholly inside it.
      expected += (1.0/3.0)*6.0*1.5*std::pow(1.0/1.5, 3);
    }
    else
    {
      name = (variant == 4) ? "13.4 an island touching the triangle's edge at a vertex" : "13.5 an island a hair off the triangle's edge";
      double d = (variant == 4) ? 0.0 : 2.0e-7;
      const double A0[2] = {1.5 + d, 1.5}, A1[2] = {2.5, 1.5}, A2[2] = {2.5, 2.0}, Aap[3] = {2.0 + d, 2.0, 0.5};
      AddTetrahedron(s, A0, A1, A2, -0.5, Aap);
      expected += (1.0/3.0)*0.25*1.0/8.0;
    }
    s.numSheetTriangles = (ll)(s.triangles.size()/3);
    Envelope e;
    Report r;
    RunAndCheckClosed(name, s, e, r, 0);
    double v = EnclosedVolume(e);
    char what[160];
    snprintf(what, sizeof(what), "volume %.6f of %.6f", v, expected);
    Check(std::abs(v - expected) < 1.0e-6*expected, what);
  }

  // 14. The jittered valley at the resolution of a real wall: 200 x 150
  // cells, the seeds that once faulted. Near the fold the extrusion folds
  // back onto itself within a hair, and the arrangement there has crease
  // points a hundred thousand times closer together than the triangles
  // are wide: ear clipping has to take needle ears and bridged holes, the
  // pieces of two layers a hair apart have to come out consistently, and
  // a point on an edge has to be on it.
  {
    const int seeds[3] = {1, 2, 6};
    for (int k = 0; k < 3; k++)
    {
      Surface inner;
      AddJitteredValley(inner, 200, 150, seeds[k]);
      Surface s = inner;
      Extrude(s, 0.5);
      s.numSheetTriangles = (ll)(s.triangles.size()/3);
      SheetNormals(inner, s);
      Envelope e;
      Report r;
      std::string error;
      int rc = svenvelope::BuildOuterEnvelope(s, e, r, error);
      ll boundary = 0, bad = 0;
      std::vector<unsigned char> crossing;
      double at[3];
      ll numCrossing = (rc == 0) ? svenvelope::CountCrossingTriangles(e.points, e.triangles, crossing, at) : -1;
      char what[240];
      snprintf(what, sizeof(what), "14.%d jittered valley 200 x 150 seed %d: built %d, faults %lld, watertight %d, crossings %lld, pockets %lld, shreds %lld", k, seeds[k], rc == 0, r.numArrangementFaults, rc == 0 && Watertight(e, 0, boundary, bad), numCrossing, r.numPockets, r.numShreds);
      Check(rc == 0 && r.numArrangementFaults == 0 && Watertight(e, 0, boundary, bad) && numCrossing == 0, what);
    }
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
