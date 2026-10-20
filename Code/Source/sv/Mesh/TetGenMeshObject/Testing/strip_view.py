#!/usr/bin/env python3
"""The zipper strips of a hybrid shell dump, strip by strip, and a picture of one neighbourhood.

The dump is what test_hybrid.cxx (and the drivers built on it) write for SV_HYBRID_DUMP: length-prefixed
arrays of the shell's points (double, 3 per point), triangles (int64, 3 per triangle), markers (int64 per
triangle: 1 zone interface, 2 outer surface, 100+k level k piece and its strips, 300+k prism walls) and
prismPoint (int64 per point: -1 for a piece point, else the prism mesh's id). A strip triangle has one or
two prism corners. Per strip (a connected set of strip triangles) it prints the fan apexes (vertices on
five or more strip triangles: a corner of the loop zipped to one chain point), the loop and chain point
counts and mean edges, and the position of the worst apex, for the picture.

usage: strip_view.py shell.bin [iface.bin | -] [cx cy cz radius out.svg]
The interface file (points, normals, thickness, triangles with faces, as the drivers read it) only adds
the faces by area; pass - without one. The picture is the outer surface (marker 2) within the radius,
projected on its mean normal: strip triangles blue (two loop corners) or orange (two chain corners),
the loop's edges red. Convert with a headless browser, e.g. chromium --headless --screenshot=out.png out.svg.
Stdlib only (no numpy), from a measurement on the user's model of 2026-10-16 (docs/wall-mesh-purpose.md, section 6).
"""
import sys, struct, array, math
def vec(f, code):
    n = struct.unpack('<q', f.read(8))[0]; a = array.array(code); a.frombytes(f.read(n*8)); return a
with open(sys.argv[1], 'rb') as f:
    pts = vec(f, 'd'); tris = vec(f, 'q'); markers = vec(f, 'q'); prism = vec(f, 'q')
nt = len(tris)//3
def P(v): return pts[3*v:3*v+3]
def dist(a, b): return math.sqrt(sum((a[k]-b[k])**2 for k in range(3)))
# interface faces: per face area and centroid
faces = {}
if len(sys.argv) > 2 and sys.argv[2] != '-':
  with open(sys.argv[2], 'rb') as f:
      npts, ntri = struct.unpack('<qq', f.read(16))
      ip = array.array('d'); ip.frombytes(f.read(npts*24)); f.read(npts*24); f.read(npts*24)
      faces = {}
      for t in range(ntri):
          a, b, c, fc = struct.unpack('<qqqi', f.read(28))
          A = ip[3*a:3*a+3]; B = ip[3*b:3*b+3]; C = ip[3*c:3*c+3]
          u = [B[k]-A[k] for k in range(3)]; v = [C[k]-A[k] for k in range(3)]
          n = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
          ar = 0.5*math.sqrt(sum(x*x for x in n))
          e = faces.setdefault(fc, [0.0, 0.0, 0.0, 0.0, 0])
          e[0] += ar; e[4] += 1
          for k in range(3): e[1+k] += ar*(A[k]+B[k]+C[k])/3
if faces: print("interface faces by area:")
for fc, e in sorted(faces.items(), key=lambda x: -x[1][0])[:8]:
    print("  face %d: area %.1f, %d triangles, centroid (%.1f, %.1f, %.1f)" % (fc, e[0], e[4], e[1]/e[0], e[2]/e[0], e[3]/e[0]))
# strip triangles and their components
strip = []
for t in range(nt):
    T = tris[3*t:3*t+3]; k = sum(1 for v in T if prism[v] >= 0)
    if 0 < k < 3: strip.append(t)
parent = {}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
for t in strip:
    for v in tris[3*t:3*t+3]: parent.setdefault(v, v)
    r = find(tris[3*t])
    for v in tris[3*t+1:3*t+3]:
        rv = find(v)
        if rv != r: parent[rv] = r
comp = {}
for t in strip: comp.setdefault(find(tris[3*t]), []).append(t)
rows = []
for root, ts in comp.items():
    loopPts = set(); chainPts = set(); val = {}; loopE = []; chainE = []; mixE = []
    for t in ts:
        T = tris[3*t:3*t+3]
        for j in range(3):
            v = T[j]; w = T[(j+1)%3]
            (loopPts if prism[v] >= 0 else chainPts).add(v)
            val[v] = val.get(v, 0) + 1
            if v < w:
                d = dist(P(v), P(w))
                if prism[v] >= 0 and prism[w] >= 0: loopE.append(d)
                elif prism[v] < 0 and prism[w] < 0: chainE.append(d)
                else: mixE.append(d)
    c = [0.0, 0.0, 0.0]
    for v in loopPts:
        for k in range(3): c[k] += P(v)[k]/len(loopPts)
    fans = sum(1 for v, n in val.items() if n >= 5)
    fanLoop = sum(1 for v, n in val.items() if n >= 5 and prism[v] >= 0)
    worst = max(val.items(), key=lambda x: x[1])[0]; c = list(P(worst))   # the view centre: the worst fan apex
    rows.append((fans, markers[ts[0]], len(ts), len(loopPts), len(chainPts), sum(loopE)/max(len(loopE),1), sum(chainE)/max(len(chainE),1), sum(mixE)/max(len(mixE),1), max(val.values()), fanLoop, c, root))
rows.sort(key=lambda r: -r[0])
print("strips (by fan apexes): fans marker tris loopPts chainPts loopEdge chainEdge across maxVal fansOnLoop centroid")
for r in rows[:24]:
    print("  %4d  m%-3d %5d %5d %5d  %.3f  %.3f  %.3f  %d  %d  (%.1f, %.1f, %.1f)" % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9], r[10][0], r[10][1], r[10][2]))
if len(sys.argv) > 7:
    cx, cy, cz, rad = map(float, sys.argv[3:7]); out = sys.argv[7]
    C = [cx, cy, cz]
    sel = []; stripSet = set(strip)
    for t in range(nt):
        if markers[t] != 2: continue
        T = tris[3*t:3*t+3]
        if dist(P(T[0]), C) < rad: sel.append(t)
    # frame: mean normal of the selected triangles
    nm = [0.0, 0.0, 0.0]
    for t in sel:
        A, B, Cc = [P(v) for v in tris[3*t:3*t+3]]
        u = [B[k]-A[k] for k in range(3)]; v = [Cc[k]-A[k] for k in range(3)]
        n = [u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1]-u[1]*v[0]]
        for k in range(3): nm[k] += n[k]
    l = math.sqrt(sum(x*x for x in nm)); nm = [x/l for x in nm]
    ax = [1.0, 0.0, 0.0] if abs(nm[0]) < 0.9 else [0.0, 1.0, 0.0]
    e1 = [ax[1]*nm[2]-ax[2]*nm[1], ax[2]*nm[0]-ax[0]*nm[2], ax[0]*nm[1]-ax[1]*nm[0]]; l = math.sqrt(sum(x*x for x in e1)); e1 = [x/l for x in e1]
    e2 = [nm[1]*e1[2]-nm[2]*e1[1], nm[2]*e1[0]-nm[0]*e1[2], nm[0]*e1[1]-nm[1]*e1[0]]
    def proj(p):
        d = [p[k]-C[k] for k in range(3)]
        return (sum(d[k]*e1[k] for k in range(3)), sum(d[k]*e2[k] for k in range(3)))
    W = 1200; sc = W/(2*rad)
    def xy(p):
        x, y = proj(p); return (W/2 + sc*x, W/2 - sc*y)
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" style="background:white">' % (W, W)]
    for t in sel:
        T = tris[3*t:3*t+3]
        pp = [xy(P(v)) for v in T]
        if t in stripSet:
            k = sum(1 for v in T if prism[v] >= 0)
            fill = "#9ecae1" if k == 2 else "#fdd0a2"
        else:
            fill = "#f0f0f0"
        svg.append('<polygon points="%s" fill="%s" stroke="#555" stroke-width="0.6"/>' % (" ".join("%.1f,%.1f" % q for q in pp), fill))
    for t in sel:
        T = tris[3*t:3*t+3]
        for j in range(3):
            v, w = T[j], T[(j+1)%3]
            if prism[v] >= 0 and prism[w] >= 0:
                a, b = xy(P(v)), xy(P(w))
                svg.append('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="#d62728" stroke-width="2"/>' % (a[0], a[1], b[0], b[1]))
    svg.append('</svg>')
    open(out, 'w').write("\n".join(svg))
    print("wrote %s: %d triangles (blue: strip with two loop corners, orange: with two chain corners, red: loop edges)" % (out, len(sel)))
