#!/usr/bin/env python3
"""The creases of a hybrid shell dump's outer surface (marker 2), and a shaded picture of one place.

The dump is what test_hybrid.cxx writes for SV_HYBRID_DUMP (see strip_view.py for the layout). Edges whose two
triangles turn by more than 35 degrees are listed by 6-unit cell with their mean turn: a valley between two
vessel walls at a branch root makes a chain of them. The picture is the outer surface within the radius,
viewed along the given direction (or its mean normal), flat-shaded as a headlight would, strips tinted warm;
convert with a headless browser, e.g. chromium --headless --screenshot=out.png out.svg. Used on 2026-10-20 to
see the teeth along the valleys (docs/wall-mesh-purpose.md, section 6) before and after the crease cut.

usage: valley_view.py shell.bin [cx cy cz radius out.svg [vx vy vz]]
"""
import sys, struct, array, math
def vec(f, code):
    n = struct.unpack('<q', f.read(8))[0]; a = array.array(code); a.frombytes(f.read(n*8)); return a
with open(sys.argv[1], 'rb') as f:
    pts = vec(f, 'd'); tris = vec(f, 'q'); markers = vec(f, 'q'); prism = vec(f, 'q')
nt = len(tris)//3
def P(v): return pts[3*v:3*v+3]
def sub(a, b): return [a[k]-b[k] for k in range(3)]
def cross(a, b): return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]
def dot(a, b): return sum(a[k]*b[k] for k in range(3))
def norm(a):
    l = math.sqrt(dot(a, a)); return [x/l for x in a] if l > 0 else [0.0, 0.0, 0.0]
outer = [t for t in range(nt) if markers[t] == 2]
normal = {}
for t in outer:
    T = tris[3*t:3*t+3]; normal[t] = norm(cross(sub(P(T[1]), P(T[0])), sub(P(T[2]), P(T[0]))))
edgeT = {}
for t in outer:
    T = tris[3*t:3*t+3]
    for j in range(3):
        e = (min(T[j], T[(j+1)%3]), max(T[j], T[(j+1)%3])); edgeT.setdefault(e, []).append(t)
crease = []   # (edge, angle)
for e, ts in edgeT.items():
    if len(ts) != 2: continue
    c = max(-1.0, min(1.0, dot(normal[ts[0]], normal[ts[1]])))
    a = math.degrees(math.acos(c))
    if a > 35.0: crease.append((e, a))
print("outer surface: %d triangles; edges with a dihedral turn over 35 degrees: %d" % (len(outer), len(crease)))
# clusters by proximity (grid cells of 6 units)
cell = {}
for e, a in crease:
    m = [(P(e[0])[k]+P(e[1])[k])/2 for k in range(3)]
    key = tuple(int(math.floor(x/6.0)) for x in m); cell.setdefault(key, []).append((m, a))
rows = sorted(cell.items(), key=lambda kv: -len(kv[1]))
print("crease edges by 6-unit cell (count, mean turn, centre):")
for key, lst in rows[:15]:
    c = [sum(m[k] for m, a in lst)/len(lst) for k in range(3)]
    print("  %4d edges, mean turn %.0f, around (%.1f, %.1f, %.1f)" % (len(lst), sum(a for m, a in lst)/len(lst), *c))
if len(sys.argv) > 6:
    C = list(map(float, sys.argv[2:5])); rad = float(sys.argv[5]); out = sys.argv[6]
    sel = []
    for t in outer:
        T = tris[3*t:3*t+3]; c = [(P(T[0])[k]+P(T[1])[k]+P(T[2])[k])/3 for k in range(3)]
        if math.sqrt(sum((c[k]-C[k])**2 for k in range(3))) < rad: sel.append((t, c))
    if len(sys.argv) > 9: view = norm(list(map(float, sys.argv[7:10])))
    else:
        nm = [0.0, 0.0, 0.0]
        for t, c in sel:
            for k in range(3): nm[k] += normal[t][k]
        view = norm(nm)
    ax = [1.0, 0.0, 0.0] if abs(view[0]) < 0.9 else [0.0, 1.0, 0.0]
    e1 = norm(cross(ax, view)); e2 = cross(view, e1)
    W = 1000; sc = W/(2*rad)
    def xy(p):
        d = sub(p, C); return (W/2 + sc*dot(d, e1), W/2 - sc*dot(d, e2))
    sel.sort(key=lambda tc: dot(tc[1], view))   # far first: the view vector points at the camera
    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" style="background:#303030">' % (W, W)]
    for t, c in sel:
        T = tris[3*t:3*t+3]; pp = [xy(P(v)) for v in T]
        shade = dot(normal[t], view)
        g = int(40 + 215*max(0.0, shade)) if shade > 0 else int(40 + 60*(-shade))
        strip = 0 < sum(1 for v in T if prism[v] >= 0) < 3
        col = "rgb(%d,%d,%d)" % ((g, g, g) if not strip else (g, int(g*0.85), int(g*0.7)))
        svg.append('<polygon points="%s" fill="%s" stroke="#204080" stroke-width="0.5"/>' % (" ".join("%.1f,%.1f" % q for q in pp), col))
    svg.append('</svg>'); open(out, 'w').write("\n".join(svg))
    print("wrote %s: %d triangles, view (%.2f %.2f %.2f); strips tinted warm" % (out, len(sel), *view))
