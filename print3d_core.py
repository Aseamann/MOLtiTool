"""
print3d_core.py  -  geometry engine for print3d.py (PyMOL 3D-print export).

Pure numpy, no PyMOL imports, so it can be tested on its own.

Pipeline:  PyMOL VRML text  ->  parse primitives + triangle sets
           -> tessellate cylinders (as round-ended capsules) and spheres
           -> weld, fix triangle orientation, split into connected pieces
           -> cap open ends  ->  optional boolean union per color (manifold3d)
           -> write OBJ with per-vertex colors, or one STL per color.
"""

import re
import numpy as np


# ---------------------------------------------------------------------------
# VRML parsing (PyMOL's get_vrml() output)
# ---------------------------------------------------------------------------

_FLOAT = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'


def _floats(s):
    return np.array(re.findall(_FLOAT, s), dtype=np.float64)


def _block(text, start):
    """Return text of the {...} block whose '{' is at or after start."""
    i = text.index('{', start)
    depth = 0
    for j in range(i, len(text)):
        c = text[j]
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return text[i + 1:j], j + 1
    return text[i + 1:], len(text)


def _bracket_list(block, key):
    m = re.search(re.escape(key) + r'\s*\[', block)
    if not m:
        return None
    end = block.index(']', m.end())
    return block[m.end():end]


def parse_vrml(text):
    """Return (primitives, meshes).
    primitives: list of dicts {type: 'cylinder'|'sphere', center, axis (unit,
                for cylinder), height, radius, color (rgb 0-1)}
    meshes: list of (va (n,3), ta (m,3), vc (n,3) rgb 0-1, vn (n,3) or None)
    """
    prims, meshes = [], []
    pos = 0
    while True:
        m = re.search(r'\b(Transform|Shape)\s*\{', text[pos:])
        if not m:
            break
        kind = m.group(1)
        body, end = _block(text, pos + m.start())
        if kind == 'Transform':
            p = _parse_transform(body)
            if p:
                prims.append(p)
        else:
            mesh = _parse_shape_mesh(body)
            if mesh:
                meshes.append(mesh)
        pos = end
    return prims, meshes


def _axis_angle_matrix(ax, ang):
    ax = np.asarray(ax, float)
    n = np.linalg.norm(ax)
    if n == 0 or ang == 0:
        return np.eye(3)
    x, y, z = ax / n
    c, s, C = np.cos(ang), np.sin(ang), 1 - np.cos(ang)
    return np.array([[c + x*x*C, x*y*C - z*s, x*z*C + y*s],
                     [y*x*C + z*s, c + y*y*C, y*z*C - x*s],
                     [z*x*C - y*s, z*y*C + x*s, c + z*z*C]])


def _parse_transform(body):
    t = re.search(r'translation\s+(' + _FLOAT + r'\s+' + _FLOAT + r'\s+' + _FLOAT + ')', body)
    r = re.search(r'rotation\s+(' + (r'\s+'.join([_FLOAT] * 4)) + ')', body)
    g = re.search(r'geometry\s+(Cylinder|Sphere)\s*\{([^}]*)\}', body)
    col = re.search(r'diffuseColor\s+(' + _FLOAT + r'\s+' + _FLOAT + r'\s+' + _FLOAT + ')', body)
    if not g:
        return None
    center = _floats(t.group(1)) if t else np.zeros(3)
    color = _floats(col.group(1)) if col else np.ones(3)
    radius = float(re.search(r'radius\s+(' + _FLOAT + ')', g.group(2)).group(1))
    if g.group(1) == 'Sphere':
        return dict(type='sphere', center=center, radius=radius, color=color)
    h = re.search(r'height\s+(' + _FLOAT + ')', g.group(2))
    height = float(h.group(1)) if h else 2.0
    rot = _floats(r.group(1)) if r else np.array([0, 0, 1, 0.0])
    axis = _axis_angle_matrix(rot[:3], rot[3]) @ np.array([0.0, 1.0, 0.0])  # VRML cylinders lie along Y
    return dict(type='cylinder', center=center, axis=axis / np.linalg.norm(axis),
                height=height, radius=radius, color=color)


def _parse_shape_mesh(body):
    if 'IndexedFaceSet' not in body:
        return None
    pts = _bracket_list(body, 'point')
    idx = _bracket_list(body, 'coordIndex')
    if pts is None or idx is None:
        return None
    va = _floats(pts).reshape(-1, 3)
    ci = _floats(idx).astype(np.int64)
    # faces separated by -1; PyMOL writes triangles only, but fan-split polygons
    tris, face = [], []
    for k in ci:
        if k < 0:
            for i in range(1, len(face) - 1):
                tris.append((face[0], face[i], face[i + 1]))
            face = []
        else:
            face.append(k)
    ta = np.array(tris, np.int64).reshape(-1, 3)
    cblock = re.search(r'color\s+Color\s*\{', body)
    if cblock:
        cl = _bracket_list(body[cblock.end():], 'color')
        vc = _floats(cl).reshape(-1, 3)
    else:
        col = re.search(r'diffuseColor\s+(' + _FLOAT + r'\s+' + _FLOAT + r'\s+' + _FLOAT + ')', body)
        vc = np.tile(_floats(col.group(1)) if col else np.ones(3), (len(va), 1))
    vn = None
    nblock = re.search(r'normal\s+Normal\s*\{', body)
    if nblock:
        nl = _bracket_list(body[nblock.end():], 'vector')
        if nl is not None:
            vn = _floats(nl).reshape(-1, 3)
            if len(vn) != len(va):
                vn = None
    if len(vc) != len(va):
        vc = np.tile(vc[0] if len(vc) else np.ones(3), (len(va), 1))
    return va, ta, vc, vn


# ---------------------------------------------------------------------------
# Tessellation of primitives (closed, watertight meshes)
# ---------------------------------------------------------------------------

def _frame(axis):
    z = axis / np.linalg.norm(axis)
    a = np.array([1.0, 0, 0]) if abs(z[0]) < 0.9 else np.array([0, 1.0, 0])
    x = np.cross(z, a)
    x /= np.linalg.norm(x)
    return x, np.cross(z, x), z


def capsule_mesh(p0, p1, radius, segments=24, cap_rings=6, round_ends=True):
    """Closed mesh of a cylinder from p0 to p1. Round ends (hemispheres) match
    how PyMOL draws sticks; flat ends otherwise."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = np.linalg.norm(d)
    if L < 1e-9:
        return sphere_mesh(p0, radius, segments)
    x, y, z = _frame(d)
    ang = np.linspace(0, 2 * np.pi, segments, endpoint=False)
    circ = np.outer(np.cos(ang), x) + np.outer(np.sin(ang), y)       # (s,3)
    rings = []                                                        # list of (center, ring radius)
    if round_ends:
        for k in range(cap_rings, 0, -1):                             # bottom hemisphere
            th = (np.pi / 2) * k / cap_rings
            rings.append((p0 - z * radius * np.sin(th), radius * np.cos(th)))
        rings.append((p0, radius))
        rings.append((p1, radius))
        for k in range(1, cap_rings + 1):
            th = (np.pi / 2) * k / cap_rings
            rings.append((p1 + z * radius * np.sin(th), radius * np.cos(th)))
    else:
        rings = [(p0, radius), (p1, radius)]
    verts, tris = [], []
    ring_idx = []
    for c, r in rings:
        if r < 1e-9:                                                  # pole
            ring_idx.append([len(verts)])
            verts.append(c)
        else:
            ring_idx.append(list(range(len(verts), len(verts) + segments)))
            verts.extend(c + r * circ)
    s = segments
    for a, b in zip(ring_idx[:-1], ring_idx[1:]):
        if len(a) == 1:
            tris += [(a[0], b[(i + 1) % s], b[i]) for i in range(s)]
        elif len(b) == 1:
            tris += [(a[i], a[(i + 1) % s], b[0]) for i in range(s)]
        else:
            for i in range(s):
                j = (i + 1) % s
                tris += [(a[i], a[j], b[j]), (a[i], b[j], b[i])]
    va = np.array(verts, np.float64)
    if not round_ends:                                               # flat caps
        c0, c1 = len(va), len(va) + 1
        va = np.vstack([va, p0, p1])
        a, b = ring_idx[0], ring_idx[-1]
        tris += [(c0, a[(i + 1) % s], a[i]) for i in range(s)]
        tris += [(c1, b[i], b[(i + 1) % s]) for i in range(s)]
    return va, np.array(tris, np.int64)


def sphere_mesh(center, radius, segments=24):
    rings = max(4, segments // 2)
    verts = [np.array([0, 0, -1.0])]
    for k in range(1, rings):
        th = -np.pi / 2 + np.pi * k / rings
        for i in range(segments):
            ph = 2 * np.pi * i / segments
            verts.append(np.array([np.cos(th) * np.cos(ph), np.cos(th) * np.sin(ph), np.sin(th)]))
    verts.append(np.array([0, 0, 1.0]))
    s, top = segments, len(verts) - 1
    tris = [(0, 1 + (i + 1) % s, 1 + i) for i in range(s)]
    for k in range(rings - 2):
        a, b = 1 + k * s, 1 + (k + 1) * s
        for i in range(s):
            j = (i + 1) % s
            tris += [(a + i, a + j, b + j), (a + i, b + j, b + i)]
    a = 1 + (rings - 2) * s
    tris += [(a + i, a + (i + 1) % s, top) for i in range(s)]
    return np.asarray(center) + radius * np.array(verts), np.array(tris, np.int64)


# ---------------------------------------------------------------------------
# Mesh cleanup
# ---------------------------------------------------------------------------

def weld_vertices(va, ta, vc, tol=1e-4):
    key = np.round(va / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    nta = inv[ta]
    ok = (nta[:, 0] != nta[:, 1]) & (nta[:, 1] != nta[:, 2]) & (nta[:, 0] != nta[:, 2])
    nta = nta[ok]
    # drop exact duplicate triangles (same vertex set)
    _, keep = np.unique(np.sort(nta, axis=1), axis=0, return_index=True)
    nta = nta[np.sort(keep)]
    return va[first], nta.astype(np.int64), vc[first]


def orient_by_normals(va, ta, vn):
    """PyMOL writes strips with alternating winding; use the exported vertex
    normals to make every triangle face outward."""
    if vn is None:
        return ta
    v0, v1, v2 = va[ta[:, 0]], va[ta[:, 1]], va[ta[:, 2]]
    gn = np.cross(v1 - v0, v2 - v0)
    sn = vn[ta[:, 0]] + vn[ta[:, 1]] + vn[ta[:, 2]]
    flip = (gn * sn).sum(axis=1) < 0
    ta = ta.copy()
    ta[flip] = ta[flip][:, [0, 2, 1]]
    return ta


def connected_components(nv, ta):
    parent = np.arange(nv)

    def find(i):
        root = i
        while parent[root] != root:
            root = parent[root]
        while parent[i] != root:
            parent[i], i = root, parent[i]
        return root
    for a, b, c in ta:
        ra, rb, rc = find(a), find(b), find(c)
        parent[rb] = ra
        parent[find(rc)] = ra
    roots = np.array([find(i) for i in range(nv)])
    labels = np.unique(roots, return_inverse=True)[1].ravel()
    return labels[ta[:, 0]], labels


def split_components(va, ta, vc):
    if len(ta) == 0:
        return []
    tl, vl = connected_components(len(va), ta)
    out = []
    for c in np.unique(tl):
        sel = ta[tl == c]
        used, remap = np.unique(sel.ravel(), return_inverse=True)
        out.append((va[used], remap.reshape(-1, 3), vc[used]))
    return out


def close_holes(va, ta, vc):
    if len(ta) == 0:
        return va, ta, vc, 0
    e = np.concatenate([ta[:, [0, 1]], ta[:, [1, 2]], ta[:, [2, 0]]])
    n = len(va)
    fwd = e[:, 0] * n + e[:, 1]
    rev = e[:, 1] * n + e[:, 0]
    boundary = e[~np.isin(fwd, rev)]
    if len(boundary) == 0:
        return va, ta, vc, 0
    nxt = {}
    for a, b in boundary:
        nxt.setdefault(int(a), int(b))
    new_v, new_c, new_t, visited, nv = [], [], [], set(), n
    for start in list(nxt):
        if start in visited:
            continue
        loop, cur = [], start
        while cur in nxt and cur not in visited:
            visited.add(cur)
            loop.append(cur)
            cur = nxt[cur]
        if cur != start or len(loop) < 3:
            continue
        new_v.append(va[loop].mean(axis=0))
        new_c.append(vc[loop[0]])
        for i, a in enumerate(loop):
            new_t.append((loop[(i + 1) % len(loop)], a, nv))
        nv += 1
    if not new_v:
        return va, ta, vc, 0
    return (np.vstack([va, new_v]), np.vstack([ta, np.array(new_t, np.int64)]),
            np.vstack([vc, new_c]), len(new_v))


def is_watertight(ta):
    if len(ta) == 0:
        return False
    e = np.sort(np.concatenate([ta[:, [0, 1]], ta[:, [1, 2]], ta[:, [2, 0]]]), axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return bool((counts == 2).all())


# ---------------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------------

def to_uint8(rgb):
    return np.clip(np.round(np.asarray(rgb) * 255), 0, 255).astype(np.uint8)


def triangle_colors(ta, vc):
    c0, c1, c2 = vc[ta[:, 0]], vc[ta[:, 1]], vc[ta[:, 2]]
    eq12 = (c1 == c2).all(axis=1)
    out = c0.copy()
    out[eq12] = c1[eq12]
    return out


def single_color_triangles(va, ta, vc):
    tc = triangle_colors(ta, vc)
    _, cidx = np.unique(tc, axis=0, return_inverse=True)
    cidx = cidx.ravel()
    ncol = cidx.max() + 1
    keys = ta * ncol + cidx[:, None]
    ukeys, inv = np.unique(keys.ravel(), return_inverse=True)
    new_ta = inv.reshape(-1, 3)
    new_vc = np.empty((len(ukeys), vc.shape[1]), vc.dtype)
    new_vc[new_ta.ravel()] = np.repeat(tc, 3, axis=0)
    return va[ukeys // ncol], new_ta, new_vc


def apply_color_tolerance(vc, tolerance):
    if tolerance <= 0 or len(vc) == 0:
        return vc
    uniq, inv, counts = np.unique(vc, axis=0, return_inverse=True, return_counts=True)
    reps, mapping = [], np.empty(len(uniq), np.int64)
    for i in np.argsort(-counts):
        c = uniq[i].astype(int)
        for ri in reps:
            if np.abs(uniq[ri].astype(int) - c).max() <= tolerance:
                mapping[i] = ri
                break
        else:
            reps.append(i)
            mapping[i] = i
    return uniq[mapping][inv.ravel()]


def hex_color(c):
    return '%02x%02x%02x' % tuple(int(x) for x in c[:3])


# ---------------------------------------------------------------------------
# Union (optional, needs manifold3d)
# ---------------------------------------------------------------------------

def have_manifold():
    try:
        import manifold3d  # noqa
        return True
    except ImportError:
        return False


def union_same_color(pieces, notes):
    """pieces: list of (va, ta, vc uint8). Single-color pieces of the same
    color are merged into one watertight solid; others pass through."""
    try:
        import manifold3d as mf
    except ImportError:
        notes.append('Merge skipped: manifold3d not installed (use the GUI button '
                     '"Install merge library" or: pip install manifold3d).')
        return pieces
    groups, out, failed = {}, [], 0
    for va, ta, vc in pieces:
        if len(vc) and (vc == vc[0]).all():
            groups.setdefault(tuple(vc[0]), []).append((va, ta))
        else:
            out.append((va, ta, vc))
    for color, plist in groups.items():
        solids, raw = [], []
        for va, ta in plist:
            m = mf.Manifold(mf.Mesh(vert_properties=np.ascontiguousarray(va, np.float32),
                                    tri_verts=np.ascontiguousarray(ta, np.uint32)))
            if m.status() == mf.Error.NoError and not m.is_empty():
                solids.append(m)
            else:
                raw.append((va, ta))
                failed += 1
        if solids:
            u = mf.Manifold.batch_boolean(solids, mf.OpType.Add)
            bodies = [b for b in u.decompose() if b.volume() > 1e-3]
            if bodies:
                u = mf.Manifold.compose(bodies)
            mesh = u.to_mesh()
            uva = np.asarray(mesh.vert_properties)[:, :3].astype(np.float64)
            uta = np.asarray(mesh.tri_verts).astype(np.int64)
            out.append((uva, uta, np.tile(np.array(color, np.uint8), (len(uva), 1))))
        for va, ta in raw:
            out.append((va, ta, np.tile(np.array(color, np.uint8), (len(va), 1))))
    notes.append('Merged sticks, spheres and struts into one solid per color'
                 + (' (%d pieces kept as-is)' % failed if failed else ''))
    return out


# ---------------------------------------------------------------------------
# Whole pipeline
# ---------------------------------------------------------------------------

def vrml_to_pieces(text, segments=24, round_ends=True, close_ends=True, notes=None,
                   info=None):
    """Parse PyMOL VRML and return list of closed pieces (va, ta, vc uint8).
    info (dict, optional) receives 'radii': list of (type, radius) of primitives."""
    notes = [] if notes is None else notes
    prims, meshes = parse_vrml(text)
    if info is not None:
        info['radii'] = [(p['type'], p['radius']) for p in prims]
    pieces = []
    for p in prims:
        col = to_uint8(p['color'])
        if p['type'] == 'sphere':
            va, ta = sphere_mesh(p['center'], p['radius'], segments)
        else:
            half = 0.5 * p['height'] * p['axis']
            va, ta = capsule_mesh(p['center'] - half, p['center'] + half, p['radius'],
                                  segments, max(3, segments // 4), round_ends)
        pieces.append((va, ta, np.tile(col, (len(va), 1))))
    capped = 0
    for va, ta, vc, vn in meshes:
        ta = orient_by_normals(va, ta, vn)
        va, ta, vc = weld_vertices(va, ta, to_uint8(vc))
        for cva, cta, cvc in split_components(va, ta, vc):
            if close_ends:
                cva, cta, cvc, n = close_holes(cva, cta, cvc)
                capped += n
            pieces.append((cva, cta, cvc))
    nprim = len(prims)
    notes.append('Read %d sticks/spheres/struts and %d surface/cartoon pieces'
                 % (nprim, len(pieces) - nprim))
    if capped:
        notes.append('Capped %d open ends in cartoon/surface meshes' % capped)
    return pieces


def combine(pieces):
    if not pieces:
        return np.zeros((0, 3)), np.zeros((0, 3), np.int64), np.zeros((0, 3), np.uint8)
    vas, tas, vcs, off = [], [], [], 0
    for va, ta, vc in pieces:
        vas.append(va)
        tas.append(ta + off)
        vcs.append(vc[:, :3])
        off += len(va)
    return np.vstack(vas), np.vstack(tas), np.vstack(vcs)


def split_by_color(va, ta, vc):
    tc = triangle_colors(ta, vc)
    uniq, inv, counts = np.unique(tc, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    for ci in np.argsort(-counts):
        sel = ta[inv == ci]
        used, remap = np.unique(sel.ravel(), return_inverse=True)
        yield uniq[ci], va[used], remap.reshape(-1, 3)


def write_obj(path, va, ta, vc, comment=''):
    with open(path, 'w') as f:
        f.write('# %s\n' % (comment or 'print3d export'))
        f.write('# vertex colors as "v x y z r g b" (0-1)\n')
        rgb = vc[:, :3].astype(np.float64) / 255.0
        np.savetxt(f, np.hstack([va.astype(np.float64), rgb]),
                   fmt='v %.4f %.4f %.4f %.4f %.4f %.4f')
        np.savetxt(f, ta.astype(np.int64) + 1, fmt='f %d %d %d')


def write_binary_stl(path, va, ta, comment=''):
    v = va[ta].astype(np.float32)
    n = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    ln = np.linalg.norm(n, axis=1)
    ln[ln == 0] = 1
    rec = np.zeros(len(ta), dtype=[('n', '<f4', (3,)), ('v', '<f4', (3, 3)), ('a', '<u2')])
    rec['n'] = (n / ln[:, None]).astype(np.float32)
    rec['v'] = v
    with open(path, 'wb') as f:
        f.write((comment or 'print3d export').encode('ascii', 'replace')[:80].ljust(80, b' '))
        f.write(np.uint32(len(ta)).tobytes())
        f.write(rec.tobytes())


# ===========================================================================
# v2 engine: placement (orientation + scale), filament palette, colored
# solids (manifold3d), printability report, base/stand, connectors, 3MF.
# Colors ride along as per-vertex properties, so booleans keep them.
# ===========================================================================

PLA_DENSITY = 1.24   # g/cm^3

# camera (x right, y up, z toward viewer) -> printer bed (X right, Y back, Z up)
CAMERA_TO_BED = np.array([[1.0, 0, 0], [0, 0, 1.0], [0, -1.0, 0]])


class Placement:
    """out = x @ L + offset   (x in Angstrom model coordinates, out in mm)."""

    def __init__(self, L, offset, scale):
        self.L, self.offset, self.scale = np.asarray(L, float), np.asarray(offset, float), float(scale)

    def apply(self, x):
        return np.asarray(x, float) @ self.L + self.offset

    def inverse(self, y):
        return (np.asarray(y, float) - self.offset) @ np.linalg.inv(self.L)


def view_rotation(view):
    """PyMOL get_view() -> 3x3 model->bed rotation ('print it as seen')."""
    R = np.array(view[:9], float).reshape(3, 3)   # camera = x @ R
    return R @ CAMERA_TO_BED


def make_placement(points, rotation=None, mm_per_A=1.0, size_mm=0.0):
    """Rotate, scale (size_mm = longest side, overrides mm_per_A) and move so
    the model is centered in X/Y and sits on Z = 0."""
    Rm = np.eye(3) if rotation is None else np.asarray(rotation, float)
    p = np.asarray(points, float) @ Rm
    ext = p.max(0) - p.min(0)
    scale = float(size_mm) / max(ext.max(), 1e-9) if size_mm and float(size_mm) > 0 else float(mm_per_A or 1.0)
    p *= scale
    lo, hi = p.min(0), p.max(0)
    offset = -np.array([(lo[0] + hi[0]) / 2, (lo[1] + hi[1]) / 2, lo[2]])
    return Placement(Rm * scale, offset, scale)


# ---------------------------------------------------------------------------
# filament palette
# ---------------------------------------------------------------------------

def _srgb_to_lab(rgb255):
    c = np.asarray(rgb255, float) / 255.0
    c = np.where(c > 0.04045, ((c + 0.055) / 1.055) ** 2.4, c / 12.92)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = c @ M.T / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > 0.008856, np.cbrt(xyz), 7.787 * xyz + 16 / 116)
    return np.stack([116 * f[..., 1] - 16, 500 * (f[..., 0] - f[..., 1]),
                     200 * (f[..., 1] - f[..., 2])], axis=-1)


def snap_to_palette(colors, palette):
    """colors (n,3) uint8, palette (k,3) uint8 -> (snapped colors, slot index)."""
    colors = np.asarray(colors, np.uint8)
    pal = np.asarray(palette, np.uint8)
    uniq, inv = np.unique(colors, axis=0, return_inverse=True)
    d = ((_srgb_to_lab(uniq)[:, None, :] - _srgb_to_lab(pal)[None, :, :]) ** 2).sum(-1)
    slot = d.argmin(1)
    return pal[slot][inv.ravel()], slot[inv.ravel()]


def palette_mapping(colors, palette):
    """Distinct input colors -> (color, slot) pairs, for the report."""
    uniq = np.unique(np.asarray(colors, np.uint8), axis=0)
    _, slots = snap_to_palette(uniq, palette)
    return [(tuple(int(x) for x in c), int(s)) for c, s in zip(uniq, slots)]


# ---------------------------------------------------------------------------
# colored solids
# ---------------------------------------------------------------------------

def colored_solid(va, ta, vc):
    """Manifold with rgb (0-1) carried as vertex properties, or None."""
    import manifold3d as mf
    vp = np.hstack([np.asarray(va, np.float64), np.asarray(vc, np.float64)[:, :3] / 255.0])
    m = mf.Manifold(mf.Mesh(vert_properties=np.ascontiguousarray(vp, np.float32),
                            tri_verts=np.ascontiguousarray(ta, np.uint32)))
    if m.status() != mf.Error.NoError or m.is_empty():
        return None
    return m


def solid_mesh(m):
    """-> (va float64 (n,3), ta int64, vc uint8 (n,3))"""
    me = m.to_mesh()
    vp = np.asarray(me.vert_properties, np.float64)
    ta = np.asarray(me.tri_verts, np.int64)
    if vp.shape[1] >= 6:
        vc = np.clip(np.round(vp[:, 3:6] * 255), 0, 255).astype(np.uint8)
    else:
        vc = np.full((len(vp), 3), 200, np.uint8)
    return vp[:, :3].copy(), ta, vc


def _mesh_with(me, vert_properties):
    import manifold3d as mf
    kw = dict(vert_properties=np.ascontiguousarray(vert_properties, np.float32),
              tri_verts=np.ascontiguousarray(np.asarray(me.tri_verts), np.uint32))
    mfrom, mto = np.asarray(me.merge_from_vert), np.asarray(me.merge_to_vert)
    if len(mfrom):
        kw['merge_from_vert'] = np.ascontiguousarray(mfrom, np.uint32)
        kw['merge_to_vert'] = np.ascontiguousarray(mto, np.uint32)
    m = mf.Manifold(mf.Mesh(**kw))
    return m if m.status() == mf.Error.NoError and not m.is_empty() else None


def tint(m, rgb255):
    """Give a solid a single color (as properties), keeping its topology."""
    me = m.to_mesh()
    vp = np.asarray(me.vert_properties, np.float64)
    col = np.tile(np.asarray(rgb255, float) / 255.0, (len(vp), 1))
    return _mesh_with(me, np.hstack([vp[:, :3], col]))


def main_color(ta, vc):
    tc = triangle_colors(ta, vc)
    u, cnt = np.unique(tc, axis=0, return_counts=True)
    return tuple(int(x) for x in u[np.argmax(cnt)])


class Part:
    """One filament/color region of a printed object."""

    def __init__(self, color, solid=None, mesh=None, name=None):
        self.name = name
        self.color = tuple(int(x) for x in color)
        self.solid = solid            # manifold3d.Manifold (colored) or None
        self.mesh = mesh              # (va, ta, vc) when no solid is available

    def get_mesh(self):
        return solid_mesh(self.solid) if self.solid is not None else self.mesh

    def volume(self):
        return self.solid.volume() if self.solid is not None else 0.0


def build_parts(pieces, notes):
    """pieces (va, ta, vc uint8) in mm -> list of Part, one per main color.
    Each part is the boolean union of its pieces (colors kept per vertex)."""
    import manifold3d as mf
    groups, loose = {}, []
    for va, ta, vc in pieces:
        col = main_color(ta, vc)
        s = colored_solid(va, ta, vc)
        if s is None:
            loose.append((col, (va, ta, vc)))
        else:
            groups.setdefault(col, []).append(s)
    parts = []
    for col, solids in sorted(groups.items(), key=lambda kv: -sum(s.volume() for s in kv[1])):
        u = mf.Manifold.batch_boolean(solids, mf.OpType.Add) if len(solids) > 1 else solids[0]
        bodies = [b for b in u.decompose() if b.volume() > 1e-6 * max(u.volume(), 1e-9)]
        if bodies and len(bodies) != len(u.decompose()):
            u = mf.Manifold.compose(bodies)
        parts.append(Part(col, solid=u))
    for col, mesh in loose:
        parts.append(Part(col, mesh=mesh))
    if loose:
        notes.append('%d pieces could not be made solid and are kept as plain meshes '
                     '(no merging/cutting for them)' % len(loose))
    return parts


def transform_pieces(pieces, placement):
    return [(placement.apply(va), ta, vc) for va, ta, vc in pieces]


# ---------------------------------------------------------------------------
# printability report
# ---------------------------------------------------------------------------

def analyze(parts, near_fn=None, touch_mm=0.2):
    """Return dict: size (mm), volume_mm3, grams, per_part volumes, bodies
    (main + floating: list of (volume, centroid)). near_fn(centroid_mm)->str."""
    import manifold3d as mf
    solids = [p.solid for p in parts if p.solid is not None]
    out = dict(parts=[(p.color, p.volume()) for p in parts])
    allv = [p.get_mesh()[0] for p in parts]
    pts = np.vstack(allv) if allv else np.zeros((1, 3))
    out['size'] = pts.max(0) - pts.min(0)
    if not solids:
        out.update(volume=0.0, grams=0.0, bodies=[])
        return out
    u = mf.Manifold.batch_boolean(solids, mf.OpType.Add) if len(solids) > 1 else solids[0]
    out['volume'] = u.volume()
    out['grams'] = u.volume() / 1000.0 * PLA_DENSITY
    # bodies closer than touch_mm count as attached (they fuse when printed)
    parts_b = u.decompose()
    n = len(parts_b)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    if n > 1:
        boxes = [np.asarray(b.bounding_box(), float) for b in parts_b]
        for i in range(n):
            for j in range(i + 1, n):
                bi, bj = boxes[i], boxes[j]
                if np.any(bi[:3] - touch_mm > bj[3:]) or np.any(bj[:3] - touch_mm > bi[3:]):
                    continue
                if parts_b[i].min_gap(parts_b[j], touch_mm) < touch_mm:
                    parent[find(i)] = find(j)
    clusters = {}
    for i, b in enumerate(parts_b):
        clusters.setdefault(find(i), []).append(b)
    bodies = []
    for bl in clusters.values():
        vol = sum(b.volume() for b in bl)
        bb = np.array([np.asarray(b.bounding_box(), float) for b in bl])
        c = (bb[:, :3].min(0) + bb[:, 3:].max(0)) / 2
        bodies.append((vol, c, near_fn(c) if near_fn else ''))
    bodies.sort(key=lambda t: -t[0])
    out['bodies'] = bodies
    try:
        out['bed_contact'] = u.slice(0.2).area()
    except Exception:
        out['bed_contact'] = None
    return out


# ---------------------------------------------------------------------------
# base: flat cut or pillar stand (Z is up; model sits on Z = 0)
# ---------------------------------------------------------------------------

def _all_solids(parts):
    import manifold3d as mf
    s = [p.solid for p in parts if p.solid is not None]
    if not s:
        return None
    return mf.Manifold.batch_boolean(s, mf.OpType.Add) if len(s) > 1 else s[0]


def _shift(parts, dz):
    for p in parts:
        if p.solid is not None:
            p.solid = p.solid.translate((0, 0, dz))
        elif p.mesh is not None:
            va, ta, vc = p.mesh
            p.mesh = (va + np.array([0, 0, dz]), ta, vc)


def flat_cut(parts, depth_mm, notes):
    """Slice off the bottom depth_mm so the model has a flat face on the bed."""
    allm = _all_solids(parts)
    if allm is None:
        notes.append('flat cut skipped (needs solid parts)')
        return parts
    zmin = np.asarray(allm.bounding_box())[2]
    zc = zmin + float(depth_mm)
    for p in parts:
        if p.solid is not None:
            p.solid = p.solid.trim_by_plane((0, 0, 1), zc)
    parts = [p for p in parts if p.solid is None or not p.solid.is_empty()]
    _shift(parts, -zc)
    area = _all_solids(parts).slice(0.1).area()
    notes.append('flat cut %.1f mm off the bottom: %.0f mm2 touches the bed%s'
                 % (depth_mm, area, ' (small, consider "stand")' if area < 80 else ''))
    return parts


def pillar_stand(parts, n_pillars=3, pillar_d=3.0, plate_t=2.0, gap=4.0, margin=4.0,
                 notes=None, host_parts=None):
    """Add a plate under the model and pillars up to its lowest points.
    Returns the stand as a solid (in final coordinates) and shifts parts."""
    import manifold3d as mf
    notes = [] if notes is None else notes
    allm = _all_solids(parts)
    host = _all_solids(host_parts) if host_parts else allm
    if allm is None or host is None:
        notes.append('stand skipped (needs solid parts)')
        return None
    va, _, _ = solid_mesh(host)
    zmin, zmax = va[:, 2].min(), va[:, 2].max()
    low = va[va[:, 2] <= zmin + 0.3 * (zmax - zmin)]
    # farthest-point sampling in XY, weighted toward low points
    picks = [low[np.argmin(low[:, 2])]]
    for _ in range(int(n_pillars) - 1):
        d = np.min(np.linalg.norm(low[:, None, :2] - np.array(picks)[None, :, :2], axis=2), axis=1)
        score = d - 0.5 * (low[:, 2] - zmin)
        picks.append(low[np.argmax(score)])
    r = float(pillar_d) / 2
    top_embed = 0.8
    pillars = []
    for p in picks:
        h = (p[2] - zmin) + gap + top_embed
        c = mf.Manifold.cylinder(h, r, r, 24).translate((p[0], p[1], zmin - gap))
        pillars.append(c)
    plate_cs = allm.project().hull().offset(float(margin), mf.JoinType.Round)
    plate = mf.Manifold.extrude(plate_cs, float(plate_t)).translate((0, 0, zmin - gap - plate_t))
    stand = mf.Manifold.batch_boolean([plate] + pillars, mf.OpType.Add)
    # pillars must not poke through the model: keep only what is below it or inside the host
    stand = stand - (allm - host) if allm is not host else stand
    shift = -(zmin - gap - plate_t)
    _shift(parts, shift)
    stand = stand.translate((0, 0, shift))
    notes.append('stand: %.0f mm plate, %d pillars %.1f mm thick' % (plate_t, len(picks), pillar_d))
    return stand


# ---------------------------------------------------------------------------
# snap-fit connectors between two printed objects
# ---------------------------------------------------------------------------

def _cyl_between(p0, p1, r, seg=40):
    import manifold3d as mf
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    d = p1 - p0
    L = float(np.linalg.norm(d))
    c = mf.Manifold.cylinder(L, r, r, seg)          # along +Z from origin
    z = d / L
    x = np.cross(z, [1, 0, 0] if abs(z[0]) < 0.9 else [0, 1, 0])
    x /= np.linalg.norm(x)
    y = np.cross(z, x)
    M = np.column_stack([x, y, z, p0])              # 3x4
    return c.transform(M)


def dilate(m, dist):
    """Approximate outward offset by moving vertices along their normals
    (normals computed on the welded topology, so seams stay closed)."""
    if dist <= 0:
        return m
    me = m.to_mesh()
    vp = np.asarray(me.vert_properties, np.float64).copy()
    ta = np.asarray(me.tri_verts, np.int64)
    canon = np.arange(len(vp))
    mfrom, mto = np.asarray(me.merge_from_vert, np.int64), np.asarray(me.merge_to_vert, np.int64)
    if len(mfrom):
        canon[mfrom] = mto
    tc = canon[ta]
    va = vp[:, :3]
    n = np.zeros_like(va)
    fn = np.cross(va[tc[:, 1]] - va[tc[:, 0]], va[tc[:, 2]] - va[tc[:, 0]])
    for k in range(3):
        np.add.at(n, tc[:, k], fn)
    n = n[canon]
    ln = np.linalg.norm(n, axis=1)
    ln[ln == 0] = 1
    vp[:, :3] = va + n / ln[:, None] * dist
    out = _mesh_with(me, vp)
    return out if out is not None else m


def _first_hit(m, a, b, entering=True):
    hits = m.ray_cast(tuple(map(float, a)), tuple(map(float, b)))
    if not hits:
        return None
    seg = np.asarray(b, float) - np.asarray(a, float)
    for h in hits:
        inward = float(np.dot(np.asarray(h.normal), seg)) < 0
        if inward == entering:
            return np.asarray(a, float) + seg * h.distance
    return None


def snapfit(host_parts, guest_parts, host_atoms=None, guest_atoms=None, kind='peg', diameter=4.0,
            length=6.0, clearance=0.2, wall=1.6, magnet_d=6.0, magnet_h=3.0,
            contact_mm=None, notes=None, debug=None):
    """Modify parts in place: guest gets a seat (cut to the host's shape plus
    clearance) and either a peg (host gets boss + hole) or magnet pockets in
    both. host/guest_atoms are atom positions in mm."""
    import manifold3d as mf
    notes = [] if notes is None else notes
    host = _all_solids(host_parts)
    guest = _all_solids(guest_parts)
    if host is None or guest is None:
        notes.append('connector skipped (needs solid parts)')
        return
    # contact from the displayed shapes (atoms may be hidden under a cartoon)
    hv = solid_mesh(host)[0]
    gv = solid_mesh(guest)[0]
    rng = np.random.default_rng(0)
    if len(hv) > 20000:
        hv = hv[rng.choice(len(hv), 20000, replace=False)]
    if len(gv) > 6000:
        gv = gv[rng.choice(len(gv), 6000, replace=False)]
    best = (1e30, 0, 0)
    for k in range(0, len(gv), 400):
        d2 = ((gv[k:k + 400, None, :] - hv[None, :, :]) ** 2).sum(-1)
        i = np.unravel_index(np.argmin(d2), d2.shape)
        if d2[i] < best[0]:
            best = (d2[i], k + i[0], i[1])
    g, h = gv[best[1]], hv[best[2]]
    C = (g + h) / 2
    R = max(4.0 * diameter, 8.0)
    h_loc = hv[np.linalg.norm(hv - C, axis=1) < R].mean(0)
    g_cen = gv[np.linalg.norm(gv - C, axis=1) < 2 * R].mean(0)
    u = g_cen - h_loc
    if np.linalg.norm(u) < 1e-6:
        u = g - h
    u /= np.linalg.norm(u)
    if np.sqrt(best[0]) > 2 * diameter:
        notes.append('connector: parts are %.1f mm apart at their closest point' % np.sqrt(best[0]))

    # 1. seat: the guest yields to the host (+ clearance)
    host_d = dilate(host, clearance)
    for p in guest_parts:
        if p.solid is not None:
            p.solid = p.solid - host_d
    guest = _all_solids(guest_parts)
    if guest is None or guest.is_empty():
        notes.append('connector skipped: nothing left of the guest after the seat cut')
        return
    far = 4 * (np.ptp(np.vstack([hv, gv]), axis=0).max() + 50)
    start = C - u * far                      # outside everything, on the host side
    Pg = _first_hit(guest, start, C + u * far, entering=True)
    if Pg is None:
        notes.append('connector skipped: could not find the guest surface')
        return
    # last point where the ray leaves the host before reaching the guest
    Ph = None
    seg = (C + u * far) - start
    for hit in host.ray_cast(tuple(map(float, start)), tuple(map(float, Pg))):
        if float(np.dot(np.asarray(hit.normal), seg)) > 0:
            Ph = start + (Pg - start) * hit.distance
    gap = float(np.linalg.norm(Pg - Ph)) if Ph is not None else 0.0
    if debug is not None:
        debug.update(C=C, u=u, Pg=Pg, Ph=Ph, gap=gap, g=g, h=h)
    color_h = max(host_parts, key=lambda p: p.volume()).color
    color_g = max(guest_parts, key=lambda p: p.volume()).color

    def add_to(parts, solid, color):
        if solid is None or solid.is_empty():
            notes.append('connector: nothing to add (shape was empty)')
            return
        target = max(parts, key=lambda p: p.volume())
        t = tint(solid, color)
        if t is not None:
            target.solid = target.solid + t

    def cut_all(parts, solid):
        for p in parts:
            if p.solid is not None:
                p.solid = p.solid - solid

    if kind == 'peg':
        r = diameter / 2.0
        root = Pg + u * (length / 2.0)              # half the peg inside the guest
        tip = Pg - u * (gap + length / 2.0)         # half reaches into the host
        peg = _cyl_between(root, tip, r)
        boss = _cyl_between(Pg, tip - u * (clearance + wall), r + clearance + wall)
        boss = boss - dilate(guest, clearance)
        add_to(host_parts, boss, color_h)
        cut_all(host_parts, _cyl_between(Pg + u * length, tip - u * clearance, r + clearance))
        add_to(guest_parts, peg, color_g)
        half = _cyl_between(Pg, root, r)
        anchored = (half ^ guest).volume() / max(half.volume(), 1e-9)
        if anchored < 0.8:
            notes.append('warning: only %.0f%% of the peg root is inside the guest; '
                         'use a shorter/thinner peg or a larger scale' % (100 * anchored))
        notes.append('peg %.1f x %.1f mm (+%.2f mm clearance) at %.1f mm gap'
                     % (diameter, length, clearance, gap))
    else:
        r = magnet_d / 2.0 + clearance
        depth = magnet_h + clearance
        boss = _cyl_between(Pg, Pg - u * (gap + depth + wall), r + wall)
        boss = boss - dilate(guest, clearance)
        add_to(host_parts, boss, color_h)
        pocket_h = _cyl_between(Pg + u * 0.5, Pg - u * (gap + depth), r)
        pocket_g = _cyl_between(Pg - u * 0.5, Pg + u * depth, r)
        cut_all(host_parts, pocket_h)
        shell = _cyl_between(Pg, Pg + u * (depth + wall * 0.75), r + wall * 0.75)
        frac = (shell ^ guest).volume() / max(shell.volume(), 1e-9)
        cut_all(guest_parts, pocket_g)
        notes.append('magnet pockets %.1f x %.1f mm in both parts (+%.2f mm clearance)'
                     % (magnet_d, magnet_h, clearance))
        if frac < 0.85:
            notes.append('warning: the guest is thin around its magnet pocket (%.0f%% '
                         'filled); use a smaller magnet or a larger scale' % (100 * frac))
    return u, Pg


# ---------------------------------------------------------------------------
# 3MF (OrcaSlicer / Anycubic Slicer Next / Bambu Studio)
# ---------------------------------------------------------------------------

def paint_code(slot):
    """Orca/Bambu per-triangle filament paint code for filament number slot."""
    s = int(slot)
    if s <= 0:
        return ''
    if s < 3:
        return '%X' % (s << 2)
    if s < 18:
        return '%XC' % (s - 3)
    return '%XFC' % (s - 18)


def _xml_attr(s):
    return (str(s).replace('&', '&amp;').replace('"', '&quot;')
            .replace('<', '&lt;').replace('>', '&gt;'))


def write_3mf(path, object_name, parts3mf):
    """parts3mf: list of dict(name, va, ta, extruder, paint=None or per-triangle
    filament numbers). One object, one part per entry, each part assigned to
    its filament; painted triangles override it."""
    import zipfile
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<model unit="millimeter" xml:lang="en-US" '
             'xmlns="http://schemas.microsoft.com/3dmanufacturing/core/2015/02">',
             ' <metadata name="Title">%s</metadata>' % _xml_attr(object_name),
             ' <metadata name="Application">print3d (PyMOL)</metadata>',
             ' <resources>']
    ids = []
    for i, p in enumerate(parts3mf, start=1):
        ids.append(i)
        va, ta = np.asarray(p['va'], float), np.asarray(p['ta'], np.int64)
        lines.append('  <object id="%d" type="model" name="%s"><mesh><vertices>'
                     % (i, _xml_attr(p['name'])))
        lines.append('\n'.join('<vertex x="%.4f" y="%.4f" z="%.4f"/>' % tuple(v) for v in va))
        lines.append('</vertices><triangles>')
        paint = p.get('paint')
        if paint is None:
            lines.append('\n'.join('<triangle v1="%d" v2="%d" v3="%d"/>' % tuple(t) for t in ta))
        else:
            rows = []
            for t, s in zip(ta, paint):
                code = paint_code(s) if int(s) != int(p['extruder']) else ''
                rows.append('<triangle v1="%d" v2="%d" v3="%d"%s/>' % (
                    t[0], t[1], t[2], ' paint_color="%s"' % code if code else ''))
            lines.append('\n'.join(rows))
        lines.append('</triangles></mesh></object>')
    asm = len(parts3mf) + 1
    lines.append('  <object id="%d" type="model" name="%s"><components>' % (asm, _xml_attr(object_name)))
    lines += ['   <component objectid="%d"/>' % i for i in ids]
    lines += ['  </components></object>', ' </resources>',
              ' <build><item objectid="%d"/></build>' % asm, '</model>']
    cfg = ['<?xml version="1.0" encoding="UTF-8"?>', '<config>',
           '  <object id="%d">' % asm,
           '    <metadata key="name" value="%s"/>' % _xml_attr(object_name),
           '    <metadata key="extruder" value="%d"/>' % int(parts3mf[0]['extruder'] if parts3mf else 1)]
    for i, p in zip(ids, parts3mf):
        cfg += ['    <part id="%d" subtype="normal_part">' % i,
                '      <metadata key="name" value="%s"/>' % _xml_attr(p['name']),
                '      <metadata key="extruder" value="%d"/>' % int(p['extruder']),
                '    </part>']
    cfg += ['  </object>', '</config>']
    ctypes = ('<?xml version="1.0" encoding="UTF-8"?>\n'
              '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
              '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
              '<Default Extension="model" ContentType="application/vnd.ms-package.3dmanufacturing-3dmodel+xml"/>'
              '<Default Extension="config" ContentType="text/xml"/></Types>')
    rels = ('<?xml version="1.0" encoding="UTF-8"?>\n'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Target="/3D/3dmodel.model" Id="rel0" '
            'Type="http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel"/></Relationships>')
    with zipfile.ZipFile(path, 'w', zipfile.ZIP_DEFLATED) as z:
        z.writestr('[Content_Types].xml', ctypes)
        z.writestr('_rels/.rels', rels)
        z.writestr('3D/3dmodel.model', '\n'.join(lines))
        z.writestr('Metadata/model_settings.config', '\n'.join(cfg))


def _closest_pair(a, b, chunk=400, max_a=6000, max_b=20000):
    rng = np.random.default_rng(1)
    if len(a) > max_a:
        a = a[rng.choice(len(a), max_a, replace=False)]
    if len(b) > max_b:
        b = b[rng.choice(len(b), max_b, replace=False)]
    best = (1e30, None, None)
    for k in range(0, len(a), chunk):
        d2 = ((a[k:k + chunk, None, :] - b[None, :, :]) ** 2).sum(-1)
        i = np.unravel_index(np.argmin(d2), d2.shape)
        if d2[i] < best[0]:
            best = (d2[i], a[k + i[0]], b[i[1]])
    return np.sqrt(best[0]), best[1], best[2]


def attach_floating(parts, max_mm=4.0, rod_d=1.6, touch_mm=0.2, notes=None):
    """Connect loose pieces (e.g. a side chain whose stick ends just off the
    cartoon) to the rest of the model with a short rod. Pieces farther than
    max_mm are left alone (reported as floating)."""
    import manifold3d as mf
    notes = [] if notes is None else notes
    allm = _all_solids(parts)
    if allm is None:
        return 0
    bodies = allm.decompose()
    if len(bodies) < 2:
        return 0
    # cluster touching bodies
    n = len(bodies)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    boxes = [np.asarray(b.bounding_box(), float) for b in bodies]
    for i in range(n):
        for j in range(i + 1, n):
            bi, bj = boxes[i], boxes[j]
            if np.any(bi[:3] - touch_mm > bj[3:]) or np.any(bj[:3] - touch_mm > bi[3:]):
                continue
            if bodies[i].min_gap(bodies[j], touch_mm) < touch_mm:
                parent[find(i)] = find(j)
    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(i)
    cl = sorted(clusters.values(), key=lambda idx: -sum(bodies[i].volume() for i in idx))
    main_v = np.vstack([solid_mesh(bodies[i])[0] for i in cl[0]])
    attached = 0
    r = rod_d / 2.0
    for idx in cl[1:]:
        fv = np.vstack([solid_mesh(bodies[i])[0] for i in idx])
        d, pf, pm = _closest_pair(fv, main_v)
        if d > max_mm:
            continue
        axis = pm - pf
        L = np.linalg.norm(axis)
        axis = axis / L if L > 1e-6 else np.array([0, 0, 1.0])
        rod = _cyl_between(pf - axis * r, pm + axis * r, r, 16)
        # rod takes the loose piece's color
        fb_mesh = solid_mesh(bodies[idx[0]])
        col = main_color(fb_mesh[1], fb_mesh[2])
        target = next((p for p in parts if p.color == col and p.solid is not None),
                      max(parts, key=lambda p: p.volume()))
        t = tint(rod, target.color)
        if t is not None:
            target.solid = target.solid + t
            attached += 1
            main_v = np.vstack([main_v, fv])
    if attached:
        notes.append('attached %d loose piece(s) to the model with %.1f mm rods' % (attached, rod_d))
    return attached


def separate_colors(parts, clearance=0.15, notes=None):
    """For single-color printers: one piece per color that can be printed
    separately and assembled. Same-color parts are merged; bigger pieces get
    sockets (with clearance) where smaller pieces plug in."""
    import manifold3d as mf
    notes = [] if notes is None else notes
    merged = {}
    for p in parts:
        if p.solid is None:
            merged.setdefault(p.color, []).append(p)
            continue
        merged.setdefault(p.color, []).append(p)
    out = []
    for col, plist in merged.items():
        solids = [p.solid for p in plist if p.solid is not None]
        if solids:
            u = mf.Manifold.batch_boolean(solids, mf.OpType.Add) if len(solids) > 1 else solids[0]
            out.append(Part(col, solid=tint(u, col) or u))
        for p in plist:
            if p.solid is None:
                out.append(p)
    out.sort(key=lambda p: -p.volume())
    for i, big in enumerate(out):
        if big.solid is None:
            continue
        for small in out[i + 1:]:
            if small.solid is not None:
                big.solid = big.solid - dilate(small.solid, clearance)
    pieces = []
    for p in out:
        n = len(p.solid.decompose()) if p.solid is not None else 1
        pieces.append(n)
    notes.append('single-color printing: %d colors as separate pieces; smaller pieces '
                 'plug into sockets (%.2f mm clearance)' % (len(out), clearance))
    return out


def part_output(p):
    """Welded mesh for writing: (va, ta, triangle colors uint8). Uses the
    solid's own vertex merge info so color seams don't open the surface."""
    if p.solid is None:
        va, ta, vc = p.mesh
        tc = triangle_colors(ta, vc)
        va2, ta2, _ = weld_vertices(va, ta, vc)
        if len(ta2) == len(ta):
            return va2, ta2, tc
        return va, ta, tc
    me = p.solid.to_mesh()
    vp = np.asarray(me.vert_properties, np.float64)
    ta = np.asarray(me.tri_verts, np.int64)
    if vp.shape[1] >= 6:
        vc = np.clip(np.round(vp[:, 3:6] * 255), 0, 255).astype(np.uint8)
    else:
        vc = np.tile(np.asarray(p.color, np.uint8), (len(vp), 1))
    tc = triangle_colors(ta, vc)
    canon = np.arange(len(vp))
    mfrom, mto = np.asarray(me.merge_from_vert, np.int64), np.asarray(me.merge_to_vert, np.int64)
    if len(mfrom):
        canon[mfrom] = mto
    tcan = canon[ta]
    used, inv = np.unique(tcan.ravel(), return_inverse=True)
    return vp[used, :3], inv.reshape(-1, 3), tc
