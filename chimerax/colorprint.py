"""
colorprint: export displayed ChimeraX models for multi-material 3D printing
with their ChimeraX colors carried over.

Registers one ChimeraX command:

  colorprint <file.obj | file-stem> [models <spec>] [split true|false]
             [colorTolerance <0-255>] [singleColorTriangles true|false]
             [closeEnds true|false] [union true|false]

  * Default (split false): writes an OBJ with per-vertex colors
    ("v x y z r g b"), which OrcaSlicer / Anycubic Slicer Next / Bambu Studio
    read on import and map to filament slots.
  * split true: writes one binary STL per color (<stem>_1_0000ff.stl, ...).
    Load them together in the slicer as ONE object with multiple parts and
    give each part a filament. Gives perfectly sharp color boundaries.

Printability (on by default):
  * closeEnds: ChimeraX draws bond and strut cylinders as open tubes (no end
    caps). These are capped so every piece is a closed solid.
  * union: all single-color pieces of the same color (atom spheres, bond
    cylinders, struts) are boolean-merged into one watertight solid per color.
    Needs the manifold3d package: in ChimeraX run  pip install manifold3d
    (without it, union is skipped and capped pieces are exported).

Only what is displayed is exported (like ChimeraX's own STL export).
Coordinates are in Angstroms, same as ChimeraX's STL export.

Install:  in ChimeraX run   open /path/to/colorprint.py
To load it automatically: Favorites > Settings > Startup >
"Execute these commands at startup":  open /path/to/colorprint.py
"""

import os
import numpy as np


# ---------------------------------------------------------------------------
# Pure geometry / color helpers (no ChimeraX imports; unit tested)
# ---------------------------------------------------------------------------

def triangle_colors(ta, vc):
    """Majority color of each triangle's 3 vertices (ties -> first vertex)."""
    c0, c1, c2 = vc[ta[:, 0]], vc[ta[:, 1]], vc[ta[:, 2]]
    eq12 = (c1 == c2).all(axis=1)
    out = c0.copy()
    out[eq12] = c1[eq12]          # vertices 1 and 2 agree, outvote vertex 0
    return out


def single_color_triangles(va, ta, vc):
    """Duplicate vertices only where a vertex touches triangles of different
    colors, so every triangle has one color but the mesh stays welded inside
    each color region. Returns (va, ta, vc)."""
    tc = triangle_colors(ta, vc)
    uniq, cidx = np.unique(tc.view(np.dtype((np.void, 4))).ravel(),
                           return_inverse=True)
    ncol = len(uniq)
    corner_keys = ta.astype(np.int64) * ncol + np.repeat(cidx.reshape(-1, 1), 3, axis=1)
    keys, inverse = np.unique(corner_keys.ravel(), return_inverse=True)
    new_va = va[keys // ncol]
    new_ta = inverse.reshape(-1, 3).astype(np.int32)
    new_vc = np.empty((len(keys), 4), np.uint8)
    new_vc[new_ta.ravel()] = np.repeat(tc, 3, axis=0)
    return new_va, new_ta, new_vc


def merge_similar_colors(colors, counts, tolerance):
    """Greedy clustering of RGB colors (most used first). Returns an index
    array mapping each input color to its representative color."""
    order = np.argsort(-counts)
    reps = []
    mapping = np.empty(len(colors), np.int64)
    for i in order:
        c = colors[i, :3].astype(int)
        for ri in reps:
            if np.abs(colors[ri, :3].astype(int) - c).max() <= tolerance:
                mapping[i] = ri
                break
        else:
            reps.append(i)
            mapping[i] = i
    return mapping


def apply_color_tolerance(vc, tolerance):
    if tolerance <= 0 or len(vc) == 0:
        return vc
    uniq, inv, counts = np.unique(vc, axis=0, return_inverse=True, return_counts=True)
    mapping = merge_similar_colors(uniq, counts, tolerance)
    return uniq[mapping][inv.ravel()]


def write_obj(path, va, ta, vc, comment=''):
    with open(path, 'w') as f:
        f.write('# %s\n' % comment if comment else '# colorprint export\n')
        f.write('# vertex colors as "v x y z r g b" (0-1)\n')
        rgb = vc[:, :3].astype(np.float64) / 255.0
        np.savetxt(f, np.hstack([va.astype(np.float64), rgb]),
                   fmt='v %.4f %.4f %.4f %.4f %.4f %.4f')
        np.savetxt(f, ta.astype(np.int64) + 1, fmt='f %d %d %d')


def write_binary_stl(path, va, ta, comment=''):
    v = va[ta].astype(np.float32)                         # (nt,3,3)
    n = np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0])
    ln = np.linalg.norm(n, axis=1)
    ln[ln == 0] = 1
    n = (n / ln[:, None]).astype(np.float32)
    rec = np.zeros(len(ta), dtype=[('n', '<f4', (3,)), ('v', '<f4', (3, 3)),
                                   ('a', '<u2')])
    rec['n'] = n
    rec['v'] = v
    header = (comment or 'colorprint export').encode('ascii', 'replace')[:80]
    with open(path, 'wb') as f:
        f.write(header.ljust(80, b' '))
        f.write(np.uint32(len(ta)).tobytes())
        f.write(rec.tobytes())


def split_by_color(va, ta, vc):
    """Yield (rgba, va, ta) per distinct triangle color, largest first."""
    tc = triangle_colors(ta, vc)
    uniq, inv, counts = np.unique(tc, axis=0, return_inverse=True, return_counts=True)
    inv = inv.ravel()
    for ci in np.argsort(-counts):
        sel = ta[inv == ci]
        used, remap = np.unique(sel.ravel(), return_inverse=True)
        yield uniq[ci], va[used], remap.reshape(-1, 3).astype(np.int32)


def weld_vertices(va, ta, vc, tol=1e-4):
    """Merge vertices at the same position (ChimeraX duplicates vertices at
    seams for shading), so open edges can be detected."""
    key = np.round(va / tol).astype(np.int64)
    _, first, inv = np.unique(key, axis=0, return_index=True, return_inverse=True)
    inv = inv.ravel()
    nta = inv[ta]
    # drop triangles that collapsed to a line or point
    ok = ((nta[:, 0] != nta[:, 1]) & (nta[:, 1] != nta[:, 2]) & (nta[:, 0] != nta[:, 2]))
    return va[first], nta[ok].astype(np.int32), vc[first]


def close_holes(va, ta, vc):
    """Cap open boundary loops (e.g. ChimeraX bond cylinders, which are drawn
    without end caps) with a triangle fan to a center point, giving closed
    solids. Returns (va, ta, vc, number_of_loops_capped)."""
    if len(ta) == 0:
        return va, ta, vc, 0
    e = np.concatenate([ta[:, [0, 1]], ta[:, [1, 2]], ta[:, [2, 0]]])
    n = len(va)
    fwd = e[:, 0].astype(np.int64) * n + e[:, 1]
    rev = e[:, 1].astype(np.int64) * n + e[:, 0]
    boundary = e[~np.isin(fwd, rev)]          # directed edges with no partner
    if len(boundary) == 0:
        return va, ta, vc, 0
    nxt = {}
    for a, b in boundary:
        nxt.setdefault(int(a), int(b))
    new_v, new_c, new_t = [], [], []
    visited = set()
    nv = n
    for start in list(nxt):
        if start in visited:
            continue
        loop, cur = [], start
        while cur in nxt and cur not in visited:
            visited.add(cur)
            loop.append(cur)
            cur = nxt[cur]
        if cur != start or len(loop) < 3:
            continue                           # not a clean closed loop
        center = va[loop].mean(axis=0)
        new_v.append(center)
        new_c.append(vc[loop[0]])
        for i, a in enumerate(loop):
            b = loop[(i + 1) % len(loop)]
            new_t.append((b, a, nv))           # reversed to match orientation
        nv += 1
    if not new_v:
        return va, ta, vc, 0
    va = np.concatenate([va, np.array(new_v, va.dtype)])
    vc = np.concatenate([vc, np.array(new_c, vc.dtype)])
    ta = np.concatenate([ta, np.array(new_t, np.int32)])
    return va, ta, vc, len(new_v)


def is_watertight(ta):
    e = np.sort(np.concatenate([ta[:, [0, 1]], ta[:, [1, 2]], ta[:, [2, 0]]]), axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return bool((counts == 2).all())


def union_same_color_pieces(pieces, log_lines):
    """Boolean-union all single-color pieces of the same color into one
    watertight solid using manifold3d. Multi-color pieces (ribbons, surfaces)
    are passed through unchanged. Returns new piece list."""
    try:
        import manifold3d as mf
    except ImportError:
        log_lines.append('Union skipped: manifold3d is not installed. In ChimeraX run '
                         '"pip install manifold3d", then rerun colorprint.')
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
            # Drop zero-volume slivers left where stick and sphere surfaces
            # touch exactly (same radius), which slicers flag as errors.
            bodies = [b for b in u.decompose() if b.volume() > 1e-3]
            if bodies:
                u = mf.Manifold.compose(bodies)
            mesh = u.to_mesh()
            uva = np.asarray(mesh.vert_properties)[:, :3].astype(np.float32)
            uta = np.asarray(mesh.tri_verts).astype(np.int32)
            out.append((uva, uta, np.tile(np.array(color, np.uint8), (len(uva), 1))))
        for va, ta in raw:
            out.append((va, ta, np.tile(np.array(color, np.uint8), (len(va), 1))))
    log_lines.append('Merged sticks/struts into one solid per color'
                     + (' (%d pieces could not be merged and were kept as-is)' % failed
                        if failed else ''))
    return out


def combine_colored_pieces(pieces):
    """pieces: list of (va, ta, vc). Returns combined (va, ta, vc)."""
    if not pieces:
        return (np.zeros((0, 3), np.float32), np.zeros((0, 3), np.int32),
                np.zeros((0, 4), np.uint8))
    vas, tas, vcs, off = [], [], [], 0
    for va, ta, vc in pieces:
        vas.append(va)
        tas.append(ta + off)
        vcs.append(vc)
        off += len(va)
    return np.concatenate(vas), np.concatenate(tas), np.concatenate(vcs)


def hex_color(rgba):
    return '%02x%02x%02x' % tuple(int(x) for x in rgba[:3])


# ---------------------------------------------------------------------------
# ChimeraX scene traversal
# ---------------------------------------------------------------------------

def drawing_colored_geometry(d):
    """Expand one drawing (all displayed instances) into scene-coordinate
    pieces [(va, ta, vc)], carrying vertex, instance or single colors."""
    va, ta = d.vertices, d.masked_triangles
    if va is None or ta is None or len(ta) == 0:
        return []
    positions = d.get_scene_positions(displayed_only=True)
    npos = len(positions)
    if npos == 0:
        return []
    vcol = d.vertex_colors
    icol = d.get_colors(displayed_only=True) if vcol is None else None
    nic = len(icol) if icol is not None else 0
    ta = np.asarray(ta, np.int32)
    pieces = []
    for i, p in enumerate(positions):
        pva = va if p.is_identity() else p * va
        if vcol is not None:
            pvc = np.asarray(vcol, np.uint8)
        else:
            # Parent instances multiply child instances with child index
            # varying fastest, so instance color index is i % nic.
            c = icol[i % nic] if nic else d.color
            pvc = np.tile(np.asarray(c, np.uint8), (len(va), 1))
        pieces.append((np.asarray(pva, np.float32), ta, pvc))
    return pieces


def collect_pieces(models):
    drawings = []
    seen = set()
    for m in models:
        if not m.visible:
            continue
        for d in m.all_drawings(displayed_only=True):
            if d not in seen:
                seen.add(d)
                drawings.append(d)
    pieces, skipped = [], []
    for d in drawings:
        if not (d.display and d.parents_displayed):
            continue
        if d.vertices is None or d.num_masked_triangles == 0:
            continue
        if getattr(d, 'texture', None) is not None or getattr(d, 'multitexture', None):
            skipped.append('%s (texture: labels/images are not exported)' % d.name)
            continue
        if d.display_style != d.Solid:
            skipped.append('%s (mesh/dot style, only solid surfaces are exported)' % d.name)
            continue
        pieces.extend(drawing_colored_geometry(d))
    return pieces, skipped


def prepare_pieces(pieces, close=True, union=True, log_lines=None):
    """Weld, cap open tubes and optionally union same-color pieces."""
    log_lines = [] if log_lines is None else log_lines
    if close:
        fixed, capped = [], 0
        for va, ta, vc in pieces:
            va, ta, vc = weld_vertices(va, ta, vc)
            va, ta, vc, n = close_holes(va, ta, vc)
            capped += n
            fixed.append((va, ta, vc))
        pieces = fixed
        if capped:
            log_lines.append('Capped %d open ends (bond/strut cylinders) to make closed solids'
                             % capped)
    if union:
        pieces = union_same_color_pieces(pieces, log_lines)
    return pieces


def collect_scene(models, close=True, union=False, log_lines=None):
    pieces, skipped = collect_pieces(models)
    pieces = prepare_pieces(pieces, close, union, log_lines)
    return combine_colored_pieces(pieces), skipped


# ---------------------------------------------------------------------------
# Command
# ---------------------------------------------------------------------------

def colorprint(session, path, models=None, split=False, color_tolerance=0,
               single_color_triangles_=True, close_ends=True, union=True):
    from chimerax.core.errors import UserError
    log = session.logger
    if models is None:
        models = session.models.list()
    prep_lines = []
    (va, ta, vc), skipped = collect_scene(models, close_ends, union, prep_lines)
    for s in skipped:
        log.warning('colorprint skipped %s' % s)
    if len(ta) == 0:
        raise UserError('colorprint: nothing displayed to export.')

    # Ignore transparency: e.g. semi-transparent ribbon tethers would otherwise
    # be counted as a separate color from the opaque atoms of the same RGB.
    vc = vc.copy()
    vc[:, 3] = 255
    vc = apply_color_tolerance(vc, color_tolerance)
    path = os.path.expanduser(path)
    stem, ext = os.path.splitext(path)

    colors = list(split_by_color(va, ta, vc))
    summary = ', '.join('#%s (%d tri)' % (hex_color(c), len(t)) for c, _, t in colors)
    lines = ['colorprint: %d triangles, %d colors: %s' % (len(ta), len(colors), summary)]
    lines.extend(prep_lines)

    if split:
        if ext.lower() in ('.stl', '.obj'):
            path = stem
        written = []
        for i, (c, cva, cta) in enumerate(colors, start=1):
            fn = '%s_%d_%s.stl' % (path, i, hex_color(c))
            write_binary_stl(fn, cva, cta, 'colorprint color #%s' % hex_color(c))
            written.append(os.path.basename(fn))
        lines.append('Wrote %d STL files to %s: %s'
                     % (len(written), os.path.dirname(fn) or '.', ', '.join(written)))
        lines.append('Load them together and answer Yes to "load as a single object '
                     'with multiple parts", then assign a filament per part.')
    else:
        if ext.lower() != '.obj':
            path = stem + '.obj'
        if single_color_triangles_:
            va, ta, vc = single_color_triangles(va, ta, vc)
        write_obj(path, va, ta, vc, 'colorprint export from ChimeraX (Angstrom units)')
        lines.append('Wrote %s' % path)
    if len(colors) > 4:
        lines.append('Note: %d colors exceeds a 4-filament setup; recolor or use '
                     'colorTolerance to merge near-identical colors.' % len(colors))
    log.info('\n'.join(lines))


def register_colorprint(logger):
    from chimerax.core.commands import (CmdDesc, register, SaveFileNameArg,
                                        ModelsArg, BoolArg, IntArg)

    def _cmd(session, path, models=None, split=False, color_tolerance=0,
             single_color_triangles=True, close_ends=True, union=True):
        colorprint(session, path, models=models, split=split,
                   color_tolerance=color_tolerance,
                   single_color_triangles_=single_color_triangles,
                   close_ends=close_ends, union=union)

    desc = CmdDesc(required=[('path', SaveFileNameArg)],
                   keyword=[('models', ModelsArg),
                            ('split', BoolArg),
                            ('color_tolerance', IntArg),
                            ('single_color_triangles', BoolArg),
                            ('close_ends', BoolArg),
                            ('union', BoolArg)],
                   synopsis='export displayed models with colors for multi-material 3D printing')
    register('colorprint', desc, _cmd, logger=logger)


# Register when opened in ChimeraX ("open colorprint.py" provides `session`).
if 'session' in globals():
    register_colorprint(session.logger)
    session.logger.info('colorprint command registered. Usage: colorprint ~/Desktop/protein.obj  '
                        '(or add "split true" for one STL per color)')
