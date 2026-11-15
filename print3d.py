"""
print3d.py  -  prepare PyMOL scenes for 3D printing (multi- or single-color).

Load in PyMOL:   run /path/to/print3d.py
(or use the Script Hub panel: Plugin > Script Hub)

Commands
--------
printstyle [helix_width, helix_thickness, strand_width, strand_thickness,
            loop_thickness, stick_radius, sampling, hide_h]
    Thicken cartoon and sticks. Full widths/thicknesses in Angstrom.
    Defaults: helix 3.0 x 1.5, strand 3.0 x 1.0, loops 1.25, sticks 0.7.

split_parts [selection, by, show, color, quality, prefix]
    Separate a structure into objects: by = ligands (polymer + each connected
    ligand, e.g. each glycan), ligands_together, chains or objects; show =
    surface_ligands, surface or keep. Each part is surfaced on its own.

struts [selection, length, loop, radius, color, name, replace]
    Support rods between atoms close in space but far apart along the chain
    (same rules as ChimeraX "struts"). Default: name CA, 8, 60, 0.5.

snapfit host, guest [, kind, diameter, length, clearance, magnet_d, magnet_h, wall]
    Join two separately printed objects at their contact: the guest gets a
    seat shaped to the host, plus a peg (host gets a walled socket) or magnet
    pockets in both. Sizes in mm. Applied by colorprint ..., each=1.
snapfit_clear

colorprint filename [, objects, format, each, size_mm, mm_per_A, orient,
                    filaments, base, single_color, preset, ...]
    Export what is shown, sized, checked and colored:
      format=3mf   parts already assigned to filament slots (and painted where
                   colors mix): opens ready to slice in OrcaSlicer, Anycubic
                   Slicer Next, Bambu Studio. Also obj (vertex colors) or stl
                   (one file per filament).
      size_mm      longest side in mm (or mm_per_A, default 1)
      orient=view  printed the way it looks on screen (screen-down = bed)
      filaments    colors in slots 1..n, e.g. "cyan orange blue green";
                   every PyMOL color snaps to the closest filament
      each=1       one file per object; parts that must fit together (and
                   snapfit connectors) need to be exported in one call
      single_color=1  one file per color for single-filament printers; pieces
                   plug into sockets in the bigger pieces
      base=cut|stand  flat bottom, or a plate with pillars
    It also joins loose pieces (e.g. side chains ending just off the cartoon)
    with short rods and prints a report: size, filament grams, parts thinner
    than min_feature (default 1.5 mm), anything still floating, bed contact.

print3d_preset list | save NAME | load NAME | delete NAME
    Named presets (thickness, struts, split, export, snapfit options) in
    print3d_presets.json. "colorprint file, preset=NAME" uses the export ones.

print3d_install_merge
    Install manifold3d into PyMOL's Python. Needed for merging, connectors,
    bases, single-color sockets and the full report.

Example (protein + glycan, 4-color printer)
    printstyle 3.0, 1.5, 3.0, 1.0, 1.25, 0.7
    color cyan, polymer
    set cartoon_color, cyan
    show sticks, resi 13+14+56+74 and not name N+C+O
    color orange, resi 13+14+56+74 and elem C
    color blue, resi 13+14+56+74 and elem N
    split_parts DN3, ligands_together, surface_ligands, green
    struts DN3_polymer and name CA, 8, 60, 0.5, cyan
    snapfit DN3_polymer, DN3_ligands, peg
    colorprint ~/Desktop/DN3.3mf, each=1, size_mm=120, filaments=cyan orange blue green, base=stand
"""

import os
import sys
import heapq

import numpy as np
from pymol import cmd, cgo

# Find this file's folder. PyMOL's "run" sets __file__ to PyMOL's own package,
# but the code object always knows the real path.
import inspect as _inspect
_HERE = os.path.dirname(os.path.abspath(_inspect.currentframe().f_code.co_filename))
if not os.path.exists(os.path.join(_HERE, 'print3d_core.py')):
    _HERE = os.path.expanduser('~/MOLtiTool')
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import print3d_core as pc  # noqa: E402

# cross-section size per unit of the PyMOL setting (measured from PyMOL's
# exported cartoon geometry)
_HELIX_FACTOR = 2.0      # cartoon_oval_length/width are half sizes
_STRAND_FACTOR = 1.414   # cartoon_rect_length/width
_LOOP_FACTOR = 2.0       # cartoon_loop_radius is a radius

_struts_registry = {}    # name -> dict(struts=[(i1, i2, xyz1, xyz2, obj1, obj2)], radius, rgb)


def _f(x):
    return float(x)


def _b(x):
    if isinstance(x, str):
        return x.strip().lower() not in ('0', 'false', 'no', 'off', '')
    return bool(x)


def _msg(lines):
    for line in lines:
        print(' print3d: ' + line)


# ---------------------------------------------------------------------------
# printstyle
# ---------------------------------------------------------------------------

def printstyle(helix_width=3.0, helix_thickness=1.5, strand_width=3.0,
               strand_thickness=1.0, loop_thickness=1.25, stick_radius=0.7,
               sampling=14, hide_h=1, _self=cmd):
    """
DESCRIPTION

    Thicken cartoon and sticks for 3D printing. Sizes are full widths and
    thicknesses in Angstroms.

USAGE

    printstyle [helix_width [, helix_thickness [, strand_width [,
               strand_thickness [, loop_thickness [, stick_radius [,
               sampling [, hide_h ]]]]]]]]
    """
    s = _self
    s.set('cartoon_oval_length', _f(helix_width) / _HELIX_FACTOR)
    s.set('cartoon_oval_width', _f(helix_thickness) / _HELIX_FACTOR)
    s.set('cartoon_rect_length', _f(strand_width) / _STRAND_FACTOR)
    s.set('cartoon_rect_width', _f(strand_thickness) / _STRAND_FACTOR)
    s.set('cartoon_loop_radius', _f(loop_thickness) / _LOOP_FACTOR)
    s.set('stick_radius', _f(stick_radius))
    s.set('cartoon_sampling', int(float(sampling)))
    s.set('cartoon_fancy_helices', 0)
    if _b(hide_h):
        s.hide('everything', 'hydro')
    s.rebuild()
    _msg(['cartoon helix %.2f x %.2f, strand %.2f x %.2f, loops %.2f, sticks radius %.2f A'
          % (_f(helix_width), _f(helix_thickness), _f(strand_width),
             _f(strand_thickness), _f(loop_thickness), _f(stick_radius))])


# ---------------------------------------------------------------------------
# struts (algorithm follows ChimeraX struts.py)
# ---------------------------------------------------------------------------

def _bond_graph(selection, s):
    """Adjacency over all atoms of the objects containing selection."""
    model = s.get_model('byobject (%s)' % selection)
    xyz = np.array([a.coord for a in model.atom], float)
    keys = [(a.model, a.index) for a in model.atom]
    adj = [[] for _ in model.atom]
    for b in model.bond:
        i, j = b.index
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        adj[i].append((j, d))
        adj[j].append((i, d))
    return keys, xyz, adj


def _connections(sel_idx, adj, dmax):
    """Through-bond distance (<= dmax) between selected atoms."""
    sel_set = set(sel_idx)
    con = []
    for a in sel_idx:
        dist = {a: 0.0}
        heap = [(0.0, a)]
        while heap:
            d, i = heapq.heappop(heap)
            if d > dist.get(i, 1e30):
                continue
            for j, dj in adj[i]:
                nd = d + dj
                if nd <= dmax and nd < dist.get(j, 1e30):
                    dist[j] = nd
                    heapq.heappush(heap, (nd, j))
        for j, d in dist.items():
            if j != a and j in sel_set:
                con.append((a, j, d))
    return con


def _short_connection(i1, i2, dmax, sc):
    adist = {i1: 0.0}
    bndry = {i1}
    while bndry:
        i = bndry.pop()
        d = adist[i]
        for n, dc in sc.get(i, []):
            dn = d + dc
            if dn <= dmax:
                if n == i2:
                    return True
                if n not in adist or dn < adist[n]:
                    adist[n] = dn
                    bndry.add(n)
    return False


def _strut_cgo(struts, radius, rgb):
    obj = []
    for st in struts:
        p1, p2 = st[2], st[3]
        obj += [cgo.SAUSAGE] + list(map(float, p1)) + list(map(float, p2)) + \
               [float(radius)] + list(rgb) + list(rgb)
    return obj


def _draw_struts(name, s):
    entry = _struts_registry.get(name)
    s.delete(name)
    if entry and entry['struts']:
        s.load_cgo(_strut_cgo(entry['struts'], entry['radius'], entry['rgb']), name)


def struts(selection='name CA', length=8.0, loop=60.0, radius=0.5,
           color='grey70', name='struts', replace=1, _self=cmd):
    """
DESCRIPTION

    Add support rods (struts) for 3D printing between atoms in selection that
    are within "length" Angstroms in space but more than "loop" Angstroms
    apart through bonds (or through struts already added).

USAGE

    struts [selection [, length [, loop [, radius [, color [, name [, replace ]]]]]]]
    """
    s = _self
    length, loop, radius = _f(length), _f(loop), _f(radius)
    keys, xyz, adj = _bond_graph(selection, s)
    key_index = {k: i for i, k in enumerate(keys)}
    sel_model = s.get_model(selection)
    sel_idx = [key_index[(a.model, a.index)] for a in sel_model.atom
               if (a.model, a.index) in key_index]
    if len(sel_idx) < 2:
        print(' struts: need at least 2 atoms in selection "%s"' % selection)
        return
    sxyz = xyz[sel_idx]

    # candidate pairs within length, shortest first
    pairs = []
    for a in range(len(sel_idx)):
        d2 = ((sxyz[a + 1:] - sxyz[a]) ** 2).sum(axis=1)
        for b in np.nonzero(d2 <= length * length)[0]:
            pairs.append((d2[b], sel_idx[a], sel_idx[a + 1 + b]))
    pairs.sort()

    sc = {}
    for i1, i2, d in _connections(sel_idx, adj, loop):
        sc.setdefault(i1, []).append((i2, d))
    previous = [] if _b(replace) else list(_struts_registry.get(name, {}).get('struts', []))
    for i1, i2, p1, p2, _m1, _m2 in previous:
        d = float(np.linalg.norm(np.subtract(p1, p2)))
        sc.setdefault(i1, []).append((i2, d))
        sc.setdefault(i2, []).append((i1, d))

    new = []
    for d2, i1, i2 in pairs:
        if not _short_connection(i1, i2, loop, sc):
            d = float(np.sqrt(d2))
            # remember which object each end belongs to (for per-object export)
            new.append((i1, i2, xyz[i1].copy(), xyz[i2].copy(), keys[i1][0], keys[i2][0]))
            sc.setdefault(i1, []).append((i2, d))
            sc.setdefault(i2, []).append((i1, d))

    _last['struts'] = dict(selection=selection, length=length, loop=loop, radius=radius,
                           color=color if isinstance(color, str) else 'grey70')
    all_struts = previous + new
    rgb = s.get_color_tuple(color) if not isinstance(color, (tuple, list)) else tuple(color)
    _struts_registry[name] = dict(struts=all_struts, radius=radius, rgb=rgb)
    _draw_struts(name, s)
    _msg(['created %d struts for "%s" (max length %.3g, loop %.3g), object "%s"'
          % (len(new), selection, length, loop, name)])
    return len(new)


def struts_delete(name='struts', _self=cmd):
    """
DESCRIPTION

    Delete a struts object made by "struts".
    """
    _struts_registry.pop(name, None)
    _self.delete(name)


# ---------------------------------------------------------------------------
# glycan surface
# ---------------------------------------------------------------------------

GLYCAN_DEFAULT = 'organic and not solvent'


def glycan_surface(selection=GLYCAN_DEFAULT, name='glycan_surf', color='green',
                   quality=1, hide_source=1, _self=cmd):
    """
DESCRIPTION

    Copy selection (default: all non-polymer organic molecules, e.g. glycans)
    into its own object and show its surface alone, so it is not cut into by
    the protein surface.

USAGE

    glycan_surface [selection [, name [, color [, quality [, hide_source ]]]]]
    """
    s = _self
    if s.count_atoms(selection) == 0:
        print(' glycan_surface: no atoms in "%s"' % selection)
        return
    s.delete(name)
    s.create(name, selection)
    s.hide('everything', name)
    s.set('surface_quality', int(float(quality)), name)
    s.color(color, name)
    s.show('surface', name)
    if _b(hide_source):
        s.hide('everything', '(%s) and not %s' % (selection, name))
    _msg(['surface for %d atoms in object "%s"' % (s.count_atoms(name), name)])


# ---------------------------------------------------------------------------
# split into parts (ligands / chains / objects)
# ---------------------------------------------------------------------------

SPLIT_MODES = ('ligands', 'ligands_together', 'chains', 'objects')
SHOW_MODES = ('keep', 'surface_ligands', 'surface')


def _copy_part(name, sel, src, s):
    s.delete(name)
    s.create(name, sel)
    # create copies atom colors/representations, not object settings
    try:
        for idx, _typ, val in s.get_object_settings(src):
            s.set(idx, val, name)
    except Exception:
        pass
    return name


def _ligand_groups(sel, s):
    """Connected groups of atoms in sel, using only bonds inside sel (so an
    N-glycan bonded to Asn is still separated from the protein)."""
    m = s.get_model(sel)
    n = len(m.atom)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    for b in m.bond:
        a, c = b.index
        parent[find(a)] = find(c)
    groups = {}
    for i, a in enumerate(m.atom):
        groups.setdefault(find(i), []).append(a)
    out = []
    for atoms in sorted(groups.values(), key=lambda g: (g[0].chain, int(''.join(
            ch for ch in g[0].resi if ch.isdigit()) or 0))):
        out.append(atoms)
    return out


def _reassign_struts(part_names, s):
    """Struts made before splitting point at the original object; move each
    strut end to the part that now holds an atom at that position."""
    if not _struts_registry or not part_names:
        return
    def key(xyz):   # integer key, same for float32 and float64 coordinates
        return tuple(int(round(float(v) * 100)) for v in xyz)
    lookup = {}
    for p in part_names:
        for xyz in s.get_coords(p) if s.count_atoms(p) else []:
            lookup[key(xyz)] = p
    for entry in _struts_registry.values():
        fixed = []
        for i1, i2, p1, p2, m1, m2 in entry['struts']:
            m1 = lookup.get(key(p1), m1)
            m2 = lookup.get(key(p2), m2)
            fixed.append((i1, i2, p1, p2, m1, m2))
        entry['struts'] = fixed


def _safe(name):
    return ''.join(c if (c.isalnum() or c in '_-') else '_' for c in name)


def split_parts(selection='all', by='ligands', show='surface_ligands', color='',
                quality=1, prefix='', _self=cmd):
    """
DESCRIPTION

    Separate a structure into parts, each its own object, so each can be
    surfaced alone (not cut by the others) and exported as its own file.
    The original object is disabled; the parts keep its colors and look.

USAGE

    split_parts [selection [, by [, show [, color [, quality [, prefix ]]]]]]

ARGUMENTS

    by = ligands: polymer + one object per connected ligand (e.g. each glycan)
         ligands_together: polymer + one object with all ligands
         chains: one object per chain
         objects: no copying; use the objects as they are
    show = keep: keep the current look
           surface_ligands: show ligand parts as surfaces, keep the rest
           surface: show every part as a surface
    color = color for the surfaced parts (blank keeps current colors)

EXAMPLES

    split_parts                       # protein + each ligand, ligands surfaced
    split_parts all, chains, keep     # one object per chain
    split_parts DN3, ligands_together, surface_ligands, green
    """
    s = _self
    by, show = str(by).strip().lower(), str(show).strip().lower()
    if by not in SPLIT_MODES:
        print(' split_parts: by must be one of %s' % ', '.join(SPLIT_MODES))
        return
    if show not in SHOW_MODES:
        print(' split_parts: show must be one of %s' % ', '.join(SHOW_MODES))
        return
    _last['split'] = dict(by=by, show=show, color=color, quality=quality)
    objs = [o for o in s.get_object_list('(%s)' % selection)
            if o not in _struts_registry]
    if not objs:
        print(' split_parts: no atoms in "%s"' % selection)
        return
    parts = []          # (name, is_ligand)
    for obj in objs:
        base = _safe(prefix or obj)
        within = '%s and (%s)' % (obj, selection)
        if by == 'objects':
            parts.append((obj, s.count_atoms(within + ' and polymer') == 0))
            continue
        made = []
        if by == 'chains':
            for ch in s.get_chains(within):
                csel = '%s and chain "%s"' % (within, ch)
                made.append((_copy_part('%s_%s' % (base, _safe(ch) or 'nochain'), csel, obj, s),
                             s.count_atoms(csel + ' and polymer') == 0))
        else:
            psel = within + ' and polymer'
            if s.count_atoms(psel):
                made.append((_copy_part(base + '_polymer', psel, obj, s), False))
            lsel = within + ' and organic'
            if s.count_atoms(lsel):
                if by == 'ligands_together':
                    made.append((_copy_part(base + '_ligands', lsel, obj, s), True))
                else:
                    for k, atoms in enumerate(_ligand_groups(lsel, s), start=1):
                        ids = '+'.join(str(a.index) for a in atoms)
                        made.append((_copy_part('%s_lig%d' % (base, k),
                                                '%s and index %s' % (obj, ids), obj, s), True))
        if made:
            s.disable(obj)
            parts += made
    _reassign_struts([n for n, _ in parts], s)
    for name, is_lig in parts:
        if show == 'surface' or (show == 'surface_ligands' and is_lig):
            s.hide('everything', name)
            s.set('surface_quality', int(float(quality)), name)
            if str(color).strip():
                s.color(color, name)
            s.show('surface', name)
    _msg(['parts: ' + ', '.join(n for n, _ in parts),
          'tip: export with "One file per object" (colorprint ..., each=1)'])
    return [n for n, _ in parts]


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

_connectors = []          # list of dict(host, guest, kind, ...)
_last = {}                # last-used options per section (for presets)

EXPORT_DEFAULTS = dict(objects='', format='3mf', each=0, size_mm=0.0, mm_per_A=1.0,
                       orient='view', filaments='', base='none', base_depth=1.0,
                       pillars=3, pillar_d=3.0, min_feature=1.5, attach=1, attach_mm=4.0,
                       single_color=0, socket_clearance=0.15,
                       union=1, close_ends=1,
                       round_ends=1, segments=24, color_tolerance=0)


def _scene_vrml(s):
    """VRML of everything shown, in model coordinates (independent of view)."""
    view = s.get_view()
    try:
        ext = np.array(s.get_extent('enabled'), float)
        far = float(np.abs(ext).max() + np.linalg.norm(ext[1] - ext[0]) + 100.0)
    except Exception:
        far = 1000.0
    try:
        s.set_view([1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, -far, 0, 0, 0, 1.0, 3 * far, 1])
        return s.get_vrml()
    finally:
        s.set_view(view)


def _vrml_of(s, names, temp_cgo=None):
    restore = None
    if temp_cgo:
        s.load_cgo(temp_cgo[1], temp_cgo[0])
        names = list(names) + [temp_cgo[0]]
    if names:
        restore = s.get_names('objects', enabled_only=1)
        s.disable('all')
        for n in names:
            s.enable(n)
    try:
        return _scene_vrml(s)
    finally:
        if restore is not None:
            s.disable('all')
            for n in restore:
                s.enable(n)
        if temp_cgo:
            s.delete(temp_cgo[0])


def _parse_colors(text, s):
    """'cyan orange 0x3366ff #00ff00' -> list of (name, rgb uint8)."""
    out = []
    for tok in str(text).replace(',', ' ').split():
        t = tok.strip()
        if not t:
            continue
        if t.startswith('#'):
            t = '0x' + t[1:]
        try:
            rgb = s.get_color_tuple(t)
            if rgb is None or rgb == -1:
                raise ValueError
        except Exception:
            print(' print3d: unknown color "%s" ignored' % tok)
            continue
        out.append((tok, pc.to_uint8(rgb)))
    return out


_named = {}


def _color_name(rgb, palette_names):
    key = tuple(int(x) for x in rgb)
    for name, c in palette_names:
        if tuple(int(x) for x in c) == key:
            return name
    if not _named:
        try:
            for name, idx in cmd.get_color_indices():
                c = tuple(int(round(v * 255)) for v in cmd.get_color_tuple(idx))
                if c not in _named and not name[:1].isdigit() and '_' not in name:
                    _named[c] = name
        except Exception:
            _named[None] = None
    return _named.get(key) or '#' + pc.hex_color(rgb)


def _near_atom(s, names, placement):
    """centroid in mm -> 'TRP 56 A' of the nearest atom (for the report)."""
    sel = ' or '.join(n for n in names if n not in _struts_registry) or 'none'
    try:
        m = s.get_model(sel)
    except Exception:
        return None
    if not m.atom:
        return None
    xyz = np.array([a.coord for a in m.atom], float)

    def near(c_mm):
        c = placement.inverse(c_mm)
        a = m.atom[int(np.argmin(((xyz - c) ** 2).sum(1)))]
        return '%s %s %s' % (a.resn, a.resi, a.chain)
    return near


def _thin_features(s, names, radii, scale, min_mm):
    """Warn about rods/cartoon thinner than min_mm once printed."""
    warn = []
    obj = next((n for n in names if n not in _struts_registry), None)
    stick_r = float(s.get('stick_radius', obj)) if obj else None
    strut_r = {round(e['radius'], 3) for e in _struts_registry.values()}
    seen = {}
    for typ, r in radii:
        key = round(r, 3)
        if key in seen:
            continue
        if key in strut_r:
            label = 'struts'
        elif stick_r is not None and abs(key - stick_r) < 1e-3:
            label = 'sticks'
        elif typ == 'sphere':
            label = 'spheres'
        else:
            label = 'rods (r=%.2f A)' % r
        seen[key] = label
        d = 2 * r * scale
        if d < min_mm:
            warn.append('%s %.2f mm thick (< %.1f): radius %.2f A or %.2f mm/A would reach it'
                        % (label, d, min_mm, min_mm / (2 * scale), min_mm / (2 * r)))
    if obj and s.count_atoms('(%s) and rep cartoon' % ' or '.join(
            n for n in names if n not in _struts_registry)):
        for label, val in (('helix thickness', float(s.get('cartoon_oval_width', obj)) * _HELIX_FACTOR),
                           ('strand thickness', float(s.get('cartoon_rect_width', obj)) * _STRAND_FACTOR),
                           ('loop thickness', float(s.get('cartoon_loop_radius', obj)) * _LOOP_FACTOR)):
            if val * scale < min_mm:
                warn.append('cartoon %s %.2f mm (< %.1f): set it to %.2f A or scale %.2f mm/A'
                            % (label, val * scale, min_mm, min_mm / scale, min_mm / val))
    return warn


def _parts_fallback(pieces):
    """Without manifold3d: one Part per main color holding the plain mesh."""
    groups = {}
    for va, ta, vc in pieces:
        groups.setdefault(pc.main_color(ta, vc), []).append((va, ta, vc))
    return [pc.Part(col, mesh=pc.combine(pl)) for col, pl in groups.items()]


def _file_for(stem, ext, label, each):
    if not each:
        return stem + ext
    base, folder = os.path.basename(stem), os.path.dirname(stem)
    oname = _safe(label)
    return os.path.join(folder, oname if oname.startswith(base + '_') else
                        '%s_%s' % (base, oname)) + ext


def _write_group(parts, fmt, path_stem, ext, label, slot_of, palette_names, lines):
    """Write one printed object in the chosen format; return written names."""
    written = []
    if fmt == 'obj':
        meshes = [p.get_mesh() for p in parts]
        va, ta, vc = pc.combine([(m[0], m[1], m[2]) for m in meshes])
        fn = path_stem + '.obj'
        pc.write_obj(fn, *pc.single_color_triangles(va, ta, vc),
                     comment='print3d export from PyMOL (mm)')
        written.append(fn)
    elif fmt == 'stl':
        for i, p in enumerate(parts, start=1):
            va, ta, _ = pc.part_output(p)
            slot = slot_of(p.color)
            fn = '%s_%d_%s%s.stl' % (path_stem, slot, _safe(_color_name(p.color, palette_names)),
                                     '_' + p.name if p.name else '')
            pc.write_binary_stl(fn, va, ta, 'print3d filament %d' % slot)
            written.append(fn)
    else:
        p3 = []
        for p in parts:
            va, ta, tc = pc.part_output(p)
            paint = [slot_of(tuple(c)) for c in tc]
            ext_slot = slot_of(p.color)
            p3.append(dict(name='%s%s (filament %d)' % (p.name + ' ' if p.name else '',
                                                          _color_name(p.color, palette_names), ext_slot),
                           va=va, ta=ta, extruder=ext_slot,
                           paint=paint if any(x != ext_slot for x in paint) else None))
        fn = path_stem + '.3mf'
        pc.write_3mf(fn, label, p3)
        written.append(fn)
    return written


def colorprint(filename, objects='', format='', each=-1, size_mm=-1.0, mm_per_A=-1.0,
               orient='', filaments=None, base='', base_depth=-1.0, pillars=-1,
               pillar_d=-1.0, min_feature=-1.0, attach=-1, attach_mm=-1.0,
               single_color=-1, socket_clearance=-1.0, union=-1,
               close_ends=-1, round_ends=-1,
               segments=-1, color_tolerance=-1, split=None, preset='', _self=cmd):
    """
DESCRIPTION

    Export shown geometry with colors, sized and checked for 3D printing.

USAGE

    colorprint filename [, objects [, format [, each [, size_mm [, mm_per_A
               [, orient [, filaments [, base ... ]]]]]]]]

ARGUMENTS

    objects = space-separated objects to export (default: all shown)
    format = 3mf (parts pre-assigned to filament slots, opens ready to slice
             in OrcaSlicer / Anycubic Slicer Next / Bambu Studio), obj (vertex
             colors) or stl (one file per filament). Default: 3mf
    each = 1: one file per object (struts go with the object they brace)
    size_mm = longest side in mm (overrides mm_per_A); mm_per_A default 1
    orient = view: printed the way you see it (screen-down = bed) | model
    filaments = colors loaded in slots 1, 2, 3... e.g. "cyan orange blue green";
                every color is snapped to the closest filament
    base = none | cut (flat bottom, base_depth mm) | stand (plate + pillars)
    min_feature = warn below this thickness in mm (default 1.5)
    single_color = 1: for printers with one filament: every color becomes its
             own file of separate pieces (print each in its color, then glue;
             bigger pieces get sockets for the smaller ones)
    attach = 1: join loose pieces within attach_mm (default 4 mm) to the model
             with a short rod (e.g. side chains ending just off the cartoon)
    preset = name of a saved preset whose export options fill unset arguments
    """
    s = _self
    opts = dict(EXPORT_DEFAULTS)
    if preset:
        opts.update(_load_presets().get(str(preset), {}).get('export', {}))
    given = dict(objects=objects, format=format, each=each, size_mm=size_mm, mm_per_A=mm_per_A,
                 orient=orient, filaments=filaments, base=base, base_depth=base_depth,
                 pillars=pillars, pillar_d=pillar_d, min_feature=min_feature, attach=attach,
                 attach_mm=attach_mm, single_color=single_color,
                 socket_clearance=socket_clearance, union=union,
                 close_ends=close_ends, round_ends=round_ends, segments=segments,
                 color_tolerance=color_tolerance)
    for k, v in given.items():
        if v is None or (isinstance(v, str) and v == '' and k != 'objects') or \
                (not isinstance(v, str) and float(v) == -1):
            continue
        if k == 'objects' and v == '':
            continue
        opts[k] = v
    if split is not None and _b(split):
        opts['format'] = 'stl'
    o = opts
    fmt = str(o['format']).lower().strip('.')
    filename = os.path.abspath(os.path.expanduser(str(filename)))
    stem, ext = os.path.splitext(filename)
    if ext.lower() in ('.obj', '.stl', '.3mf') and not format and not preset:
        fmt = ext.lower().strip('.')
    if fmt not in ('3mf', 'obj', 'stl'):
        print(' colorprint: format must be 3mf, obj or stl')
        return
    each = _b(o['each'])
    size_mm, mm_per_A = float(o['size_mm']), float(o['mm_per_A'])
    segments, tol = int(float(o['segments'])), int(float(o['color_tolerance']))
    names = str(o['objects']).replace(',', ' ').split()
    missing = [n for n in names if n not in s.get_names('all')]
    if missing:
        print(' colorprint: unknown object(s): %s' % ', '.join(missing))
        return
    _last['export'] = {k: o[k] for k in EXPORT_DEFAULTS}
    solid_ok = pc.have_manifold()

    # ---- what to export: one group, or one per object (struts follow) ----
    groups = []                     # (label, names, temp_cgo)
    if each:
        targets = names or s.get_names('objects', enabled_only=1)
        strut_objs = [n for n in targets if n in _struts_registry] if names else \
            [n for n in _struts_registry if n in targets]
        for obj in [n for n in targets if n not in _struts_registry]:
            sub, radius, rgb = [], 0.5, (0.7, 0.7, 0.7)
            for sname in strut_objs:
                e = _struts_registry[sname]
                radius, rgb = e['radius'], e['rgb']
                sub += [st for st in e['struts'] if st[4] == obj and st[5] == obj]
            temp = ('_p3d_struts_tmp', _strut_cgo(sub, radius, rgb)) if sub else None
            groups.append((obj, [obj], temp))
    else:
        label = names[0] if len(names) == 1 else os.path.basename(stem)
        groups.append((label, names, None))

    notes_all = []
    raw = []
    for label, gnames, temp in groups:
        notes, info = [], {}
        pieces = pc.vrml_to_pieces(_vrml_of(s, gnames, temp), segments, _b(o['round_ends']),
                                   _b(o['close_ends']), notes, info)
        if pieces:
            raw.append(dict(label=label, names=gnames or s.get_names('objects', enabled_only=1),
                            pieces=pieces, radii=info.get('radii', []), notes=notes))
    if not raw:
        print(' colorprint: nothing shown to export')
        return

    # ---- shared placement so separately printed parts still match ----
    rot = pc.view_rotation(s.get_view()) if str(o['orient']).lower() == 'view' else None
    allpts = np.vstack([va for g in raw for va, _, _ in g['pieces']])
    pl = pc.make_placement(allpts, rot, mm_per_A, size_mm)

    # ---- filament palette ----
    pal = _parse_colors(o['filaments'], s) if o['filaments'] else []
    pal_rgb = np.array([c for _, c in pal], np.uint8) if pal else None
    orig_map = []
    if pal_rgb is not None:
        allc = np.vstack([vc[:, :3] for g in raw for _, _, vc in g['pieces']])
        orig_map = pc.palette_mapping(allc, pal_rgb)
    for g in raw:
        pieces = pc.transform_pieces(g['pieces'], pl)
        if tol and pal_rgb is None:
            pieces = [(va, ta, pc.apply_color_tolerance(vc, tol)) for va, ta, vc in pieces]
        if pal_rgb is not None:
            pieces = [(va, ta, pc.snap_to_palette(vc, pal_rgb)[0]) for va, ta, vc in pieces]
        g['pieces'] = pieces
        g['parts'] = pc.build_parts(pieces, g['notes']) if (solid_ok and _b(o['union'])) \
            else _parts_fallback(pieces)

    # slot numbers: palette order, else by total volume/triangles
    if pal_rgb is not None:
        def slot_of(c):
            return int(pc.snap_to_palette(np.array([c], np.uint8), pal_rgb)[1][0]) + 1
        pal_names = [(n, c) for n, c in pal]
    else:
        order = []
        for g in raw:
            for p in g['parts']:
                if p.color not in order:
                    order.append(p.color)
        for g in raw:
            for p in g['parts']:
                va, ta, vc = p.get_mesh()
                for c in np.unique(pc.triangle_colors(ta, vc), axis=0):
                    c = tuple(int(x) for x in c)
                    if c not in order:
                        order.append(c)

        def slot_of(c, _order=order):
            c = tuple(int(x) for x in c)
            return _order.index(c) + 1 if c in _order else 1
        pal_names = []

    # ---- connectors (between separately printed objects) ----
    by_label = {g['label']: g for g in raw}
    for con in _connectors:
        gh, gg = by_label.get(con['host']), by_label.get(con['guest'])
        if gh is None or gg is None:
            if each:
                notes_all.append('connector %s/%s skipped (both objects must be exported)'
                                 % (con['host'], con['guest']))
            else:
                notes_all.append('connectors need "One file per object" (each=1)')
            continue
        if not solid_ok:
            notes_all.append('connectors need manifold3d')
            continue
        pc.snapfit(gh['parts'], gg['parts'], None, None, con['kind'], float(con['diameter']),
                   float(con['length']), float(con['clearance']), float(con['wall']),
                   float(con['magnet_d']), float(con['magnet_h']), notes=gh['notes'])

    # ---- join loose pieces (side chains just off the cartoon, etc.) ----
    if solid_ok and _b(o['attach']):
        for g in raw:
            pc.attach_floating(g['parts'], float(o['attach_mm']),
                               max(float(o['min_feature']), 1.2), notes=g['notes'])

    # ---- base on the biggest object ----
    main = max(raw, key=lambda g: sum(len(p.get_mesh()[1]) for p in g['parts']))
    bmode = str(o['base']).lower()
    if bmode in ('cut', 'stand') and not solid_ok:
        notes_all.append('base needs manifold3d')
    elif bmode == 'cut':
        main['parts'] = pc.flat_cut(main['parts'], float(o['base_depth']), main['notes'])
    elif bmode == 'stand':
        stand = pc.pillar_stand(main['parts'], int(float(o['pillars'])), float(o['pillar_d']),
                                notes=main['notes'])
        if stand is not None:
            c = max(main['parts'], key=lambda p: p.volume()).color
            main['parts'].append(pc.Part(c, solid=pc.tint(stand, c), name='stand'))

    single = _b(o['single_color'])
    if single:
        if solid_ok:
            for g in raw:
                g['parts'] = pc.separate_colors(g['parts'], float(o['socket_clearance']),
                                                g['notes'])
        else:
            notes_all.append('single-color pieces without sockets (manifold3d not installed)')

    # ---- report + write ----
    out_files, report = [], []
    report.append('scale %.2f mm per A, %s' % (pl.scale, 'oriented as on screen'
                                                if rot is not None else 'molecule frame'))
    if pal:
        report.append('filaments: ' + ', '.join('%d=%s' % (i + 1, n) for i, (n, _) in enumerate(pal)))
        report.append('   colors mapped: ' + ', '.join(
            '#%s->%d' % (pc.hex_color(c), sl + 1) for c, sl in orig_map))
    for g in raw:
        parts = g['parts']
        near = _near_atom(s, g['names'], pl)
        a = pc.analyze(parts, near) if solid_ok else dict(size=None, bodies=[], grams=None)
        head = '[%s]' % g['label']
        if a.get('size') is not None:
            sz = a['size']
            head += ' %.0f x %.0f x %.0f mm' % (sz[0], sz[1], sz[2])
        if a.get('grams'):
            head += ', about %.0f g if printed solid' % a['grams']
        report.append(head)
        per = []
        for p in parts:
            if p.solid is not None:
                per.append('%s%s on filament %d: %.1f g' % (
                    _color_name(p.color, pal_names), ' (%s)' % p.name if p.name else '',
                    slot_of(p.color), p.volume() / 1000 * pc.PLA_DENSITY))
        if per:
            report.append('   ' + '; '.join(per))
        for w in _thin_features(s, g['names'], g['radii'], pl.scale, float(o['min_feature'])):
            report.append('   thin: ' + w)
        bodies = a.get('bodies') or []
        if len(bodies) > 1:
            fl = bodies[1:]
            report.append('   floating: %d piece(s) not attached to the main body: %s'
                          % (len(fl), '; '.join('%.0f mm3 near %s' % (v, n) for v, _, n in fl[:6])
                                               + (' ...' if len(fl) > 6 else '')))
        if a.get('bed_contact') is not None and bmode != 'stand':
            report.append('   bed contact %.0f mm2%s' % (a['bed_contact'], ' (small: supports '
                          'or base=cut/stand recommended)' if a['bed_contact'] < 30 else ''))
        report += ['   ' + n for n in g['notes'] if not n.startswith('Read ')]
        path_stem = _file_for(stem, '', g['label'], each)
        if single:
            for prt in parts:
                cname = _safe(_color_name(prt.color, pal_names).lstrip('#'))
                nb = len(prt.solid.decompose()) if prt.solid is not None else 1
                if fmt == 'stl':
                    va_, ta_, _ = pc.part_output(prt)
                    fn_ = '%s_%s.stl' % (path_stem, cname)
                    pc.write_binary_stl(fn_, va_, ta_, 'print3d %s' % cname)
                    files = [fn_]
                else:
                    files = _write_group([prt], fmt, '%s_%s' % (path_stem, cname), ext,
                                         '%s %s' % (g['label'], cname), lambda c: 1,
                                         pal_names, report)
                report.append('   print in %s: %s (%d piece%s)' % (
                    _color_name(prt.color, pal_names), os.path.basename(files[0]),
                    nb, '' if nb == 1 else 's'))
                out_files += files
        else:
            out_files += _write_group(parts, fmt, path_stem, ext, g['label'], slot_of,
                                      pal_names, report)
    report += notes_all
    report.append('wrote: ' + ', '.join(os.path.basename(f) for f in out_files))
    if fmt == 'stl' and not single:
        report.append('load the STLs together as ONE object with multiple parts')
    _msg(report)
    _last['report'] = report
    return ', '.join(os.path.basename(f) for f in out_files)


# ---------------------------------------------------------------------------
# snap-fit connectors
# ---------------------------------------------------------------------------

def snapfit(host, guest, kind='peg', diameter=4.0, length=6.0, clearance=0.2,
            magnet_d=6.0, magnet_h=3.0, wall=1.6, _self=cmd):
    """
DESCRIPTION

    Join two separately printed objects (e.g. protein and glycan) at their
    contact: the guest gets a seat cut to the host's shape, plus either a peg
    (host gets a socket) or magnet pockets in both. Applied by
    "colorprint ..., each=1". Sizes in mm (after scaling).

USAGE

    snapfit host, guest [, kind [, diameter [, length [, clearance
            [, magnet_d [, magnet_h [, wall ]]]]]]]

    kind = peg | magnet
    """
    s = _self
    for n in (host, guest):
        if n not in s.get_names('objects'):
            print(' snapfit: unknown object %s' % n)
            return
    global _connectors
    _connectors = [c for c in _connectors if not (c['host'] == host and c['guest'] == guest)]
    con = dict(host=host, guest=guest, kind=str(kind).lower(), diameter=_f(diameter),
               length=_f(length), clearance=_f(clearance), magnet_d=_f(magnet_d),
               magnet_h=_f(magnet_h), wall=_f(wall))
    _connectors.append(con)
    _last['snapfit'] = {k: con[k] for k in con if k not in ('host', 'guest')}
    _msg(['%s between %s (host) and %s: applied when exporting with each=1'
          % (con['kind'], host, guest)])


def snapfit_clear(_self=cmd):
    """Remove all snap-fit connectors."""
    del _connectors[:]
    _msg(['connectors cleared'])


def snapfit_list():
    return ['%s: %s -> %s' % (c['kind'], c['host'], c['guest']) for c in _connectors]


# ---------------------------------------------------------------------------
# presets
# ---------------------------------------------------------------------------

PRESET_FILE = os.path.join(_HERE, 'print3d_presets.json')


def _load_presets():
    import json
    try:
        with open(PRESET_FILE) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _save_presets(data):
    import json
    os.makedirs(os.path.dirname(PRESET_FILE), exist_ok=True)
    with open(PRESET_FILE, 'w') as f:
        json.dump(data, f, indent=2, sort_keys=True)


def current_settings(s=cmd):
    """Snapshot of style (from PyMOL settings) and last-used options."""
    style = dict(helix_width=float(s.get('cartoon_oval_length')) * _HELIX_FACTOR,
                 helix_thickness=float(s.get('cartoon_oval_width')) * _HELIX_FACTOR,
                 strand_width=float(s.get('cartoon_rect_length')) * _STRAND_FACTOR,
                 strand_thickness=float(s.get('cartoon_rect_width')) * _STRAND_FACTOR,
                 loop_thickness=float(s.get('cartoon_loop_radius')) * _LOOP_FACTOR,
                 stick_radius=float(s.get('stick_radius')),
                 sampling=int(float(s.get('cartoon_sampling'))))
    data = dict(style={k: round(v, 3) for k, v in style.items()})
    for k in ('struts', 'split', 'export', 'snapfit'):
        if k in _last:
            data[k] = dict(_last[k])
    return data


def print3d_preset(action='list', name='', _self=cmd):
    """
DESCRIPTION

    Save, load, list or delete named print presets (thickness, struts,
    split, export and snap-fit options), stored in print3d_presets.json.

USAGE

    print3d_preset list
    print3d_preset save, DN3      # current thickness + last-used options
    print3d_preset load, DN3      # applies thickness; export options are
                                  # used by: colorprint file, preset=DN3
    print3d_preset delete, DN3
    """
    s = _self
    action, name = str(action).lower(), str(name)
    data = _load_presets()
    if action == 'list':
        _msg(['presets: ' + (', '.join(sorted(data)) or '(none)')])
        return sorted(data)
    if not name:
        print(' print3d_preset: give a name')
        return
    if action == 'save':
        data[name] = current_settings(s)
        _save_presets(data)
        _msg(['saved preset "%s" to %s' % (name, PRESET_FILE)])
    elif action == 'load':
        p = data.get(name)
        if p is None:
            print(' print3d_preset: no preset "%s"' % name)
            return
        if 'style' in p:
            printstyle(_self=s, **p['style'])
        _msg(['loaded preset "%s"' % name])
        return p
    elif action == 'delete':
        if data.pop(name, None) is not None:
            _save_presets(data)
            _msg(['deleted preset "%s"' % name])


def print3d_install_merge(_self=cmd):
    """
DESCRIPTION

    Install manifold3d (used by colorprint union=1) into PyMOL's Python.
    """
    import subprocess
    candidates = [sys.executable,
                  os.path.join(sys.prefix, 'bin', 'python'),
                  os.path.join(sys.prefix, 'bin', 'python3')]
    for py in candidates:
        if py and os.path.exists(py) and 'python' in os.path.basename(py).lower():
            print(' print3d: installing manifold3d with %s ...' % py)
            r = subprocess.run([py, '-m', 'pip', 'install', 'manifold3d'],
                               capture_output=True, text=True)
            print(r.stdout[-1500:] + r.stderr[-1500:])
            if r.returncode == 0:
                print(' print3d: manifold3d installed. Restart PyMOL if import still fails.')
                return True
    print(' print3d: could not find PyMOL\'s Python. In a terminal, run the pip of '
          'the Python PyMOL uses:  <that python> -m pip install manifold3d')
    return False


cmd.extend('printstyle', printstyle)
cmd.extend('struts', struts)
cmd.extend('struts_delete', struts_delete)
cmd.extend('glycan_surface', glycan_surface)
cmd.extend('split_parts', split_parts)
cmd.extend('colorprint', colorprint)
cmd.extend('snapfit', snapfit)
cmd.extend('snapfit_clear', snapfit_clear)
cmd.extend('print3d_preset', print3d_preset)
cmd.extend('print3d_install_merge', print3d_install_merge)
cmd.auto_arg[0]['struts'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['glycan_surface'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['split_parts'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['snapfit'] = cmd.auto_arg[0]['disable']
cmd.auto_arg[1]['snapfit'] = cmd.auto_arg[0]['disable']
