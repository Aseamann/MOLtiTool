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

surfacestyle [selection, quality, probe, solvent, smooth_edges, color, show]
    Surface resolution (quality -1..4), probe radius, molecular or solvent-
    accessible surface, smooth edges; shows the surface.

print3d_printers
    Printer profiles (build volume, filament slots) for colorprint printer=...

auto_filaments [objects, slots]
    Common PLA colors closest to the scene's colors (colorprint filaments=auto).

ligand_surface [selection, name, color, quality, hide_source]
    Shortcut: copy the ligand(s) into one object and surface it on its own.

split_parts [selection, by, show, color, quality, prefix]
    Separate a structure into objects: by = ligands (polymer + each connected
    ligand, e.g. each substrate or cofactor), ligands_together, chains or
    objects; show =
    surface_ligands, surface or keep. Each part is surfaced on its own.

struts [selection, length, loop, radius, color, name, replace]
    Support rods between atoms close in space but far apart along the chain
    (same rules as ChimeraX "struts"). Default: name CA, 8, 60, 0.5. Color
    blank = each strut takes its chain's color (half and half between chains).
strut_add atom1, atom2 [, name, radius, color]
strut_remove atom1 [, atom2, name]
    Add or remove one strut (an atom, or a residue = its CA). With only
    atom1, strut_remove removes every strut at that atom.
strut_remove_selected [selection]
    Click a strut in the viewer (it is selected), then remove it with this.
struts_label [on]
    Show an ID label on every strut, e.g. "S3 A:LYS45-A:GLU112".

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
      nozzle, layer_height  what to tune for (default 0.4 / 0.2 mm): sets the
                   thinnest-feature warning, fit and socket clearances, pillar
                   size, base sizes in whole layers and rod smoothness
      preview=1    build it all but show it in PyMOL (preview_<name> objects,
                   in place over the molecule) instead of writing files
    It also joins loose pieces (e.g. side chains ending just off the cartoon)
    with short rods and prints a report: size, filament grams, parts thinner
    than min_feature (default 1.5 mm), anything still floating, bed contact.

print3d_preset list | save NAME | load NAME | delete NAME
    Named presets (thickness, struts, split, export, snapfit options) in
    print3d_presets.json. "colorprint file, preset=NAME" uses the export ones.

print3d_install_merge
    Install manifold3d into PyMOL's Python. Needed for merging, connectors,
    bases, single-color sockets and the full report.

Example (protein + bound ligand, 4-color printer)
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
if _HERE not in sys.path:
    sys.path.insert(0, _HERE)

import print3d_core as pc  # noqa: E402

# cross-section size per unit of the PyMOL setting (measured from PyMOL's
# exported cartoon geometry)
_HELIX_FACTOR = 2.0      # cartoon_oval_length/width are half sizes
_STRAND_FACTOR = 1.414   # cartoon_rect_length/width
_LOOP_FACTOR = 2.0       # cartoon_loop_radius is a radius

_struts_registry = {}    # name -> dict(struts=[dict(p1, p2, m1, m2)], radius, color=None|rgb)


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
               sampling=14, hide_h=1, cylinder_helices=0, helix_radius=2.25, show_joins=1,
               _self=cmd):
    """
DESCRIPTION

    Thicken cartoon and sticks for 3D printing. Sizes are full widths and
    thicknesses in Angstroms. Double and aromatic bonds are drawn as one
    full-thickness stick (PyMOL's thin double-bond lines don't print).

USAGE

    printstyle [helix_width [, helix_thickness [, strand_width [,
               strand_thickness [, loop_thickness [, stick_radius [,
               sampling [, hide_h [, cylinder_helices [, helix_radius ]]]]]]]]]]

    cylinder_helices = 1: draw helices as solid cylinders (helix_radius A)
             so no support material can get trapped inside them
    show_joins = 1 (default): show where export joins loops into beta
             strands (object "strand_joins", display only; 0 removes it)
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
    s.set('cartoon_cylindrical_helices', 1 if _b(cylinder_helices) else 0)
    s.set('cartoon_helix_radius', _f(helix_radius))
    _single_bonds(s)
    if _b(hide_h):
        s.hide('everything', 'hydro')
    s.rebuild()
    nj = show_strand_joins(on=show_joins, _self=s)
    _msg(['cartoon helix %.2f x %.2f, strand %.2f x %.2f, loops %.2f, sticks radius %.2f A'
          % (_f(helix_width), _f(helix_thickness), _f(strand_width),
             _f(strand_thickness), _f(loop_thickness), _f(stick_radius))] +
         (['helices as solid cylinders, radius %.2f A' % _f(helix_radius)]
          if _b(cylinder_helices) else []) +
         (['showing %d strand-loop joins (object %s)' % (nj, JOINS_OBJECT)] if nj else []))


# ---------------------------------------------------------------------------
# struts (algorithm follows ChimeraX struts.py)
# ---------------------------------------------------------------------------

def _bond_graph(selection, s):
    """Adjacency over all atoms of the objects containing selection, plus a
    short residue label per atom (e.g. 'A:LYS45')."""
    model = s.get_model('byobject (%s)' % selection)
    xyz = np.array([a.coord for a in model.atom], float)
    keys = [(a.model, a.index) for a in model.atom]
    labels = ['%s:%s%s' % (a.chain, a.resn, a.resi) for a in model.atom]
    adj = [[] for _ in model.atom]
    for b in model.bond:
        i, j = b.index
        d = float(np.linalg.norm(xyz[i] - xyz[j]))
        adj[i].append((j, d))
        adj[j].append((i, d))
    return keys, xyz, adj, labels


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


def _chain_color_lookup(s):
    """Returns color_of(model, xyz) -> rgb: the color of the chain the atom at
    xyz belongs to (the chain's cartoon_color if set, else its most common CA
    color), so struts match their chain."""
    cache = {}

    def model_info(model):
        if model in cache:
            return cache[model]
        from pymol import stored
        stored._p3d = []
        s.iterate(model, 'stored._p3d.append((chain, name, color))')
        rows = stored._p3d
        xyz = s.get_coords(model)
        obj_cc = str(s.get('cartoon_color', model))
        chain_rgb = {}
        by_chain = {}
        for ch, nm, col in rows:
            if nm == 'CA':
                by_chain.setdefault(ch, []).append(col)
        for ch, cols in by_chain.items():
            if obj_cc not in ('default', '-1'):
                chain_rgb[ch] = tuple(s.get_color_tuple(obj_cc))
            else:
                best = max(set(cols), key=cols.count)
                chain_rgb[ch] = tuple(s.get_color_tuple(best))
        info = (np.asarray(xyz, float) if xyz is not None else np.zeros((0, 3)),
                [r[0] for r in rows], [r[2] for r in rows], chain_rgb)
        cache[model] = info
        return info

    def color_of(model, p):
        try:
            xyz, chains, cols, chain_rgb = model_info(model)
        except Exception:
            return (0.7, 0.7, 0.7)
        if not len(xyz):
            return (0.7, 0.7, 0.7)
        i = int(np.argmin(((xyz - np.asarray(p, float)) ** 2).sum(1)))
        ch = chains[i]
        if ch in chain_rgb:
            return chain_rgb[ch]
        return tuple(s.get_color_tuple(cols[i]))
    return color_of


def _ends(st, radius=0.0):
    """Where a strut is drawn. Ends anchored on the cartoon (d1/d2, points on
    its center line) are pulled back by the strut radius, so the round tip
    stops at the middle of the ribbon instead of poking out the far side.
    Other ends sit at the atoms (p1/p2). p1/p2 always stay the atom positions,
    so struts keep matching their atoms."""
    a = np.asarray(st.get('d1', st['p1']), float)
    b = np.asarray(st.get('d2', st['p2']), float)
    d = b - a
    n = float(np.linalg.norm(d))
    r = float(radius)
    if r > 0 and n > 1e-6:
        u = d / n
        pull = min(r, 0.25 * n)          # never shrink a short strut away
        if 'd1' in st:
            a = a + u * pull
        if 'd2' in st:
            b = b - u * pull
    return a, b


def _pkey(xyz):
    return tuple(int(round(float(v) * 100)) for v in xyz)


def _cartoon_anchor_map(models, s, meshes=None):
    """{(model, atom-position key): point on the cartoon's center line} for
    every CA shown as cartoon. PyMOL smooths beta strands flat, so a strand's
    CA atoms zig-zag up to ~1.5 A off the ribbon; a strut ending on the
    center line meets the strand squarely instead of poking out of one face."""
    out = {}
    for model in models:
        if model not in s.get_names('objects') or not s.count_atoms(
                '(%s) and name CA and rep cartoon' % model):
            continue
        reps = {}
        s.iterate(model, 'R[index] = reps', space={'R': reps})
        idx = [k for k, v in reps.items() if v & 32]          # 32 = cartoon bit
        try:
            s.hide('everything', model)          # VRML of the cartoon alone
            for k in range(0, len(idx), 400):
                s.show('cartoon', '%s and index %s' % (model, '+'.join(map(str, idx[k:k + 400]))))
            vr = _vrml_of(s, [model])
        finally:
            s.alter(model, 'reps = R[index]', space={'R': reps})
            s.rebuild(model)
        pieces = pc.vrml_to_pieces(vr, 6, True, False, [], {})
        if not pieces:
            continue
        V = np.vstack([va for va, _, _ in pieces])
        if meshes is not None:
            meshes[model] = V
        rows = []
        s.iterate_state(-1, '(%s) and name CA and rep cartoon' % model,
                        'rows.append((segi, chain, int(resv), x, y, z))',
                        space={'rows': rows, 'int': int})
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        xyz = np.array([r[3:] for r in rows], float)
        for i, r in enumerate(rows):
            prev = i - 1 if i > 0 and rows[i - 1][:2] == r[:2] and rows[i - 1][2] == r[2] - 1 else None
            nxt = i + 1 if i + 1 < len(rows) and rows[i + 1][:2] == r[:2] and \
                rows[i + 1][2] == r[2] + 1 else None
            if prev is None and nxt is None:
                continue
            t = xyz[nxt if nxt is not None else i] - xyz[prev if prev is not None else i]
            nt = np.linalg.norm(t)
            if nt < 1e-6:
                continue
            t /= nt
            d = V - xyz[i]
            along = d @ t
            lat = np.linalg.norm(d - np.outer(along, t), axis=1)
            sel = (np.abs(along) < 0.6) & (lat < 3.0)
            if sel.sum() < 4:
                continue
            c = V[sel].mean(0)
            if np.linalg.norm(c - xyz[i]) < 2.5:
                out[(model, _pkey(xyz[i]))] = c
    return out


def _anchor_struts(struts, s):
    amap = _cartoon_anchor_map(sorted({st['m1'] for st in struts} | {st['m2'] for st in struts}), s)
    moved = 0
    for st in struts:
        for pk, dk, mk in (('p1', 'd1', 'm1'), ('p2', 'd2', 'm2')):
            a = amap.get((st[mk], _pkey(st[pk])))
            if a is not None:
                st[dk] = a
                moved += 1
    return moved


def _strut_cgo(struts, radius, color, s=cmd):
    """CGO rods. color None = each half takes the color of its end's chain."""
    color_of = _chain_color_lookup(s) if color is None else None
    obj = []
    for st in struts:
        p1, p2 = _ends(st, radius)
        if color is not None:
            c1 = c2 = tuple(color)
        else:
            c1, c2 = color_of(st['m1'], st['p1']), color_of(st['m2'], st['p2'])
        if c1 == c2:
            segs = [(p1, p2, c1)]
        else:                          # two halves, one per chain color
            mid = (p1 + p2) / 2
            segs = [(p1, mid, c1), (mid, p2, c2)]
        for a, b, c in segs:
            obj += [cgo.SAUSAGE] + list(map(float, a)) + list(map(float, b)) + \
                   [float(radius)] + list(c) + list(c)
    return obj


def _strut_label(k, st):
    return 'S%d %s-%s' % (k, st.get('l1', '?'), st.get('l2', '?'))


def _draw_struts(name, s):
    """Draw struts as a small molecule object: strut k is residue STR k of
    chain S, two pseudo-atoms (E1, E2) joined by a stick. So clicking a strut
    selects it (and PyMOL prints which one), and it can be labeled."""
    from chempy import models, Atom, Bond
    entry = _struts_registry.get(name)
    was_enabled = name not in s.get_names('objects') or \
        name in s.get_names('objects', enabled_only=1)
    s.delete(name)
    if not entry or not entry['struts']:
        return
    color_of = _chain_color_lookup(s) if entry.get('color') is None else None
    m = models.Indexed()
    colors = []
    for k, st in enumerate(entry['struts'], start=1):
        for j, key in enumerate(('p1', 'p2')):
            a = Atom()
            a.name = 'E%d' % (j + 1)
            a.resn = 'STR'
            a.resi = str(k)
            a.resi_number = k
            a.chain = 'S'
            a.symbol = 'PS'          # pseudo-atom: never counted as ligand/organic
            a.hetatm = 1
            a.coord = [float(x) for x in _ends(st, entry['radius'])[j]]
            m.add_atom(a)
            mdl = st['m1'] if j == 0 else st['m2']
            colors.append(tuple(entry['color']) if color_of is None else color_of(mdl, st[key]))
        b = Bond()
        b.index = [2 * (k - 1), 2 * (k - 1) + 1]
        b.order = 1
        m.add_bond(b)
    s.load_model(m, name)
    if name not in s.get_names('objects'):
        s.load_model(m, name)                # once more (seen after undo in PyMOL 3)
    if name not in s.get_names('objects'):
        print(' struts: could not draw "%s"; run struts again' % name)
        return
    s.hide('everything', name)
    s.show('sticks', name)
    s.set('stick_radius', float(entry['radius']), name)
    s.flag('ignore', name, 'set')            # keep them out of surfaces
    cidx = {}
    for i, c in enumerate(colors):
        hx = '0x%02x%02x%02x' % tuple(int(round(v * 255)) for v in c)
        cidx[(i // 2 + 1, 'E%d' % (i % 2 + 1))] = s.get_color_index(hx)
    s.alter(name, 'color = C.get((int(resv), name), color)', space={'C': cidx, 'int': int})
    s.rebuild(name)
    if entry.get('labels'):
        s.set('label_size', 14, name)
        for k, st in enumerate(entry['struts'], start=1):
            s.label('%s and resi %d and name E1' % (name, k), repr(_strut_label(k, st)))
    if not was_enabled:
        s.disable(name)


def _sync_struts(s, names=None):
    """Struts whose atoms were removed in the viewer (e.g. clicked and the
    Delete key pressed) are dropped from the registry, so they stay gone."""
    dropped = 0
    for name in (list(_struts_registry) if names is None else names):
        entry = _struts_registry.get(name)
        if not entry or name not in s.get_names('objects'):
            continue
        have = set()
        s.iterate(name, 'H.add((int(resv), name))', space={'H': have, 'int': int})
        n = len(entry['struts'])
        if len(have) == 2 * n:
            continue
        keep = [st for k, st in enumerate(entry['struts'], start=1)
                if (k, 'E1') in have and (k, 'E2') in have]
        dropped += n - len(keep)
        entry['struts'] = keep
        _draw_struts(name, s)
    if dropped:
        _msg(['%d strut(s) removed in the viewer were dropped' % dropped])
    return dropped


def _rgb_or_none(color, s):
    """'' / 'auto' / 'chain' -> None (chain colors), else an rgb tuple."""
    if color is None or (isinstance(color, str) and color.strip().lower() in ('', 'auto', 'chain')):
        return None
    if isinstance(color, (tuple, list)):
        return tuple(float(x) for x in color)
    c = str(color).strip()
    if c.startswith('#'):
        c = '0x' + c[1:]
    return tuple(s.get_color_tuple(c))


def struts(selection='name CA', length=8.0, loop=60.0, radius=0.5,
           color='', name='struts', replace=1, anchor='cartoon', _self=cmd):
    """
DESCRIPTION

    Add support rods (struts) for 3D printing between atoms in selection that
    are within "length" Angstroms in space but more than "loop" Angstroms
    apart through bonds (or through struts already added).

    color: blank (default) = each strut takes the color of its chain; give a
    color (e.g. grey70) to use one color for all.
    anchor: cartoon (default) = ends of struts on CA atoms shown as cartoon
    sit on the cartoon's center line (beta strands are smoothed, so their CA
    atoms sit off the ribbon); atom = exactly at the atoms.

USAGE

    struts [selection [, length [, loop [, radius [, color [, name [, replace
           [, anchor ]]]]]]]]

SEE ALSO

    strut_add, strut_remove, struts_delete
    """
    s = _self
    _sync_struts(s)
    length, loop, radius = _f(length), _f(loop), _f(radius)
    keys, xyz, adj, labels = _bond_graph(selection, s)
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

    def key(p):
        return tuple(int(round(float(v) * 100)) for v in p)
    pos_index = {key(xyz[i]): i for i in sel_idx}
    for st in previous:     # existing struts count as connections
        i1, i2 = pos_index.get(key(st['p1'])), pos_index.get(key(st['p2']))
        if i1 is not None and i2 is not None:
            d = float(np.linalg.norm(np.subtract(st['p1'], st['p2'])))
            sc.setdefault(i1, []).append((i2, d))
            sc.setdefault(i2, []).append((i1, d))

    new = []
    for d2, i1, i2 in pairs:
        if not _short_connection(i1, i2, loop, sc):
            d = float(np.sqrt(d2))
            # remember which object each end belongs to (for per-object export)
            new.append(dict(p1=xyz[i1].copy(), p2=xyz[i2].copy(),
                            m1=keys[i1][0], m2=keys[i2][0], l1=labels[i1], l2=labels[i2]))
            sc.setdefault(i1, []).append((i2, d))
            sc.setdefault(i2, []).append((i1, d))

    if str(anchor).lower().startswith('cart') and new:
        _anchor_struts(new, s)
    _last['struts'] = dict(selection=selection, length=length, loop=loop, radius=radius,
                           color=color if isinstance(color, str) else '',
                           anchor='cartoon' if str(anchor).lower().startswith('cart') else 'atom')
    old_labels = _struts_registry.get(name, {}).get('labels', False)
    _struts_registry[name] = dict(struts=previous + new, radius=radius,
                                  color=_rgb_or_none(color, s), labels=old_labels)
    _draw_struts(name, s)
    _msg(['created %d struts for "%s" (max length %.3g, loop %.3g), object "%s"%s'
          % (len(new), selection, length, loop, name,
             ', colored by chain' if _struts_registry[name]['color'] is None else '')])
    return len(new)


def _one_atom(sel, s):
    """(model, xyz, label) for a selection of one atom or one residue (its CA)."""
    m = s.get_model(sel)
    atoms = m.atom
    if len(atoms) != 1:
        ca = [a for a in atoms if a.name == 'CA']
        residues = {(a.model, a.chain, a.resi) for a in atoms}
        if len(ca) == 1 and len(residues) == 1:
            atoms = ca
        else:
            raise ValueError('"%s" must be one atom or one residue (has %d atoms)'
                             % (sel, len(m.atom)))
    a = atoms[0]
    return a.model, np.array(a.coord, float), '%s:%s%s' % (a.chain, a.resn, a.resi)


def strut_add(atom1, atom2, name='struts', radius=-1.0, color='', anchor='cartoon', _self=cmd):
    """
DESCRIPTION

    Add one strut between two atoms (or two residues: their CA atoms).

USAGE

    strut_add atom1, atom2 [, name [, radius [, color [, anchor ]]]]

    anchor = cartoon (default): a CA shown as cartoon is joined at the
             cartoon's center line; atom: exactly at the atom

EXAMPLE

    strut_add /DN3//A/45/CA, /DN3//A/112/CA
    strut_add chain A and resi 45, chain A and resi 112
    """
    s = _self
    _sync_struts(s)
    try:
        m1, p1, l1 = _one_atom(atom1, s)
        m2, p2, l2 = _one_atom(atom2, s)
    except ValueError as e:
        print(' strut_add: %s' % e)
        return
    entry = _struts_registry.get(name)
    if entry is None:
        entry = _struts_registry[name] = dict(struts=[], radius=0.5, color=None, labels=False)
    if float(radius) > 0:
        entry['radius'] = float(radius)
    if str(color).strip():
        entry['color'] = _rgb_or_none(color, s)
    st = dict(p1=p1, p2=p2, m1=m1, m2=m2, l1=l1, l2=l2)
    if str(anchor).lower().startswith('cart'):
        _anchor_struts([st], s)
    entry['struts'].append(st)
    _draw_struts(name, s)
    _msg(['strut S%d added %s to %s (%.1f A), object "%s"'
          % (len(entry['struts']), l1, l2, np.linalg.norm(p1 - p2), name)])
    return True


def strut_remove(atom1, atom2='', name='', _self=cmd):
    """
DESCRIPTION

    Remove the strut between atom1 and atom2, or (atom2 blank) every strut
    touching atom1. Atoms can be given as one atom or one residue (its CA).

USAGE

    strut_remove atom1 [, atom2 [, name ]]
    """
    s = _self
    _sync_struts(s)
    try:
        _, p1, l1 = _one_atom(atom1, s)
        p2 = _one_atom(atom2, s)[1] if str(atom2).strip() else None
    except ValueError as e:
        print(' strut_remove: %s' % e)
        return
    def near(a, b):
        return float(np.linalg.norm(np.asarray(a, float) - np.asarray(b, float))) < 0.05
    removed = 0
    for nm in ([name] if name else list(_struts_registry)):
        entry = _struts_registry.get(nm)
        if not entry:
            continue
        keep = []
        for st in entry['struts']:
            hit1 = near(st['p1'], p1) or near(st['p2'], p1)
            if p2 is None:
                hit = hit1
            else:
                hit = (near(st['p1'], p1) and near(st['p2'], p2)) or \
                      (near(st['p1'], p2) and near(st['p2'], p1))
            if hit:
                removed += 1
            else:
                keep.append(st)
        if len(keep) != len(entry['struts']):
            entry['struts'] = keep
            _draw_struts(nm, s)
    _msg(['removed %d strut(s) at %s' % (removed, l1)])
    return removed


def strut_remove_selected(selection='sele', _self=cmd):
    """
DESCRIPTION

    Remove the struts that are in a selection. Click a strut in the viewer
    (it selects that strut into (sele)), then run this.

USAGE

    strut_remove_selected [selection]
    """
    s = _self
    from pymol import stored
    removed = []
    for nm, entry in list(_struts_registry.items()):
        if nm not in s.get_names('objects'):
            continue
        stored._p3d = []
        try:
            s.iterate('(%s) and %s' % (selection, nm), 'stored._p3d.append(int(resv))')
        except Exception:
            continue
        have = set()
        s.iterate(nm, 'H.add((int(resv), name))', space={'H': have, 'int': int})
        gone = [k for k in range(1, len(entry['struts']) + 1)
                if (k, 'E1') not in have or (k, 'E2') not in have]   # deleted in the viewer
        idx = sorted(set(stored._p3d) | set(gone), reverse=True)
        for k in idx:
            if 1 <= k <= len(entry['struts']):
                removed.append(_strut_label(k, entry['struts'][k - 1]))
                del entry['struts'][k - 1]
        if idx:
            _draw_struts(nm, s)
    _msg(['removed %d strut(s)%s' % (len(removed), ': ' + ', '.join(removed) if removed else
                                    ' (click a strut in the viewer to select it first)')])
    return len(removed)


def struts_label(on=1, name='', _self=cmd):
    """
DESCRIPTION

    Show (on=1) or hide (on=0) an ID label on every strut, e.g.
    "S3 A:LYS45-A:GLU112". Clicking a strut also prints its ID in the log
    (S3 is residue STR 3 of the struts object).

USAGE

    struts_label [on [, name ]]
    """
    s = _self
    _sync_struts(s)
    for nm in ([name] if name else list(_struts_registry)):
        if nm in _struts_registry:
            _struts_registry[nm]['labels'] = _b(on)
            _draw_struts(nm, s)
    _msg(['strut labels %s' % ('on' if _b(on) else 'off')])


def struts_delete(name='struts', _self=cmd):
    """
DESCRIPTION

    Delete a whole struts object made by "struts".
    """
    _struts_registry.pop(name, None)
    _self.delete(name)


# ---------------------------------------------------------------------------
# ligand surface
# ---------------------------------------------------------------------------

LIGAND_DEFAULT = 'organic and not solvent'


def ligand_surface(selection=LIGAND_DEFAULT, name='ligand_surf', color='green',
                   quality=1, hide_source=1, _self=cmd):
    """
DESCRIPTION

    Copy selection (default: all non-polymer organic molecules, e.g. a
    substrate, cofactor or inhibitor)
    into its own object and show its surface alone, so it is not cut into by
    the protein surface.

USAGE

    ligand_surface [selection [, name [, color [, quality [, hide_source ]]]]]
    """
    s = _self
    if s.count_atoms(selection) == 0:
        print(' ligand_surface: no atoms in "%s"' % selection)
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
    ligand covalently bonded to a residue is still separated from the protein)."""
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
        for st in entry['struts']:
            st['m1'] = lookup.get(key(st['p1']), st['m1'])
            st['m2'] = lookup.get(key(st['p2']), st['m2'])


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

    by = ligands: polymer + one object per connected ligand (e.g. each substrate)
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
# printers, surfaces and built-in presets
# ---------------------------------------------------------------------------

# Build volumes (X, Y, Z in mm) from the manufacturers' specifications;
# slots = filaments the stock multicolor unit holds (1 = single color).
PRINTERS = {
    'Anycubic Kobra X': dict(bed=(260, 260, 260), slots=4),
    'Anycubic Kobra 3 / Kobra 3 Combo': dict(bed=(250, 250, 260), slots=4),
    'Anycubic Kobra S1 / S1 Combo': dict(bed=(250, 250, 250), slots=4),
    'Anycubic Kobra S1 Max': dict(bed=(350, 350, 350), slots=4),
    'Anycubic Kobra 3 Max': dict(bed=(420, 420, 500), slots=4),
    'Bambu Lab A1 mini': dict(bed=(180, 180, 180), slots=4),
    'Bambu Lab A1': dict(bed=(256, 256, 256), slots=4),
    'Bambu Lab A2L': dict(bed=(330, 320, 325), slots=4),
    'Bambu Lab P1S / P1P': dict(bed=(256, 256, 256), slots=4),
    'Bambu Lab X1 Carbon': dict(bed=(256, 256, 256), slots=4),
    'Prusa MINI+': dict(bed=(180, 180, 180), slots=1),
    'Prusa MK4S': dict(bed=(250, 210, 220), slots=1),
    'Prusa CORE One': dict(bed=(250, 220, 270), slots=1),
    'Creality K1C': dict(bed=(220, 220, 250), slots=1),
    'Creality K2 Plus': dict(bed=(350, 350, 350), slots=4),
}

SURFACE_QUALITY = [(-1, 'Draft (fast, coarse)'), (0, 'Normal'), (1, 'Fine'),
                   (2, 'Extra fine'), (3, 'Ultra fine'), (4, 'Maximum (slow, big files)')]


def surfacestyle(selection='all', quality=1, probe=1.4, solvent=0, smooth_edges=1,
                 color='', show=1, _self=cmd):
    """
DESCRIPTION

    Set how surfaces are built, for printing. Applies to the objects in
    selection, and shows their surface (show=1).

USAGE

    surfacestyle [selection [, quality [, probe [, solvent [, smooth_edges
                 [, color [, show ]]]]]]]

ARGUMENTS

    quality = mesh resolution: -1 draft, 0 normal, 1 fine, 2 extra fine,
              3 ultra fine, 4 maximum (finer = smoother print, bigger file)
    probe = probe radius in A (1.4 = water; larger smooths out small pockets)
    solvent = 0: molecular surface (follows the atoms), 1: solvent accessible
              surface (inflated by the probe radius, chunkier)
    smooth_edges = 1: smooth the edges where surface patches meet
    color = optional surface color
    """
    s = _self
    objs = [o for o in s.get_object_list('(%s)' % selection) if o not in _struts_registry]
    if not objs:
        print(' surfacestyle: no objects in "%s"' % selection)
        return
    for o in objs:
        s.set('surface_quality', int(float(quality)), o)
        s.set('solvent_radius', float(probe), o)
        s.set('surface_solvent', int(_b(solvent)), o)
        s.set('surface_smooth_edges', int(_b(smooth_edges)), o)
        s.set('surface_type', 0, o)          # solid triangles (needed for printing)
    if _b(show):
        s.show('surface', '(%s) and not (%s)' % (selection, ' or '.join(_struts_registry) or 'none'))
    if str(color).strip():
        s.set('surface_color', color, ' or '.join(objs))
    _last['surface'] = dict(selection=selection, quality=int(float(quality)), probe=float(probe),
                            solvent=int(_b(solvent)), smooth_edges=int(_b(smooth_edges)),
                            color=color)
    label = dict(SURFACE_QUALITY).get(int(float(quality)), str(quality))
    _msg(['surface: %s, probe %.2f A, %s, on %s' % (
        label, float(probe), 'solvent accessible' if _b(solvent) else 'molecular', ', '.join(objs))])


_CARTOON_STYLE = dict(helix_width=3.0, helix_thickness=1.5, strand_width=3.0,
                      strand_thickness=1.0, loop_thickness=1.25, stick_radius=0.7,
                      sampling=14, hide_h=1)
_STRUTS = dict(selection='name CA', length=8.0, loop=60.0, radius=0.5, color='', name='struts')
_SURFACE = dict(selection='all', quality=2, probe=1.4, solvent=0, smooth_edges=1, color='')

BUILTIN_PRESETS = {
    'Multi-color Printer Cartoon': dict(
        style=dict(_CARTOON_STYLE), struts=dict(_STRUTS),
        export=dict(format='3mf', single_color=0, base='stand', attach=1, union=1)),
    'Single Color Printer Cartoon': dict(
        style=dict(_CARTOON_STYLE), struts=dict(_STRUTS),
        export=dict(format='stl', single_color=1, base='stand', attach=1, union=1)),
    'Multi-color Printer Surface': dict(
        surface=dict(_SURFACE),
        export=dict(format='3mf', single_color=0, base='cut', base_depth=1.0, attach=1, union=1)),
    'Single Color Printer Surface': dict(
        surface=dict(_SURFACE),
        export=dict(format='stl', single_color=1, base='cut', base_depth=1.0, attach=1, union=1)),
}


def _printer_overlay_cgo(bed, pl, zshift, s=cmd):
    """Wire box of the build volume (+ a faint bed plate), in molecule space."""
    bx, by, bz = [float(v) for v in bed]
    c = np.array([[x, y, z] for z in (0, bz) for y in (-by / 2, by / 2) for x in (-bx / 2, bx / 2)])
    A = pl.inverse(c - np.array([0.0, 0.0, zshift]))
    edges = [(0, 1), (2, 3), (0, 2), (1, 3), (4, 5), (6, 7), (4, 6), (5, 7),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    r = 0.6 / pl.scale                       # 0.6 mm thick edges
    obj = []
    for i, j in edges:
        obj += [cgo.CYLINDER] + list(A[i]) + list(A[j]) + [r, 0.2, 0.5, 1.0, 0.2, 0.5, 1.0]
    P = pl.inverse(c[:4] - np.array([0.0, 0.0, zshift + 0.5]))   # bed, 0.5 mm below the print
    obj += [cgo.ALPHA, 0.2, cgo.BEGIN, cgo.TRIANGLES, cgo.COLOR, 0.6, 0.75, 1.0]
    for t in ((0, 1, 3), (0, 3, 2)):
        for k in t:
            obj += [cgo.VERTEX] + list(P[k])
    obj += [cgo.END, cgo.ALPHA, 1.0]
    return obj


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

_connectors = []          # list of dict(host, guest, kind, ...)
_last = {}                # last-used options per section (for presets)

EXPORT_DEFAULTS = dict(objects='', format='3mf', each=0, size_mm=0.0, mm_per_A=1.0,
                       nozzle=0.4, layer_height=0.2, preview=0, printer='', bed_x=0.0,
                       bed_y=0.0, bed_z=0.0, slots=0, fit=0, show_bed=1, show_thin=1,
                       helix_core=0, support_angle=30.0, show_overhang=0, join_strands=1,
                       use_manifold=1,
                       orient='view', filaments='', base='none', base_depth=1.0,
                       pillars=3, pillar_d=0.0, min_feature=0.0, attach=1, attach_mm=4.0,
                       single_color=0, socket_clearance=0.0,
                       union=1, close_ends=1,
                       round_ends=1, segments=0, color_tolerance=0)


def _cartoon_rgb(s, obj, ca_rgb):
    try:
        cc = s.get('cartoon_color', obj)
        if cc not in (None, '', 'default', '-1'):
            return np.array([int(round(255 * x)) for x in s.get_color_tuple(cc)], np.uint8)
    except Exception:
        pass
    return ca_rgb


def _helix_cores(s, names, segments=12, overlap=0.4):
    """Solid rods along the axis of every helix shown as cartoon, sized to fuse
    with the ribbon, so no support can build up inside the helix.
    Returns mesh pieces in Angstrom and the number of helices."""
    pieces, nhel = [], 0
    for obj in names:
        if obj in _struts_registry or _is_overlay(obj):
            continue
        try:
            if int(s.get('cartoon_cylindrical_helices', obj)):
                continue                                  # already solid
            thick = float(s.get('cartoon_oval_width', obj)) * _HELIX_FACTOR
        except Exception:
            continue
        rows = []
        s.iterate_state(-1, '(%s) and name CA and ss H and rep cartoon' % obj,
                        'rows.append((segi, chain, int(resv), x, y, z, color))',
                        space={'rows': rows, 'int': int})
        runs, cur = [], []
        for r in sorted(rows, key=lambda r: (r[0], r[1], r[2])):
            if cur and (r[0], r[1]) == cur[-1][:2] and r[2] == cur[-1][2] + 1:
                cur.append(r)
            else:
                if cur:
                    runs.append(cur)
                cur = [r]
        if cur:
            runs.append(cur)
        for run in runs:
            if len(run) < 5:
                continue
            ca = np.array([r[3:6] for r in run], float)
            ax = []
            for i in range(1, len(ca) - 1):
                b = (ca[i - 1] - ca[i]) + (ca[i + 1] - ca[i])
                nb = np.linalg.norm(b)
                if nb > 1e-6:
                    ax.append(ca[i] + 2.3 * b / nb)
            if len(ax) < 3:
                continue
            ax = np.array(ax)
            sm = ax.copy()
            sm[1:-1] = (ax[:-2] + ax[1:-1] + ax[2:]) / 3.0
            rad = np.median([np.min(np.linalg.norm(sm - c, axis=1)) for c in ca])
            r_core = max(0.5, rad - thick / 2.0 + overlap)
            d0 = sm[0] - sm[1]
            d1 = sm[-1] - sm[-2]
            pts = np.vstack([sm[0] + d0 / max(np.linalg.norm(d0), 1e-6) * 1.0, sm,
                             sm[-1] + d1 / max(np.linalg.norm(d1), 1e-6) * 1.0])
            col = _cartoon_rgb(s, obj, np.array([int(round(255 * x)) for x in
                                                 s.get_color_tuple(run[len(run) // 2][6])], np.uint8))
            for a, b in zip(pts[:-1], pts[1:]):
                va, ta = pc.capsule_mesh(a, b, r_core, segments, max(3, segments // 4), True)
                pieces.append((va, ta, np.tile(col, (len(va), 1))))
            nhel += 1
    return pieces, nhel


def _catmull(pts, sub=4):
    """Smooth path through pts (Catmull-Rom), sub points per span."""
    P = np.asarray(pts, float)
    if len(P) < 2:
        return P
    Q = np.vstack([2 * P[0] - P[1], P, 2 * P[-1] - P[-2]])
    out = []
    for i in range(1, len(Q) - 2):
        p0, p1, p2, p3 = Q[i - 1], Q[i], Q[i + 1], Q[i + 2]
        for t in np.linspace(0, 1, sub, endpoint=False):
            t2, t3 = t * t, t * t * t
            out.append(0.5 * (2 * p1 + (p2 - p0) * t + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t2 +
                              (3 * p1 - p0 - 3 * p2 + p3) * t3))
    out.append(P[-1])
    return np.array(out)


def _center_on_mesh(V, pts, half=0.35, reach=2.6, max_shift=1.0):
    """Move each point to the middle of the cartoon cross-section around it."""
    P = np.asarray(pts, float)
    T = np.gradient(P, axis=0)
    T /= np.maximum(np.linalg.norm(T, axis=1), 1e-12)[:, None]
    out = P.copy()
    near = V[np.linalg.norm(V - P.mean(0), axis=1) < np.linalg.norm(np.ptp(P, axis=0)) + reach + 2]
    for k, (p_, t) in enumerate(zip(P, T)):
        d = near - p_
        along = d @ t
        lat = np.linalg.norm(d - np.outer(along, t), axis=1)
        sel = (np.abs(along) < half) & (lat < reach)
        if sel.sum() >= 6:
            c = near[sel].mean(0)
            c -= t * ((c - p_) @ t)              # stay in this cross-section
            if np.linalg.norm(c - p_) < max_shift:
                out[k] = c
    sm = out.copy()                               # light smoothing
    sm[1:-1] = 0.25 * out[:-2] + 0.5 * out[1:-1] + 0.25 * out[2:]
    return sm


def _strand_joins(s, names, segments=24):
    """PyMOL ends each beta strand in an arrow whose point is the only place
    the next loop touches it, which prints as a weak point. Continue the loop
    into the strand the way ChimeraX does: one tube swept along the cartoon's
    own center line, round like the loop on the loop side and flattening to
    the strand's thickness (flush with its faces) inside the strand, so the
    arrow ends in a blunt tip the width of the loop that flows into the loop.
    Returns (pieces, n)."""
    pieces, n = [], 0
    objs = [o for o in names if o not in _struts_registry and not _is_overlay(o)]
    meshes = {}
    amap = _cartoon_anchor_map(objs, s, meshes)
    for obj in objs:
        V = meshes.get(obj)
        if V is None:
            continue
        try:
            loop_r = float(s.get('cartoon_loop_radius', obj))
            half_t = float(s.get('cartoon_rect_width', obj)) * _STRAND_FACTOR / 2.0
        except Exception:
            continue
        r_loop = 0.97 * loop_r                  # just inside the loop tube
        r_flat = 0.96 * min(half_t, loop_r)     # just inside the strand faces
        rows = []
        s.iterate_state(-1, '(%s) and name CA and rep cartoon' % obj,
                        'rows.append((segi, chain, int(resv), ss, color, x, y, z))',
                        space={'rows': rows, 'int': int})
        rows.sort(key=lambda r: (r[0], r[1], r[2]))
        runs, cur = [], []
        for r in rows:
            if cur and r[:2] == cur[-1][:2] and r[2] == cur[-1][2] + 1:
                cur.append(r)
            else:
                if cur:
                    runs.append(cur)
                cur = [r]
        if cur:
            runs.append(cur)
        for run in runs:
            ctr = [amap.get((obj, _pkey(r[5:8])), np.array(r[5:8], float)) for r in run]
            isS = [r[3] == 'S' for r in run]
            for i in range(len(run) - 1):
                if isS[i] == isS[i + 1]:
                    continue
                si = i if isS[i] else i + 1             # strand residue at the boundary
                inward = -1 if isS[i] else 1            # direction into the strand
                s2 = si + inward
                if not (0 <= s2 < len(run)) or not isS[s2]:
                    s2 = si
                # strand face normal: flattest direction of the cartoon around it
                cen = 0.5 * (ctr[si] + ctr[s2])
                near = V[np.linalg.norm(V - cen, axis=1) < 2.4]
                if len(near) < 6:
                    continue
                w_, v_ = np.linalg.eigh(np.cov((near - near.mean(0)).T))
                up = v_[:, 0]
                idx = list(range(max(i - 1, 0), min(i + 3, len(run))))
                sub = 10
                path = _catmull([ctr[j] for j in idx], sub)
                pos = np.arange(len(path)) / float(sub) + idx[0]
                # x: 0 at the strand's last residue, 1 at the first loop residue
                x = (pos - i) if isS[i] else (i + 1 - pos)
                keep = (x >= -1.0) & (x <= 1.3)
                path, pos, x = path[keep], pos[keep], x[keep]
                if len(path) < 4:
                    continue
                path = _center_on_mesh(V, path)
                # flat (strand-thick) until the arrow tip, then round like the loop
                w = np.clip((x - 0.35) / 0.55, 0, 1)
                w = w * w * (3 - 2 * w)
                thick = r_flat + (r_loop - r_flat) * w
                # bury the far end inside the loop tube (no near-coincident skin)
                e = np.clip((x - 0.85) / 0.45, 0, 1)
                bury = 1.0 - 0.25 * e * e * (3 - 2 * e)
                va, ta, ring = pc.sweep_tube(path, r_loop * bury, segments,
                                             radii2=np.minimum(thick, r_loop * bury), up=up)
                res = np.clip(np.round(pos).astype(int), 0, len(run) - 1)
                cols = np.array([_cartoon_rgb(s, obj, np.array(
                    [int(round(255 * c)) for c in s.get_color_tuple(run[k][4])], np.uint8))
                    for k in res], np.uint8)
                pieces.append((va, ta, cols[ring]))
                n += 1
    return pieces, n


def _is_overlay(name):
    """Display-only objects that are never exported."""
    return name.startswith('preview_') or name == JOINS_OBJECT


JOINS_OBJECT = 'strand_joins'


def show_strand_joins(objects='', on=1, _self=cmd):
    """
DESCRIPTION

    Show (on=1) in PyMOL the tubes that colorprint adds where loops meet beta
    strands (join_strands=1), as the display-only object "strand_joins"
    (never exported; colorprint builds its own). on=0 removes it.
    Redrawn by printstyle (show_joins=1) and after each colorprint.

USAGE

    show_strand_joins [objects [, on ]]
    """
    s = _self
    s.delete(JOINS_OBJECT)
    if not _b(on):
        return 0
    names = [n for n in str(objects).replace(',', ' ').split() if n] or \
        [n for n in s.get_names('objects', enabled_only=1)
         if not _is_overlay(n) and n not in _struts_registry]
    pieces, nj = _strand_joins(s, names, 24)
    if not pieces:
        return 0
    cg = []
    for va, ta, vc in pieces:
        cg += pc.mesh_cgo(va, ta, np.asarray(vc)[np.asarray(ta)[:, 0]])
    s.load_cgo(cg, JOINS_OBJECT)
    return nj


def _single_bonds(s, objects=None):
    """Draw double/aromatic bonds as one plain stick (valence off), globally and
    on every object (or just `objects`). Returns what to pass to _restore_bonds."""
    objs = s.get_names('objects') if objects is None else objects
    old = {'': s.get('valence')}
    for o in objs:
        try:
            old[o] = s.get('valence', o)
            s.set('valence', 0, o)
        except Exception:
            pass
    s.set('valence', 0)
    return old


def _restore_bonds(s, old):
    for o, v in old.items():
        try:
            if o:
                s.set('valence', v, o)
        except Exception:
            pass
    s.set('valence', old.get('', 'on'))


def _scene_vrml(s):
    """VRML of everything shown, in model coordinates (independent of view).
    Bond orders are flattened so double/aromatic bonds export as one full stick."""
    bonds = _single_bonds(s, s.get_names('objects', enabled_only=1))
    try:
        return _scene_vrml_raw(s)
    finally:
        _restore_bonds(s, bonds)


def _scene_vrml_raw(s):
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
    if not names:       # "everything shown", but never the preview objects
        names = [n for n in s.get_names('objects', enabled_only=1)
                 if not _is_overlay(n)]
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
               nozzle=-1.0, layer_height=-1.0, preview=-1, printer='', bed_x=-1.0,
               bed_y=-1.0, bed_z=-1.0, slots=-1, fit=-1, show_bed=-1, show_thin=-1,
               helix_core=-1, support_angle=-1.0, show_overhang=-1, join_strands=-1,
               use_manifold=-1,
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
    orient = view: printed the way you see it (screen-down = bed) | model |
             min_support: turned to need the least support material
    filaments = colors loaded in slots 1, 2, 3... e.g. "cyan orange blue green";
                every color is snapped to the closest filament.
                filaments=auto picks common PLA colors closest to the scene
                (as many as the printer's filament slots, default 4)
    use_manifold = 0: don't use the merge library (manifold3d) even if it is
             installed: plain capped meshes, no merging, bases or connectors
    join_strands = 1 (default): run each loop into its beta strand with a rod
             along the cartoon's center line, so arrow tips aren't the only
             contact (they print as a weak point)
    helix_core = 1: fill each cartoon helix with a solid rod along its axis, so
             no support builds up inside it
    support_angle = overhangs flatter than this (degrees from horizontal) need
             support, as in the slicer (default 30)
    show_overhang = 1: in preview, paint faces that will need support magenta
    base = none | cut (flat bottom, base_depth mm) | stand (plate + pillars)
    nozzle, layer_height = what to optimize for (default 0.4 / 0.2 mm). They
             set the defaults below whenever those are 0 (= auto): thinnest
             feature 3.75 x nozzle, fit clearance nozzle/2, socket clearance
             0.4 x nozzle, pillar thickness 7.5 x nozzle (>= 3 mm), base sizes in
             whole layers and rod smoothness within half a layer
    printer = a profile name from print3d_printers (sets the build volume and
             filament slots), or give bed_x, bed_y, bed_z (mm) and slots yourself
    fit = 1: scale the model to fill the build volume (overrides size_mm)
    show_bed, show_thin = in preview: draw the build volume, and highlight in
             red any rods thinner than the printer can make (default on)
    preview = 1: build everything but load it into PyMOL as objects named
             preview_<name> (in place over the molecule) instead of writing files
    min_feature = warn below this thickness in mm (0 = auto)
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
        opts.update(all_presets().get(str(preset), {}).get('export', {}))
    given = dict(objects=objects, format=format, each=each, size_mm=size_mm, mm_per_A=mm_per_A,
                 nozzle=nozzle, layer_height=layer_height, preview=preview, printer=printer,
                 bed_x=bed_x, bed_y=bed_y, bed_z=bed_z, slots=slots, fit=fit,
                 show_bed=show_bed, show_thin=show_thin, helix_core=helix_core,
                 support_angle=support_angle, show_overhang=show_overhang,
                 join_strands=join_strands, use_manifold=use_manifold,
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
    preview = _b(o['preview'])

    # ---- printer: nozzle + layer height set the auto values ----
    nz = max(float(o['nozzle']), 0.05)
    lh = max(float(o['layer_height']), 0.02)

    def layers(x, min_layers=1):            # round up to whole layers
        import math
        return max(min_layers, math.ceil(float(x) / lh - 1e-6)) * lh
    min_feature = float(o['min_feature']) or round(3.75 * nz, 2)
    socket_c = float(o['socket_clearance']) or round(0.4 * nz, 2)
    fit_c = round(0.5 * nz, 2)
    pillar_d = float(o['pillar_d']) or max(3.0, round(7.5 * nz, 1))
    base_depth = layers(o['base_depth'], 2)
    plate_t = layers(max(2.0, 8 * lh), 4)
    stand_gap = layers(4.0)

    prof = PRINTERS.get(str(o['printer']).strip()) if str(o['printer']).strip() else None
    if str(o['printer']).strip() and prof is None:
        print(' colorprint: unknown printer "%s" (see print3d_printers)' % o['printer'])
    bed = None
    if float(o['bed_x']) > 0 and float(o['bed_y']) > 0 and float(o['bed_z']) > 0:
        bed = (float(o['bed_x']), float(o['bed_y']), float(o['bed_z']))
    elif prof:
        bed = prof['bed']
    slots = int(float(o['slots'])) or (prof['slots'] if prof else 0)
    fit = _b(o['fit']) and bed is not None
    names = str(o['objects']).replace(',', ' ').split()
    missing = [n for n in names if n not in s.get_names('all')]
    if missing:
        print(' colorprint: unknown object(s): %s' % ', '.join(missing))
        return
    _last['export'] = {k: o[k] for k in EXPORT_DEFAULTS}
    solid_ok = pc.have_manifold() and _b(o['use_manifold'])

    # struts deleted in the viewer stay deleted; struts colored by chain
    # follow any recoloring done since they were made
    _sync_struts(s)
    for sname, e in list(_struts_registry.items()):
        if e.get('color') is None and sname in s.get_names('objects'):
            _draw_struts(sname, s)

    # ---- what to export: one group, or one per object (struts follow) ----
    groups = []                     # (label, names, temp_cgo)
    if each:
        targets = [n for n in (names or s.get_names('objects', enabled_only=1))
                   if not _is_overlay(n)]
        strut_objs = [n for n in targets if n in _struts_registry] if names else \
            [n for n in _struts_registry if n in targets]
        for obj in [n for n in targets if n not in _struts_registry]:
            cgo_list = []
            nsub = 0
            for sname in strut_objs:
                e = _struts_registry[sname]
                sub = [st for st in e['struts'] if st['m1'] == obj and st['m2'] == obj]
                nsub += len(sub)
                cgo_list += _strut_cgo(sub, e['radius'], e['color'], s)
            temp = ('_p3d_struts_tmp', cgo_list) if cgo_list else None
            groups.append((obj, [obj], temp))
    else:
        label = names[0] if len(names) == 1 else os.path.basename(stem)
        groups.append((label, names, None))

    notes_all = []
    vrmls = [(label, gnames, _vrml_of(s, gnames, temp)) for label, gnames, temp in groups]

    def tessellate(seg):
        out = []
        for label, gnames, vr in vrmls:
            notes, info = [], {}
            pieces = pc.vrml_to_pieces(vr, seg, _b(o['round_ends']), _b(o['close_ends']),
                                       notes, info)
            gn = gnames or [n for n in s.get_names('objects', enabled_only=1)
                            if not _is_overlay(n)]
            if pieces and _b(o['join_strands']):
                joins, nj = _strand_joins(s, gn, max(16, 2 * seg))
                pieces += joins
                if nj:
                    notes.append('joined %d strand ends to their loops' % nj)
            if pieces and _b(o['helix_core']):
                cores, nh = _helix_cores(s, gn, max(8, seg))
                pieces += cores
                if nh:
                    notes.append('filled %d helices with a solid core' % nh)
            if pieces:
                out.append(dict(label=label, names=gn,
                                pieces=pieces, radii=info.get('radii', []),
                                prims=info.get('prims', []), notes=notes))
        return out
    raw = tessellate(segments or 12)
    if not raw:
        print(' colorprint: nothing shown to export')
        return

    # ---- shared placement so separately printed parts still match ----
    omode = str(o['orient']).lower()
    sangle = float(o['support_angle'])
    rot = pc.view_rotation(s.get_view()) if omode in ('view', 'min_support', 'support') else None
    if omode in ('min_support', 'support'):
        allp = [pc_ for g in raw for pc_ in g['pieces']]
        view_up = rot[:, 2]
        up, area_best, _, _ = pc.best_up(allp, sangle, prefer=view_up)
        area_view = pc.support_scores(allp, [view_up], sangle)[0]
        rot = pc.frame_with_up(up, rot[:, 0])
        notes_all.append('turned for the least support: overhangs %.0f A2 (as on screen: %.0f A2)'
                         % (area_best, area_view))
    allpts = np.vstack([va for g in raw for va, _, _ in g['pieces']])
    if fit:
        Rm = np.eye(3) if rot is None else rot
        ext = np.ptp(allpts @ Rm, axis=0)
        room_z = bed[2] - (stand_gap + plate_t if str(o['base']).lower() == 'stand' else 0)
        sc = min(0.95 * bed[0] / ext[0], 0.95 * bed[1] / ext[1], 0.97 * room_z / ext[2])
        pl = pc.make_placement(allpts, rot, sc, 0)
    else:
        pl = pc.make_placement(allpts, rot, mm_per_A, size_mm)
    if not segments:      # auto: facets within half a layer of the true rod surface
        radii = [r for g in raw for _, r in g['radii']]
        if radii:
            segments = pc.auto_segments(float(np.median(radii)) * pl.scale, lh / 2)
            if segments != 12:
                raw = tessellate(segments)

    # ---- filament palette ----
    picked = []
    if str(o['filaments']).strip().lower() == 'auto':
        picked = pc.pick_filaments([pc_ for g in raw for pc_ in g['pieces']], slots or 4)
        o = dict(o, filaments=' '.join(h for _, h, _ in picked))
        notes_all.append('filaments picked from the scene: ' + ', '.join(
            '%d=%s (%s)' % (i + 1, n, h) for i, (n, h, _) in enumerate(picked)))
    pal = _parse_colors(o['filaments'], s) if o['filaments'] else []
    if picked and pal:
        pal = [(picked[i][0] if i < len(picked) else nm, c) for i, (nm, c) in enumerate(pal)]
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
                   float(con['length']), float(con['clearance']) or fit_c, float(con['wall']),
                   float(con['magnet_d']), float(con['magnet_h']), notes=gh['notes'])

    # ---- join loose pieces (side chains just off the cartoon, etc.) ----
    if solid_ok and _b(o['attach']):
        for g in raw:
            pc.attach_floating(g['parts'], float(o['attach_mm']),
                               max(min_feature, 2.5 * nz), notes=g['notes'])

    # ---- base on the biggest object ----
    main = max(raw, key=lambda g: sum(len(p.get_mesh()[1]) for p in g['parts']))
    bmode = str(o['base']).lower()
    if bmode in ('cut', 'stand') and not solid_ok:
        notes_all.append('base needs manifold3d')
    for g in raw:
        g['zshift'] = 0.0
    if bmode == 'cut' and solid_ok:
        main['parts'], main['zshift'] = pc.flat_cut(main['parts'], base_depth, main['notes'])
    elif bmode == 'stand' and solid_ok:
        stand, main['zshift'] = pc.pillar_stand(main['parts'], int(float(o['pillars'])), pillar_d,
                                                plate_t=plate_t, gap=stand_gap, notes=main['notes'])
        if stand is not None:
            c = max(main['parts'], key=lambda p: p.volume()).color
            main['parts'].append(pc.Part(c, solid=pc.tint(stand, c), name='stand'))

    single = _b(o['single_color'])
    if single:
        if solid_ok:
            for g in raw:
                g['parts'] = pc.separate_colors(g['parts'], socket_c, g['notes'])
        else:
            notes_all.append('single-color pieces without sockets (manifold3d not installed)')

    # ---- report + write ----
    out_files, report = [], []
    report.append('scale %.2f mm per A, %s' % (pl.scale, 'oriented as on screen'
                                                if rot is not None else 'molecule frame'))
    if bed:
        report.append('printer: %s, build volume %d x %d x %d mm%s%s' % (
            o['printer'] or 'custom', bed[0], bed[1], bed[2],
            ', %d filament slots' % slots if slots else '', ', scaled to fit' if fit else ''))
    report.append('tuned for %.2f mm nozzle, %.2f mm layers: thinnest feature %.2f mm, fit '
                  'clearance %.2f mm, socket clearance %.2f mm, %d facets per rod'
                  % (nz, lh, min_feature, fit_c, socket_c, segments))
    if pal:
        report.append('filaments: ' + ', '.join('%d=%s' % (i + 1, n) for i, (n, _) in enumerate(pal)))
        report.append('   colors mapped: ' + ', '.join(
            '#%s->%d' % (pc.hex_color(c), sl + 1) for c, sl in orig_map))
    thin_cgo, thin_tris, over_cgo = [], 0, []
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
        if bed and a.get('size') is not None:
            sz = np.asarray(a['size'], float)
            over = [('XYZ'[i], sz[i], bed[i]) for i in range(3) if sz[i] > bed[i] + 1e-6]
            if not over:
                report.append('   fits the build volume')
            elif all(np.sort(sz) <= np.sort(bed) + 1e-6):
                report.append('   too big as oriented (%s), but fits if rotated in the slicer'
                              % ', '.join('%s %.0f > %.0f mm' % t for t in over))
            else:
                report.append('   DOES NOT FIT: %s; lower the size, use fit=1, or split it'
                              % ', '.join('%s %.0f > %.0f mm' % t for t in over))
        if slots and not single:
            used = {tuple(int(round(float(x) * (255 if max(p.color) <= 1 else 1))) for x in p.color[:3])
                    for p in parts}
            if len(used) > slots:
                report.append('   uses %d colors but the printer has %d filament slot%s: set '
                              'the filaments to snap colors, recolor, or use single_color=1'
                              % (len(used), slots, '' if slots == 1 else 's'))
        per = []
        for p in parts:
            if p.solid is not None:
                per.append('%s%s on filament %d: %.1f g' % (
                    _color_name(p.color, pal_names), ' (%s)' % p.name if p.name else '',
                    slot_of(p.color), p.volume() / 1000 * pc.PLA_DENSITY))
        if per:
            report.append('   ' + '; '.join(per))
        ov = 0.0
        for prt in parts:
            va_, ta_, _ = pc.part_output(prt)
            if len(ta_):
                ov += pc.overhang_area(va_, ta_, sangle, 0.0, lh)
        if ov > 1.0:
            report.append('   needs support under about %.1f cm2 (overhangs flatter than %.0f deg)%s'
                          % (ov / 100.0, sangle, '' if omode in ('min_support', 'support') else
                             '; orient=min_support may need less'))
        for w in _thin_features(s, g['names'], g['radii'], pl.scale, min_feature):
            report.append('   thin: ' + w)
        bodies = a.get('bodies') or []
        if len(bodies) > 1:
            fl = bodies[1:]
            report.append('   floating: %d piece(s) not attached to the main body: %s'
                          % (len(fl), '; '.join('%.0f mm3 near %s' % (v, n) for v, _, n in fl[:6])
                                               + (' ...' if len(fl) > 6 else '')))
        if a.get('bed_contact') is not None and not (g is main and bmode in ('stand', 'cut')):
            report.append('   bed contact %.0f mm2%s' % (a['bed_contact'], ' (small: supports '
                          'or base=cut/stand recommended)' if a['bed_contact'] < 30 else ''))
        report += ['   ' + n for n in g['notes'] if not n.startswith('Read ')]
        path_stem = _file_for(stem, '', g['label'], each)
        if preview:
            cg = []
            for prt in parts:
                va, ta, tc = pc.part_output(prt)
                if _b(o['show_thin']) and prt.solid is not None:
                    tm = pc.thin_mask(prt.solid, va, ta, min_feature)
                    if tm.any():
                        thin_tris += int(tm.sum())
                        tri = va[ta[tm]]                       # (k, 3, 3), mm
                        nn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
                        nn /= np.maximum(np.linalg.norm(nn, axis=1), 1e-12)[:, None]
                        tri = (tri + 0.08 * nn[:, None, :]).reshape(-1, 3)   # lift over the part
                        tv = pl.inverse(tri - np.array([0.0, 0.0, g.get('zshift', 0.0)]))
                        thin_cgo += pc.mesh_cgo(tv, np.arange(len(tv)).reshape(-1, 3),
                                                np.tile([235, 25, 25], (int(tm.sum()), 1)))
                if _b(o['show_overhang']) and len(ta):
                    om = pc.overhang_mask(va, ta, sangle, 0.0, lh)
                    if om.any():
                        tri = va[ta[om]]
                        nn = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
                        nn /= np.maximum(np.linalg.norm(nn, axis=1), 1e-12)[:, None]
                        tri = (tri + 0.06 * nn[:, None, :]).reshape(-1, 3)
                        tv = pl.inverse(tri - np.array([0.0, 0.0, g.get('zshift', 0.0)]))
                        over_cgo += pc.mesh_cgo(tv, np.arange(len(tv)).reshape(-1, 3),
                                                np.tile([225, 40, 200], (int(om.sum()), 1)))
                va = pl.inverse(va - np.array([0.0, 0.0, g.get('zshift', 0.0)]))
                cg += pc.mesh_cgo(va, ta, tc)
            pname = 'preview_' + _safe(g['label'])
            s.delete(pname)
            s.load_cgo(cg, pname)
            out_files.append(pname)
            continue
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
    if JOINS_OBJECT in s.get_names('objects'):
        show_strand_joins(on=_b(o['join_strands']), _self=s)
    if preview:
        s.delete('preview_printer')
        s.delete('preview_thin')
        s.delete('preview_overhang')
        if over_cgo:
            s.load_cgo(over_cgo, 'preview_overhang')
            out_files.append('preview_overhang')
            report.append('magenta in preview_overhang = faces the slicer will support '
                          '(flatter than %.0f deg)' % sangle)
        if bed and _b(o['show_bed']):
            s.load_cgo(_printer_overlay_cgo(bed, pl, main.get('zshift', 0.0), s), 'preview_printer')
            out_files.append('preview_printer')
        if _b(o['show_thin']):
            if thin_cgo:
                s.load_cgo(thin_cgo, 'preview_thin')
                s.set('cgo_transparency', 0.0, 'preview_thin')
                out_files.append('preview_thin')
                report.append('red in preview_thin = walls or rods thinner than %.2f mm '
                              '(thicken them or scale up)' % min_feature)
            else:
                report.append('nothing thinner than %.2f mm' % min_feature)
        report.append('preview loaded: %s (nothing written; hide with: disable preview_*, '
                      'remove with: delete preview_*)' % ', '.join(out_files))
        _msg(report)
        _last['report'] = report
        return ', '.join(out_files)
    report.append('wrote: ' + ', '.join(os.path.basename(f) for f in out_files))
    if fmt == 'stl' and not single:
        report.append('load the STLs together as ONE object with multiple parts')
    _msg(report)
    _last['report'] = report
    return ', '.join(os.path.basename(f) for f in out_files)


# ---------------------------------------------------------------------------
# snap-fit connectors
# ---------------------------------------------------------------------------

def snapfit(host, guest, kind='peg', diameter=4.0, length=6.0, clearance=0.0,
            magnet_d=6.0, magnet_h=3.0, wall=1.6, _self=cmd):
    """
DESCRIPTION

    Join two separately printed objects (e.g. protein and ligand) at their
    contact: the guest gets a seat cut to the host's shape, plus either a peg
    (host gets a socket) or magnet pockets in both. Applied by
    "colorprint ..., each=1". Sizes in mm (after scaling).

USAGE

    snapfit host, guest [, kind [, diameter [, length [, clearance
            [, magnet_d [, magnet_h [, wall ]]]]]]]

    kind = peg | magnet
    clearance = gap for the fit in mm (0 = auto: half the nozzle size)
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


def all_presets():
    """Built-in presets plus the user's own (user names win)."""
    data = {k: {sk: dict(sv) for sk, sv in v.items()} for k, v in BUILTIN_PRESETS.items()}
    data.update(_load_presets())
    return data


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
    for k in ('struts', 'split', 'export', 'snapfit', 'surface'):
        if k in _last:
            data[k] = dict(_last[k])
    return data


def print3d_preset(action='list', name='', _self=cmd):
    """
DESCRIPTION

    Save, load, list or delete named print presets (thickness, surface,
    struts, split, export and snap-fit options), stored in
    print3d_presets.json. Built-in presets: "Multi-color Printer Cartoon",
    "Single Color Printer Cartoon", "Multi-color Printer Surface",
    "Single Color Printer Surface" (can be loaded, not overwritten).

USAGE

    print3d_preset list
    print3d_preset save, DN3      # current thickness + last-used options
    print3d_preset load, DN3      # applies thickness / surface settings;
                                  # export options are used by:
                                  # colorprint file, preset=DN3
    print3d_preset delete, DN3
    """
    s = _self
    action, name = str(action).lower(), str(name).strip()
    data = _load_presets()
    if action == 'list':
        _msg(['built-in: ' + ', '.join(BUILTIN_PRESETS),
              'yours: ' + (', '.join(sorted(data)) or '(none)')])
        return list(BUILTIN_PRESETS) + sorted(data)
    if not name:
        print(' print3d_preset: give a name')
        return
    if action in ('save', 'delete') and name in BUILTIN_PRESETS:
        print(' print3d_preset: "%s" is built in; use another name' % name)
        return
    if action == 'save':
        data[name] = current_settings(s)
        _save_presets(data)
        _msg(['saved preset "%s" to %s' % (name, PRESET_FILE)])
    elif action == 'load':
        p = all_presets().get(name)
        if p is None:
            print(' print3d_preset: no preset "%s"' % name)
            return
        if 'style' in p:
            printstyle(_self=s, **p['style'])
        if 'surface' in p:
            surfacestyle(_self=s, **p['surface'])
        _msg(['loaded preset "%s"' % name])
        return p
    elif action == 'delete':
        if data.pop(name, None) is not None:
            _save_presets(data)
            _msg(['deleted preset "%s"' % name])


def auto_filaments(objects='', slots=4, _self=cmd):
    """
DESCRIPTION

    Suggest filaments: the common PLA colors closest to the colors shown,
    most used first, at most `slots`. Use the result as colorprint filaments=
    (or just pass filaments=auto to colorprint).

USAGE

    auto_filaments [objects [, slots ]]
    """
    s = _self
    names = [n for n in str(objects).replace(',', ' ').split() if n]
    vr = _vrml_of(s, names)
    pieces = pc.vrml_to_pieces(vr, 8, True, True, [], {})
    picked = pc.pick_filaments(pieces, int(float(slots)) or 4)
    _msg(['filaments: ' + ', '.join('%d=%s (%s, %.0f%%)' % (i + 1, n, h, 100 * w)
                                    for i, (n, h, w) in enumerate(picked)),
          'use: colorprint file, filaments=' + ' '.join('0x' + h.lstrip('#').lower()
                                                        for _, h, _ in picked)])
    return picked


def print3d_printers(_self=cmd):
    """
DESCRIPTION

    List the printer profiles (build volume, filament slots) usable with
    colorprint ..., printer=NAME.
    """
    for n, d in PRINTERS.items():
        print(' %-34s %d x %d x %d mm, %d slot%s' % ((n,) + tuple(d['bed']) +
                                                    (d['slots'], '' if d['slots'] == 1 else 's')))
    return list(PRINTERS)


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
cmd.extend('strut_add', strut_add)
cmd.extend('strut_remove', strut_remove)
cmd.extend('strut_remove_selected', strut_remove_selected)
cmd.extend('struts_label', struts_label)
cmd.extend('ligand_surface', ligand_surface)
cmd.extend('split_parts', split_parts)
cmd.extend('surfacestyle', surfacestyle)
cmd.extend('print3d_printers', print3d_printers)
cmd.extend('auto_filaments', auto_filaments)
cmd.extend('show_strand_joins', show_strand_joins)
cmd.extend('colorprint', colorprint)
cmd.extend('snapfit', snapfit)
cmd.extend('snapfit_clear', snapfit_clear)
cmd.extend('print3d_preset', print3d_preset)
cmd.extend('print3d_install_merge', print3d_install_merge)
cmd.auto_arg[0]['struts'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['strut_add'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['strut_add'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['strut_remove'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[1]['strut_remove'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['ligand_surface'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['split_parts'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['surfacestyle'] = cmd.auto_arg[0]['zoom']
cmd.auto_arg[0]['snapfit'] = cmd.auto_arg[0]['disable']
cmd.auto_arg[1]['snapfit'] = cmd.auto_arg[0]['disable']
