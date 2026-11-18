"""
pca_noise_gradient.py  -  noise a copy of a structure, with noise growing along a chosen axis.

Usage:
    run /path/to/MOLtiTool/pca_noise_gradient.py
    pca_noise_gradient obj [, out [, sel [, sigma_min [, sigma_max ]]]] [, axis=pc1 ] [, flip=1 ] ...

The noise axis can be a principal axis of the selection (pc1/pc2/pc3), the screen
(screen_x / screen_y: turn the molecule until the gradient runs the way you want), or
the line from one selection to another. flip=1 reverses it; tilt rotates it.
atoms=backbone / ca / sidechain noises only those atoms; they are shown as spheres.
"""
from typing import Literal

import numpy as np
from pymol import cmd

AXES = ("pc1", "pc2", "pc3", "screen_x", "screen_y", "selections")
SS_MODES = ("ramp", "dss", "loop", "keep")
ATOM_SETS = {"all": "all", "backbone": "backbone", "ca": "name CA", "sidechain": "sidechain"}
SHOW_MODES = ("spheres", "split", "spheres+cartoon", "cartoon", "keep")
COMPACT_MODES = ("none", "toward_noisy_end", "toward_global")


def _truthy(v):
    return str(v).strip().lower() in ("1", "true", "yes", "on")


def _unit(v):
    n = np.linalg.norm(v)
    if n < 1e-8:
        raise ValueError("noise axis has zero length")
    return v / n


def _noise_axis(axis, X, obj, state, tilt, axis_from, axis_to, _self):
    """Unit vector (model space) along which noise grows, before flip."""
    if axis in ("pc1", "pc2", "pc3"):
        _, _, Vt = np.linalg.svd(X, full_matrices=False)
        pcs = [Vt[k] if k < len(Vt) else np.zeros(3) for k in range(3)]
        order = np.arange(len(X), dtype=float)
        for k in range(3):
            # SVD signs are arbitrary: point each PC toward the later atoms (N -> C)
            # so the noisy end is reproducible
            if np.dot(X @ pcs[k], order - order.mean()) < 0:
                pcs[k] = -pcs[k]
        k = int(axis[2]) - 1
        a, partner = pcs[k], pcs[(k + 1) % 3]
    elif axis in ("screen_x", "screen_y"):
        R = np.array(_self.get_view()[:9], dtype=float).reshape(3, 3)  # camera = xyz @ R
        sx, sy = R[:, 0], R[:, 1]
        a, partner = (sx, sy) if axis == "screen_x" else (sy, -sx)
    else:  # selections
        if not axis_from or not axis_to:
            raise ValueError("axis=selections needs axis_from and axis_to")
        c = []
        for s in (axis_from, axis_to):
            xyz = _self.get_coords(f"({obj}) and ({s})", state=state)
            if xyz is None or not len(xyz):
                raise ValueError(f"No atoms in {obj} for axis selection: {s}")
            c.append(xyz.mean(axis=0))
        a, partner = c[1] - c[0], None
    a = _unit(np.asarray(a, dtype=float))
    if tilt:
        if partner is None:
            print("[pca_noise_gradient] tilt is ignored with axis=selections")
        else:
            p = partner - np.dot(partner, a) * a
            if np.linalg.norm(p) > 1e-8:
                th = np.radians(tilt)
                a = _unit(np.cos(th) * a + np.sin(th) * _unit(p))
    return a


def pca_noise_gradient(
    obj,
    out="",
    sel="",
    sigma_min=0.0,
    sigma_max=2.0,
    seed=0,
    state=1,
    byres=0,
    store_b=1,
    ramp_start=0.0,
    ramp_power=1.0,
    atoms: Literal["all", "backbone", "ca", "sidechain"] = "all",
    axis: Literal["pc1", "pc2", "pc3", "screen_x", "screen_y", "selections"] = "pc1",
    flip=0,
    tilt=0.0,
    axis_from="",
    axis_to="",
    ss: Literal["ramp", "dss", "loop", "keep"] = "ramp",
    ss_cutoff=0.5,
    show: Literal["spheres", "split", "spheres+cartoon", "cartoon", "keep"] = "spheres",
    stretch=0.0,
    compact_mode: Literal["none", "toward_noisy_end", "toward_global"] = "none",
    compact_amp=1.0,
    compact_size=4.0,
    compact_start=0.6,
    compact_power=0.5,
    compact_offset=0.0,
    noisy_center_window=0.15,
    _self=cmd,
    **_old,
):
    """
DESCRIPTION

    Copy obj to out and add Gaussian noise to the atoms in sel. The noise SD grows
    from sigma_min (A) at the quiet end of the axis to sigma_max (A) at the noisy end.
    With store_b=1 the ramp (0-100) goes in the B-factor: blue = quiet, red = noisy.

    atoms  which atoms of sel get noise: all, backbone (N CA C O), ca, or sidechain.
           The axis and ramp still come from all of sel, so switching this does not
           move the gradient. With byres=1 the chosen atoms move together per residue.

DIRECTION

    axis   pc1 / pc2 / pc3: principal axes of sel (oriented so noise grows toward
           the later atoms, usually the C-terminal side).
           screen_x: grows left -> right on screen. screen_y: bottom -> top.
           Rotate the molecule first to point the gradient anywhere.
           selections: from the center of axis_from toward the center of axis_to.
    flip   1 = reverse the direction.
    tilt   degrees to rotate the axis toward the next one
           (pc1->pc2, pc2->pc3, pc3->pc1, screen_x->screen_y, screen_y->left).

OTHER

    ramp_start   fraction along the axis where noise starts (0-1).
    ramp_power   >1 keeps the noise low longer, then rises sharply.
    byres        1 = move whole residues rigidly (keeps side chains intact).
    show         spheres: show only the atoms that got noise, as spheres (default).
                 split: cartoon where the noise SD is below ss_cutoff, spheres of
                 the noised atoms above it (folded -> noise look).
                 spheres+cartoon: those spheres plus the cartoon of the copy.
                 cartoon: cartoon only. keep: representations copied from obj.
    ss           ramp: keep the original helix/sheet where the noise SD is below
                 ss_cutoff (A), draw the rest as loop (default; clean gradient).
                 dss: recompute secondary structure on the noisy copy.
                 loop: draw everything as loop (shows the most noise).
                 keep: keep the original helices/sheets. PyMOL smooths those
                 ribbons, which hides most of the noise.

STRETCH

    stretch      shift (A) along the noise direction, 0 at the quiet end growing to
                 this at the noisy end (follows the ramp). Pulls the noisy part past
                 the end of the original structure. Negative pushes it back.

COMPACTION (the noisy end condenses into a cluster; off with compact_mode=none)

    compact_mode  none: off.
                  toward_noisy_end: pull toward the center of the noisiest atoms
                  (the last noisy_center_window of the ramp), plus compact_offset.
                  toward_global: pull toward the center of the whole selection.
    compact_size  size (A, SD per axis) of the cluster the noisy end condenses
                  into. Set directly, so it does not grow with sigma_max or with
                  the size of the structure. ~3-6 A for a tight ball of CAs.
    compact_amp   how far atoms at the noisy end go into the cluster (0-1);
                  1 = fully, lower leaves some of the structure and noise.
    compact_start where along the ramp (0-1) the pull begins. Below it the noise
                  spreads freely; above it atoms condense into the cluster.
    compact_power how fast the pull rises after compact_start: 0.5 (default)
                  rises quickly so most of that zone condenses, 1 = linear,
                  >1 = only the very end condenses.
    compact_offset  move the cluster center this far (A) past the noisy end,
                  along the noise direction.
    """
    sigma_min, sigma_max = float(sigma_min), float(sigma_max)
    seed, state = int(seed), int(state)
    byres, store_b, flip = _truthy(byres), _truthy(store_b), _truthy(flip)
    ramp_start = min(1.0, max(0.0, float(ramp_start)))
    ramp_power = max(1e-6, float(ramp_power))
    tilt = float(tilt)
    ss_cutoff = float(ss_cutoff)
    axis = str(axis).strip().lower()
    ss = str(ss).strip().lower()
    show = str(show).strip().lower()
    stretch = float(stretch)
    compact_mode = str(compact_mode).strip().lower()
    if compact_mode == "toward_noisy_center":       # old name
        compact_mode = "toward_noisy_end"
    compact_amp = min(1.0, max(0.0, float(compact_amp)))
    compact_size = max(0.0, float(compact_size))
    compact_power = max(1e-6, float(compact_power))
    compact_start = min(0.99, max(0.0, float(compact_start)))
    compact_offset = float(compact_offset)
    noisy_center_window = min(1.0, max(1e-6, float(noisy_center_window)))
    if _old:
        print("[pca_noise_gradient] ignoring old options: " + ", ".join(sorted(_old)))

    atoms = str(atoms).strip().lower()
    for name, val, allowed in (("show", show, SHOW_MODES), ("atoms", atoms, tuple(ATOM_SETS)), ("axis", axis, AXES), ("ss", ss, SS_MODES),
                               ("compact_mode", compact_mode, COMPACT_MODES)):
        if val not in allowed:
            raise ValueError(f"{name} must be one of: {', '.join(allowed)}")

    out = out or f"{obj}_pcaNoise"
    if out in _self.get_names("objects"):
        _self.delete(out)
    _self.create(out, obj, state, 1)

    # Work on the copy itself, so atoms are matched by index (no per-atom selections)
    target = out if not sel else f"({out}) and ({sel})"
    model = _self.get_model(target, state=1)
    if not model.atom:
        _self.delete(out)
        raise ValueError(f"No atoms found for selection: {sel or obj}")
    idx = [a.index for a in model.atom]
    coords = np.array([a.coord for a in model.atom], dtype=float)
    X = coords - coords.mean(axis=0)

    u = _noise_axis(axis, X, obj, state, tilt, axis_from, axis_to, _self)
    if flip:
        u = -u

    proj = X @ u
    pmin, pmax = float(proj.min()), float(proj.max())
    t = (proj - pmin) / ((pmax - pmin) if pmax - pmin > 1e-8 else 1.0)
    t = np.clip((t - ramp_start) / (1.0 - ramp_start + 1e-12), 0.0, 1.0) ** ramp_power
    sigmas = sigma_min + t * (sigma_max - sigma_min)

    # Only these atoms move (the axis above used all of sel)
    if atoms != "all":
        chosen = set()
        _self.iterate(f"({target}) and ({ATOM_SETS[atoms]})", "_s.add(index)",
                      space={"_s": chosen})
        if not chosen:
            _self.delete(out)
            raise ValueError(f"No {atoms} atoms in the selection")
        moving = np.array([i in chosen for i in idx])
    else:
        moving = np.ones(len(idx), dtype=bool)

    # Residue groups (byres moves each group rigidly)
    groups = {}
    for i, a in enumerate(model.atom):
        if moving[i]:
            key = (a.segi, a.chain, a.resi, a.resn) if byres else i
            groups.setdefault(key, []).append(i)
    groups = list(groups.values())

    def rigid(per_atom):
        """byres: give every atom of a residue its CA's (or the mean) value."""
        if not byres:
            return per_atom
        outv = np.array(per_atom, dtype=float)
        for ii in groups:
            ca = next((i for i in ii if model.atom[i].name == "CA"), None)
            outv[ii] = per_atom[ca] if ca is not None else per_atom[ii].mean(axis=0)
        return outv

    # 1. noise
    rng = np.random.default_rng(seed)
    noise = rng.normal(0.0, 1.0, size=coords.shape) * sigmas[:, None]
    pos = coords + rigid(noise)

    # 2. stretch along the noise direction
    if stretch:
        pos = pos + rigid((stretch * t)[:, None] * u)

    # 3. compaction, after noise and stretch: blend each atom toward a spot in a
    #    Gaussian cluster (SD compact_size) around the center, by
    #    f = compact_amp * ((t - start) / (1 - start))^power
    if compact_mode != "none" and compact_amp > 0.0:
        if compact_mode == "toward_noisy_end":
            end = t >= 1.0 - noisy_center_window
            center = coords[end].mean(axis=0) if end.any() else coords.mean(axis=0)
            center = center + stretch * u          # the end moved with the stretch
        else:
            center = coords.mean(axis=0)
        center = center + compact_offset * u
        frac = compact_amp * np.clip((t - compact_start) / (1.0 - compact_start),
                                     0.0, 1.0) ** compact_power
        spot = center + compact_size * rng.normal(0.0, 1.0, size=coords.shape)
        pos = pos + rigid(frac[:, None] * (spot - pos))

    pos[~moving] = coords[~moving]
    new = {i: tuple(p) for i, p in zip(idx, pos)}
    _self.alter_state(1, target, "(x, y, z) = _new[index]", space={"_new": new})

    if store_b:
        bvals = {i: float(v * 100.0) for i, v in zip(idx, t)}
        _self.alter(target, "b = _b[index]", space={"_b": bvals})
        _self.spectrum("b", "blue_white_red", target, minimum=0, maximum=100)

    if ss == "ramp":
        sig = {i: float(v) for i, v in zip(idx, sigmas)}
        _self.alter(target, "ss = 'L' if _sig[index] >= _cut else ss",
                    space={"_sig": sig, "_cut": ss_cutoff})
    elif ss == "dss":
        _self.dss(out, state=1)
    elif ss == "loop":
        _self.alter(out, "ss = 'L'")
    _self.rebuild(out)
    noised = f"({target}) and ({ATOM_SETS[atoms]})"
    if show != "keep":
        _self.hide("everything", out)
    if show in ("spheres", "spheres+cartoon"):
        _self.show("spheres", noised)
    if show == "split":
        hi = [i for i, v in zip(idx, sigmas) if v >= ss_cutoff]
        _self.select_list("_png_hi", out, hi or [0], mode="index")
        _self.show("spheres", f"_png_hi and ({noised})")
        _self.show("cartoon", f"{out} and not (byres (_png_hi and ({noised})))")
        _self.delete("_png_hi")
    if show in ("cartoon", "spheres+cartoon"):
        _self.show("cartoon", out)

    hot = model.atom[int(np.argmax(proj))]
    print(
        f"[pca_noise_gradient] Created '{out}' from '{obj}': sigma {sigma_min}-{sigma_max} A, "
        f"atoms={atoms}, axis={axis}{' (flipped)' if flip else ''}{f' tilt={tilt:g}' if tilt else ''}, "
        f"noisiest end near {hot.chain}/{hot.resn}{hot.resi}, ss={ss}, "
        f"stretch={stretch:g}, compact_mode={compact_mode}"
        + (f" size={compact_size:g} amp={compact_amp:g} start={compact_start:g}"
           f" power={compact_power:g}"
           f" offset={compact_offset:g}"
           if compact_mode != "none" else "")
    )


cmd.extend("pca_noise_gradient", pca_noise_gradient)
cmd.auto_arg[0]["pca_noise_gradient"] = cmd.auto_arg[0]["zoom"]
