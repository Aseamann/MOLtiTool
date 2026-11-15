from pymol import cmd
import numpy as np


def _safe_sel_atom(out_obj, a):
    s = f"{out_obj} and chain {a.chain} and resi {a.resi} and name {a.name}"
    if a.alt:
        s += f" and alt {a.alt}"
    if a.segi:
        s += f" and segi {a.segi}"
    return f"({s})"


def pca_noise_gradient(
    obj,
    out=None,
    sel=None,
    sigma_min=0.0,
    sigma_max=2.0,
    seed=0,
    state=1,
    byres=0,
    store_b=1,
    ramp_start=0.0,
    ramp_power=1.0,
    # --- compaction controls (THIS is what reduces Rg) ---
    compact_mode="none",          # "none", "toward_global", "toward_noisy_center"
    compact_amp=0.0,              # Å: max drift magnitude at peak
    compact_center=0.5,           # center in ramped t-space (0..1), usually 0.5
    compact_width=0.25,           # width in t-units for where compaction is active
    compact_power=1.0,            # >1 concentrates compaction near peak
    noisy_center_window=0.15,     # for toward_noisy_center: window around compact_center
):
    """
    Add Gaussian positional noise with magnitude increasing along PC1 (primary PCA axis),
    optionally delaying onset (ramp_start) and adding *radial compaction* to lower Rg
    in/near the noisy region.

    Key parameters:
      - sigma_max: increases noise at the noisy end
      - ramp_start: noise begins after this fraction along PC1
      - compact_mode: pulls atoms inward (reduces Rg) near noisy region
      - compact_amp: strength of inward pull (Å) at peak
    """

    # Coerce types
    sigma_min = float(sigma_min)
    sigma_max = float(sigma_max)
    seed = int(seed)
    state = int(state)
    byres = int(byres)
    store_b = int(store_b)
    ramp_start = float(ramp_start)
    ramp_power = float(ramp_power)

    compact_mode = str(compact_mode).strip()
    compact_amp = float(compact_amp)
    compact_center = float(compact_center)
    compact_width = float(compact_width)
    compact_power = float(compact_power)
    noisy_center_window = float(noisy_center_window)

    # Clamp/validate
    ramp_start = max(0.0, min(1.0, ramp_start))
    ramp_power = max(1e-6, ramp_power)

    compact_center = max(0.0, min(1.0, compact_center))
    compact_width = max(1e-6, compact_width)
    compact_power = max(1e-6, compact_power)
    noisy_center_window = max(1e-6, noisy_center_window)

    if compact_mode not in ("none", "toward_global", "toward_noisy_center"):
        raise ValueError("compact_mode must be one of: none, toward_global, toward_noisy_center")

    if out is None or out == "":
        out = f"{obj}_pcaNoise"

    full_sel = obj if (sel is None or sel == "") else f"({obj}) and ({sel})"

    model = cmd.get_model(full_sel, state=state)
    if not model.atom:
        raise ValueError(f"No atoms found for selection: {full_sel}")

    coords = np.array([a.coord for a in model.atom], dtype=float)

    # Center coords and find PC1
    mean = coords.mean(axis=0)
    X = coords - mean
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    pc1 = Vt[0] / np.linalg.norm(Vt[0])

    # Project and build ramped t in 0..1
    proj = X.dot(pc1)
    pmin, pmax = float(proj.min()), float(proj.max())
    denom = (pmax - pmin) if (pmax - pmin) > 1e-8 else 1.0
    t_raw = (proj - pmin) / denom

    # delayed ramp
    t = np.clip((t_raw - ramp_start) / (1.0 - ramp_start + 1e-12), 0.0, 1.0)
    t = t ** ramp_power

    # sigma along ramp
    sigmas = sigma_min + t * (sigma_max - sigma_min)

    # Compaction weight peaked at compact_center in *ramped* t-space
    w_comp = np.exp(-0.5 * ((t - compact_center) / compact_width) ** 2)
    w_comp = w_comp ** compact_power  # sharpen if desired

    # Choose compaction center
    if compact_mode == "toward_global":
        comp_center_xyz = coords.mean(axis=0)
    elif compact_mode == "toward_noisy_center":
        # Use atoms near the noisy peak region to define a "local center"
        mask = np.abs(t - compact_center) <= noisy_center_window
        if np.any(mask):
            comp_center_xyz = coords[mask].mean(axis=0)
        else:
            # fallback if mask is empty
            comp_center_xyz = coords.mean(axis=0)
    else:
        comp_center_xyz = None  # unused

    # Create output copy
    if out in cmd.get_names("objects"):
        cmd.delete(out)
    cmd.create(out, obj, state, 1)

    rng = np.random.default_rng(seed)

    def _compaction_drift(pos_xyz, wval):
        """Drift toward chosen center to reduce Rg (radial inward pull)."""
        if compact_mode == "none" or compact_amp == 0.0:
            return np.zeros(3, dtype=float)
        v = comp_center_xyz - pos_xyz
        n = np.linalg.norm(v)
        if n < 1e-8:
            return np.zeros(3, dtype=float)
        # unit direction toward center * amplitude * weight
        return (compact_amp * float(wval)) * (v / n)

    if byres:
        atoms = model.atom
        keys = [(a.chain, a.resi, a.segi) for a in atoms]
        res_to_idxs = {}
        for i, k in enumerate(keys):
            res_to_idxs.setdefault(k, []).append(i)

        for idxs in res_to_idxs.values():
            # representative position for drift + representative sigma
            # prefer CA if present
            ca_i = None
            for i in idxs:
                if atoms[i].name == "CA":
                    ca_i = i
                    break

            if ca_i is None:
                pos_r = coords[idxs].mean(axis=0)
                sigma_r = float(sigmas[idxs].mean())
                w_r = float(w_comp[idxs].mean())
            else:
                pos_r = coords[ca_i]
                sigma_r = float(sigmas[ca_i])
                w_r = float(w_comp[ca_i])

            disp = rng.normal(0.0, sigma_r, size=3) + _compaction_drift(pos_r, w_r)

            for i in idxs:
                a = atoms[i]
                atom_sel = _safe_sel_atom(out, a)
                cmd.translate(list(disp), atom_sel, state=1, camera=0)

    else:
        for i, a in enumerate(model.atom):
            pos_i = coords[i]
            disp = rng.normal(0.0, float(sigmas[i]), size=3) + _compaction_drift(pos_i, float(w_comp[i]))
            atom_sel = _safe_sel_atom(out, a)
            cmd.translate(list(disp), atom_sel, state=1, camera=0)

    # Store ramp in B and color
    if store_b:
        for i, a in enumerate(model.atom):
            atom_sel = _safe_sel_atom(out, a)
            cmd.alter(atom_sel, f"b={float(t[i]*100.0)}")
        cmd.sort(out)
        cmd.spectrum("b", "blue_white_red", out, minimum=0, maximum=100)

    cmd.show("cartoon", out)
    print(
        f"[pca_noise_gradient] Created '{out}' from '{obj}' "
        f"sigma=[{sigma_min},{sigma_max}] ramp_start={ramp_start} ramp_power={ramp_power} "
        f"compact_mode={compact_mode} compact_amp={compact_amp} compact_center={compact_center} "
        f"compact_width={compact_width} compact_power={compact_power}"
    )


cmd.extend("pca_noise_gradient", pca_noise_gradient)

