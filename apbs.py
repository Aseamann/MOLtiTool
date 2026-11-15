from pymol import cmd
import os
import glob
import shutil
import tempfile
import subprocess


# ---------------------------------------------------------------------------
# APBS input template
# ---------------------------------------------------------------------------

APBS_TEMPLATE = """\
read
    mol pqr "{pqrfile}"
end

elec
    mg-auto
    mol 1

    fgcent {fgx:.3f} {fgy:.3f} {fgz:.3f}
    cgcent mol 1

    fglen {flx:.3f} {fly:.3f} {flz:.3f}
    cglen {clx:.3f} {cly:.3f} {clz:.3f}

    dime {dx} {dy} {dz}

    lpbe

    bcfl sdh

    pdie 2.0
    sdie 78.0

    chgm spl2
    srfm smol

    swin 0.3
    temp 310.0

    sdens 10.0

    calcenergy no
    calcforce no

    srad {srad:.3f}

    ion charge +1 conc 0.15 radius 2.0
    ion charge -1 conc 0.15 radius 1.8

    write pot dx "{mapfile}"
end

quit
"""


# ---------------------------------------------------------------------------
# Valid APBS grid sizes
# ---------------------------------------------------------------------------

_VALID_DIME = (
    33, 65, 97, 129, 161, 193, 225, 257,
    289, 321, 353, 385, 417, 449, 513,
    577, 641, 705, 769, 833, 897, 961, 1025
)


def _round_up_dime(n):
    for v in _VALID_DIME:
        if v >= n:
            return v
    return _VALID_DIME[-1]


def _compute_grid(extent, spacing=1.5, padding=10.0):

    (xmn, ymn, zmn), (xmx, ymx, zmx) = extent

    fgcent = (
        (xmn + xmx) / 2.0,
        (ymn + ymx) / 2.0,
        (zmn + zmx) / 2.0,
    )

    box = (
        xmx - xmn,
        ymx - ymn,
        zmx - zmn,
    )

    fglen = tuple(b + 2.0 * padding for b in box)

    cglen = tuple(
        max(f * 1.7, f + 30.0)
        for f in fglen
    )

    dime = tuple(
        _round_up_dime(int(round(f / spacing)) + 1)
        for f in fglen
    )

    return fgcent, fglen, cglen, dime


# ---------------------------------------------------------------------------
# Executable discovery
# ---------------------------------------------------------------------------

def _find_executables():

    import sys

    apbs = shutil.which("apbs")

    pdb2pqr = (
        shutil.which("pdb2pqr")
        or shutil.which("pdb2pqr30")
        or shutil.which("pdb2pqr.py")
    )

    if not (apbs and pdb2pqr):

        for d in [
            os.path.join(sys.prefix, "bin"),
            os.path.dirname(sys.executable),
        ]:

            if not apbs:
                hit = glob.glob(os.path.join(d, "apbs"))
                if hit:
                    apbs = hit[0]

            if not pdb2pqr:
                hit = (
                    glob.glob(os.path.join(d, "pdb2pqr"))
                    + glob.glob(os.path.join(d, "pdb2pqr30"))
                )
                if hit:
                    pdb2pqr = hit[0]

    return apbs, pdb2pqr


# ---------------------------------------------------------------------------
# Subprocess helper
# ---------------------------------------------------------------------------

def _run(cmd_list, cwd, label):

    result = subprocess.run(
        cmd_list,
        cwd=cwd,
        capture_output=True,
        text=True,
    )

    if result.returncode != 0:

        print(f"[{label}] FAILED")

        if result.stderr:
            print(result.stderr)

        if result.stdout:
            print(result.stdout)

        return False

    return True


# ---------------------------------------------------------------------------
# Main APBS function
# ---------------------------------------------------------------------------

def apbs_all(
    ph=7.0,
    ff="PARSE",
    grid_spacing=1.5,
    padding=10.0,
    probe=1.4,
    lo=-5.0,
    hi=5.0,
    selection_filter="polymer and not solvent",
):

    apbs_exe, pdb2pqr_exe = _find_executables()

    if not (apbs_exe and pdb2pqr_exe):
        print("Could not find APBS or PDB2PQR.")
        return

    objects = cmd.get_object_list("polymer")

    if not objects:
        print("No polymer objects found.")
        return

    reference = objects[0]

    print(f"Reference: {reference}")

    for obj in objects:

        print(f"\n=== {obj} ===")

        if obj != reference:

            try:
                rms, natoms, *_ = cmd.align(obj, reference)
                print(f"Aligned RMSD = {rms:.3f} Å")

            except Exception as e:
                print(e)
                continue

        prep_name = f"{obj}_prepared"
        map_name = f"{obj}_map"
        surf_name = f"{obj}_surface"
        ramp_name = f"{obj}_ramp"
        group_name = f"{obj}_apbs"

        for n in [
            prep_name,
            map_name,
            surf_name,
            ramp_name,
            group_name,
        ]:
            if n in cmd.get_names("all"):
                cmd.delete(n)

        sel = f"({obj}) and ({selection_filter})"

        with tempfile.TemporaryDirectory(
            prefix=f"apbs_{obj}_"
        ) as wd:

            pdb_in = "input.pdb"
            pqr_out = "input.pqr"
            apbs_in = "apbs.in"
            map_stem = "potential"

            cmd.save(os.path.join(wd, pdb_in), sel)

            pdb2pqr_cmd = [
                pdb2pqr_exe,
                f"--ff={ff}",
                "--with-ph",
                str(ph),
                pdb_in,
                pqr_out,
            ]

            if not _run(pdb2pqr_cmd, wd, "pdb2pqr"):
                continue

            extent = cmd.get_extent(sel)

            fgcent, fglen, cglen, dime = _compute_grid(
                extent,
                spacing=grid_spacing,
                padding=padding,
            )

            with open(os.path.join(wd, apbs_in), "w") as fh:

                fh.write(
                    APBS_TEMPLATE.format(
                        pqrfile=pqr_out,
                        mapfile=map_stem,

                        fgx=fgcent[0],
                        fgy=fgcent[1],
                        fgz=fgcent[2],

                        flx=fglen[0],
                        fly=fglen[1],
                        flz=fglen[2],

                        clx=cglen[0],
                        cly=cglen[1],
                        clz=cglen[2],

                        dx=dime[0],
                        dy=dime[1],
                        dz=dime[2],

                        srad=probe,
                    )
                )

            if not _run([apbs_exe, apbs_in], wd, "apbs"):
                continue

            dx_files = glob.glob(os.path.join(wd, "*.dx"))

            if not dx_files:
                print("No DX map generated.")
                continue

            cmd.load(dx_files[0], map_name)

        # -------------------------------------------------------------------
        # Surface rendering
        # -------------------------------------------------------------------

        cmd.create(surf_name, obj)

        cmd.hide("everything", surf_name)
        cmd.show("surface", surf_name)

        # geometry stays correct
        cmd.set("surface_mode", 1, surf_name)

        # mild smoothing only (does NOT inflate surface)
        cmd.set("surface_quality", 2, surf_name)
        cmd.set("surface_smooth_edges", 1, surf_name)

        # IMPORTANT: keep coloring smooth without changing geometry
        cmd.set("surface_ramp_above_mode", 1)

        # -------------------------------------------------------------------
        # Electrostatic coloring
        # -------------------------------------------------------------------

        cmd.ramp_new(
            ramp_name,
            map_name,
            [float(lo), 0.0, float(hi)],
            ["red", "white", "blue"]
        )

        cmd.set("surface_color", ramp_name, surf_name)

        # -------------------------------------------------------------------
        # Group outputs
        # -------------------------------------------------------------------

        cmd.group(
            group_name,
            f"{map_name} {surf_name} {ramp_name}"
        )

        print(f"Generated: {group_name}")

    print("\nDone.")


cmd.extend("apbs_all", apbs_all)