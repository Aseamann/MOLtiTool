from pymol import cmd


def residue_shell(selection, dist=5.0, name=None, color="red", width=0.5):
    """
    Draw a thin mesh outline at `dist` A from every heavy atom of `selection`.

    Usage (PyMOL command line):
        run /path/to/MOLtiTool/residue_shell.py
        residue_shell chain A and resi 120
        residue_shell chain A and resi 120, 5.0, mycolor_shell, red, 0.5

    Arguments:
        selection  any PyMOL selection (a residue, a ligand, several residues)
        dist       shell distance in angstroms (default 5.0)
        name       name of the new object (default: shell_<dist>)
        color      mesh color (default red)
        width      mesh line width (default 0.5)

    How it works: a hydrogen-free copy of the selection is made, every atom's
    radius is set to `dist`, and the surface of that copy is drawn as a mesh with
    the solvent radius at 0. The mesh is the boundary of the union of 5 A spheres,
    i.e. the set of points exactly `dist` from the nearest atom. The original
    objects are not changed. Delete the outline with: delete <name>
    """
    dist = float(dist)
    if name is None:
        name = f"shell_{dist:g}"

    cmd.delete(name)
    cmd.create(name, f"({selection}) and not hydro")
    cmd.hide("everything", name)

    # every atom gets a radius equal to the shell distance
    cmd.alter(name, f"vdw={dist}")
    cmd.rebuild(name)

    # surface at exactly the atom radii (no extra solvent probe)
    cmd.set("solvent_radius", 0.0, name)
    cmd.set("surface_quality", 1, name)
    cmd.set("mesh_color", color, name)
    cmd.set("mesh_width", float(width), name)
    cmd.show("mesh", name)
    cmd.enable(name)


cmd.extend("residue_shell", residue_shell)
