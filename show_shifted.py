from pymol import cmd


def show_shifted(obj1, obj2, cutoff=1.5, do_align=1):
    """
    Show sticks for residues whose CA moves more than `cutoff` A between obj1 and obj2.
    Cartoon and existing colors are kept; hydrogens are hidden.

    Usage (PyMOL command line):
        run /path/to/MOLtiTool/show_shifted.py
        show_shifted obj1, obj2, 1.5
        show_shifted obj1, obj2, 1.5, 0     # skip alignment (already superposed)

    Residues are matched by chain and residue number.
    """
    cutoff = float(cutoff)
    if int(do_align):
        # cycles=0 -> plain least-squares fit, no outlier rejection,
        # so the big movers don't get down-weighted out of the superposition
        cmd.align(f"{obj1} and name CA", f"{obj2} and name CA", cycles=0)

    ca1, ca2 = {}, {}
    cmd.iterate_state(1, f"{obj1} and name CA", "d[(chain, resi)] = (x, y, z)", space={"d": ca1})
    cmd.iterate_state(1, f"{obj2} and name CA", "d[(chain, resi)] = (x, y, z)", space={"d": ca2})

    moved = {}
    for key in ca1.keys() & ca2.keys():
        a, b = ca1[key], ca2[key]
        dist = sum((p - q) ** 2 for p, q in zip(a, b)) ** 0.5
        if dist > cutoff:
            moved.setdefault(key[0], []).append(key[1])

    if not moved:
        print(f"No residues shift more than {cutoff} A")
        return

    clauses = [f"(chain {ch} and resi {'+'.join(resis)})" for ch, resis in moved.items()]
    sel = " or ".join(clauses)
    cmd.select("shifted", f"({obj1} or {obj2}) and ({sel})")

    # keep existing cartoon and colors; just add sticks for shifted residues
    cmd.hide("everything", "hydro")
    cmd.show("cartoon", f"{obj1} or {obj2}")
    cmd.show("sticks", "shifted and not hydro")
    cmd.zoom("shifted")
    print(f"{sum(len(v) for v in moved.values())} residues shift more than {cutoff} A")


cmd.extend("show_shifted", show_shifted)
