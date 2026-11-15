"""
align_chains.py

Superpose one chain from every loaded object onto the same chain of a
reference object. The whole object moves, but only the selected chain
drives the fit.

Load it inside PyMOL with:

    run /path/to/align_chains.py

then call `align_chains` from the command line or from Python.
"""

from fnmatch import fnmatch

from pymol import cmd

try:
    from pymol import CmdException
except ImportError:  # very old PyMOL
    CmdException = RuntimeError

__all__ = ["align_chains"]


def _as_bool(value):
    if isinstance(value, str):
        return value.strip().lower() not in ("0", "false", "no", "off", "none", "")
    return bool(value)


def _as_int(value):
    return int(float(value))


def _blank(value):
    return value is None or (isinstance(value, str) and value.strip() in ("", "None"))


def align_chains(chain="A", ref=None, mobile="*", ref_chain=None, method="align",
                 ca_only=1, cycles=5, extra="", transform=1, aln_object="",
                 sort=1, quiet=0):
    """
DESCRIPTION

    Superpose one chain of every loaded object onto the same chain of a
    reference object. Only the chosen chain drives the fit, but the entire
    object is moved, so other chains and ligands follow along.

USAGE

    align_chains [ chain [, ref [, mobile [, ref_chain [, method [, ca_only
                 [, cycles [, extra [, transform [, aln_object [, sort
                 [, quiet ]]]]]]]]]]]]

ARGUMENTS

    chain       chain ID taken from each mobile object       {default: A}
    ref         reference object, first loaded if omitted    {default: None}
    mobile      space separated glob patterns of objects
                to move, e.g. "design_* native_*"            {default: *}
    ref_chain   chain ID on the reference, if different      {default: same as chain}
    method      align, super, or cealign                     {default: align}
    ca_only     restrict the fit to CA atoms                 {default: 1}
    cycles      outlier rejection cycles, align/super only   {default: 5}
    extra       extra selection text ANDed into both sides,
                e.g. "resi 20-140"                           {default: }
    transform   0 reports RMSD without moving anything       {default: 1}
    aln_object  prefix for per pair alignment objects        {default: }
    sort        sort the printed table by RMSD               {default: 1}
    quiet       suppress the printed table                   {default: 0}

RETURNS

    dict mapping object name to RMSD.

EXAMPLES

    align_chains
    align_chains chain=B
    align_chains chain=A, ref=native, method=super
    align_chains chain=A, mobile=design_*, extra=resi 20-140
    align_chains chain=H, ref_chain=A, method=cealign, transform=0
    align_chains chain=A, aln_object=aln

    # from Python
    rms = align_chains("A", ref="native", method="super", quiet=1)
    worst = max(rms, key=rms.get)

SEE ALSO

    align, super, cealign, alignto
    """
    ca_only = _as_bool(ca_only)
    transform = _as_bool(transform)
    sort = _as_bool(sort)
    quiet = _as_bool(quiet)
    cycles = _as_int(cycles)
    extra = "" if _blank(extra) else str(extra).strip()
    aln_object = "" if _blank(aln_object) else str(aln_object).strip()
    method = str(method).strip().lower()

    if method not in ("align", "super", "cealign"):
        raise CmdException("method must be one of: align, super, cealign")

    if _blank(ref_chain):
        ref_chain = chain

    everything = cmd.get_object_list("(all)")
    if not everything:
        print(" align_chains: nothing is loaded")
        return {}

    # Resolve the reference, allowing a pattern as a convenience.
    if _blank(ref):
        ref = everything[0]
    elif ref not in everything:
        hits = [o for o in everything if fnmatch(o, ref)]
        if not hits:
            raise CmdException("no object matches ref=%s" % ref)
        ref = hits[0]

    patterns = str(mobile).split() or ["*"]
    targets = [o for o in everything
               if o != ref and any(fnmatch(o, p) for p in patterns)]
    if not targets:
        print(" align_chains: no mobile objects matched %r" % mobile)
        return {}

    def _sel(obj, ch):
        parts = ["(%s)" % obj, "polymer"]
        if not _blank(ch) and str(ch).strip().lower() not in ("*", "all"):
            parts.append("chain %s" % ch)
        if ca_only:
            parts.append("name CA")
        if extra:
            parts.append("(%s)" % extra)
        return " and ".join(parts)

    ref_sel = _sel(ref, ref_chain)
    if cmd.count_atoms(ref_sel) == 0:
        raise CmdException("reference selection is empty: %s" % ref_sel)

    results = []
    for obj in targets:
        obj_sel = _sel(obj, chain)
        if cmd.count_atoms(obj_sel) == 0:
            print(" align_chains: skipping %s, empty selection (%s)" % (obj, obj_sel))
            continue

        pair_object = "%s_%s" % (aln_object, obj) if aln_object else None

        try:
            if method == "cealign":
                out = cmd.cealign(ref_sel, obj_sel, transform=int(transform),
                                  object=pair_object)
                rmsd, natoms = out["RMSD"], out["alignment_length"]
            else:
                fn = cmd.align if method == "align" else cmd.super
                out = fn(obj_sel, ref_sel, cycles=cycles, transform=int(transform),
                         object=pair_object)
                rmsd, natoms = out[0], out[1]
        except Exception as err:  # keep going if one structure fails
            print(" align_chains: %s failed (%s)" % (obj, err))
            continue

        results.append((obj, rmsd, natoms))

    if not quiet and results:
        rows = sorted(results, key=lambda row: row[1]) if sort else results
        width = max(len(row[0]) for row in rows)
        note = "" if transform else "   [RMSD only, nothing moved]"
        print("")
        print(" reference: %s chain %s   method: %s%s"
              % (ref, ref_chain, method, note))
        print(" %s  %8s  %7s" % ("object".ljust(width), "RMSD", "atoms"))
        print(" %s" % ("-" * (width + 19)))
        for obj, rmsd, natoms in rows:
            print(" %s  %8.3f  %7d" % (obj.ljust(width), rmsd, natoms))
        print("")

    return dict((obj, rmsd) for obj, rmsd, _ in results)


cmd.extend("align_chains", align_chains)

# Tab completion for the ref argument (positional index 1).
try:
    cmd.auto_arg[1]["align_chains"] = cmd.auto_arg[0]["align"]
except Exception:
    pass
