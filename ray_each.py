"""
ray_each.py  -  Ray trace every object in the current PyMOL session, one at a time,
using the exact view you have set, and save each as <object_name>.png in the
current working directory.

Usage (in PyMOL):
    run /path/to/MOLtiTool/ray_each.py
    ray_each                          # uses current window size
    ray_each 2400, 1800, 300          # width, height, dpi
    ray_each prefix=fig_, suffix=_v1  # optional filename prefix/suffix

Notes:
  * Selections such as (sele) are skipped; only objects are rendered.
  * Group objects are skipped as images, but their members are rendered.
  * Your camera view and which objects were visible are restored afterward.
  * Check/change the output folder with:  pwd  /  cd <path>
"""

import os
import re
from pymol import cmd


def _safe_filename(name):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name)


def ray_each(width=0, height=0, dpi=-1, prefix="", suffix=""):
    width = int(width)
    height = int(height)
    dpi = float(dpi)

    # Remember the camera and what is currently shown
    view = cmd.get_view()
    previously_enabled = cmd.get_names("all", enabled_only=1)

    # Objects only (no selections), excluding group containers
    objects = [
        n for n in cmd.get_names("objects", enabled_only=0)
        if cmd.get_type(n) != "object:group"
    ]

    if not objects:
        print("ray_each: no objects found.")
        return

    outdir = os.getcwd()
    saved = []

    try:
        for name in objects:
            cmd.disable("all")
            cmd.enable(name, parents=1)
            cmd.set_view(view)  # keep your exact pose

            out = os.path.join(outdir, f"{prefix}{_safe_filename(name)}{suffix}.png")
            cmd.png(out, width=width, height=height, dpi=dpi, ray=1, quiet=1)
            saved.append(out)
            print(f"ray_each: saved {out}")
    finally:
        # Restore original visibility and view
        cmd.disable("all")
        for n in previously_enabled:
            cmd.enable(n, parents=1)
        cmd.set_view(view)

    print(f"ray_each: done, {len(saved)} image(s) written to {outdir}")


cmd.extend("ray_each", ray_each)
