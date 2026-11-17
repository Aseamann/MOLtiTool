# AGENTS.md: working on MOLtiTool

Guide for AI coding agents (and humans) adding tools to this repo. Read this before editing.

## Layout

```
script_hub.py      Script Hub panel (PyMOL Qt dock): "3D Print" tab + "Scripts" tab
3d_printing/
  print3d.py       PyMOL commands for 3D printing (thin layer that talks to PyMOL)
  print3d_core.py  geometry engine, pure numpy (+ optional manifold3d), no PyMOL imports
  print3d_presets.json  user presets (created on first save, git-ignored)
*.py (others)      standalone tools; each registers commands with cmd.extend
chimerax/          ChimeraX tools (separate from the PyMOL hub)
docs/images/       README images
```

`script_hub.py` puts `3d_printing/` first on `sys.path` before `import print3d`; `print3d.py` adds its own folder too, so `run 3d_printing/print3d.py` also works alone. The Scripts tab lists only top-level `.py` files, so put new standalone tools in the root (or give a new subfolder its own tab).

Files listed in `.gitignore` (personal scripts, structures, presets, exports) stay local. Never add or un-ignore them.

## Adding a standalone tool (the common case)

Drop a `.py` file in the repo root. Script Hub finds it by **static parsing** (`scan_script` in `script_hub.py`, using `ast`; the file is not executed until the user presses Load/Run). It recognizes:

- `cmd.extend("name", func)` and `cmd.extend(func)`
- `@cmd.extend` decorators

For each command it builds a form from the function signature: one field per argument, defaults pre-filled, arguments without defaults marked "required", and arguments starting with `_` (such as `_self`) hidden. The function docstring is shown in the panel. Commands registered in loops or otherwise invisible to the parser are picked up after the user clicks **Run script/Load script** (`runtime_params` + diffing `cmd.keyword`).

Template:

```python
"""
my_tool.py  -  one-line summary (first docstring line is the tooltip in Script Hub).

Usage:
    run /path/to/MOLtiTool/my_tool.py
    my_tool selection [, cutoff ]
"""
from pymol import cmd


def my_tool(selection="all", cutoff=4.0, _self=cmd):
    """
DESCRIPTION

    What it does, in a sentence or two.

USAGE

    my_tool [selection [, cutoff ]]
    """
    cutoff = float(cutoff)      # command-line arguments arrive as strings
    ...


cmd.extend("my_tool", my_tool)
cmd.auto_arg[0]["my_tool"] = cmd.auto_arg[0]["zoom"]   # optional: tab-complete selections
```

Rules:
- Keep top-level code to definitions and `cmd.extend`. Anything that changes the scene at import time runs as soon as the user loads the file.
- Convert argument types yourself (`float`, `int`, booleans from "0"/"1"/"true").
- Use generic paths in docs (`/path/to/MOLtiTool/...`), never personal ones.
- Add a row to the "Other PyMOL scripts" table in README.md.
- No new hard dependencies. PyMOL ships numpy; anything else must be optional with a clear message.

## Adding a section to the 3D Print tab

All GUI code is in `make_panel()` in `script_hub.py` (class `Hub`). `_print_tab` builds numbered `QGroupBox` sections. Helpers:

- `_dspin(value, lo, hi, step, decimals)`, `_ispin(value, lo, hi)`: number fields
- `_color_row(default)` returns `(row_widget, line_edit)`, with a color picker writing `0xRRGGBB`
- `_combo([(value, label), ...])`: combo box whose stored value is `value`
- `_reg(section, key, widget)`: registers a field so presets save/load it. **`key` must equal the keyword argument of the function you call**, so `self._values(section)` can be passed straight through: `print3d.my_cmd(**self._values('my_section'))`
- `_run(label, fn, *args, **kwargs)`: calls a command, shows "done"/errors in the status line and prints tracebacks to the PyMOL log
- `_say(text, error=False)`: status line
- Reset buttons are added automatically by `_add_resets` after the tab is built: every registered field gets one, with the value it had at build time as its default. Set a widget's default before `_print_tab` returns and it just works.
- `_help(layout, key)`: fold-out "About" text from the `HELP` dict at the top of the file; add one to every new section

Steps: write the command in `print3d.py` (or a new module) and register it with `cmd.extend`; add a `QGroupBox` in `_print_tab` with fields registered under a new section name; connect a button to `self._run(...)`; if the command should remember options for presets, store them in `print3d._last['section']`.

Presets (`print3d_presets.json`) are `{name: {section: {key: value}}}` and are shared by the GUI and the `print3d_preset` command. Sections: `style` (printstyle), `surface` (surfacestyle), `struts`, `split`, `snapfit`, `export` (colorprint, including the printer fields). Built-in presets live in `BUILTIN_PRESETS` in `print3d.py`; `all_presets()` merges them with the user's file, and saving or deleting a built-in name is refused in both the command and the GUI.

Printer profiles are the `PRINTERS` dict in `print3d.py` (`bed` = build volume X, Y, Z in mm from the maker's spec sheet, `slots` = filaments the stock multicolor unit holds). Add a printer there and it appears in the GUI combo and in `print3d_printers`.

Qt: import from `pymol.Qt` (never PyQt5/PySide directly) and use `_enum(QtCore.Qt, 'Name', 'Enum.Name')` for enums, so both Qt5 and Qt6 builds of PyMOL work.

## How print3d export works

`colorprint` in `print3d.py`:
1. `_vrml_of` / `_scene_vrml`: PyMOL's `get_vrml()` with the view temporarily set to identity, so coordinates are model coordinates.
2. `pc.vrml_to_pieces`: sticks/struts (VRML cylinders) become round-ended capsules, spheres become spheres, cartoon/surface triangles are welded, re-oriented, split into connected pieces and capped.
3. `pc.make_placement` (+ `pc.view_rotation` for "as on screen"): rotate, scale to mm, put on the bed.
4. Optional `pc.snap_to_palette` to the user's filaments.
5. `pc.build_parts`: one `Part` per color, the boolean union of its pieces. Colors live in manifold3d vertex properties (channels 3-5, rgb 0-1), so booleans keep them.
6. `pc.snapfit`, `pc.attach_floating`, `pc.flat_cut` / `pc.pillar_stand`, `pc.separate_colors`.
7. `pc.analyze` + `_thin_features` give the report; `_write_group` writes 3MF (`pc.write_3mf`), OBJ or STL using `pc.part_output` (welded meshes). With `preview=1` the parts are mapped back to Angstrom (`Placement.inverse`, undoing the base's z shift) and loaded as CGO (`pc.mesh_cgo`) named `preview_*`; exports always skip `preview_*` objects.

Printer checks: with `printer=` or `bed_x/y/z`, the report compares each group's size with the build volume and its colors with `slots`; `fit=1` replaces `size_mm` with the largest scale that fits (leaving room for a stand). In preview, `_printer_overlay_cgo` draws the build volume (`preview_printer`) and `pc.thin_mask` ray-casts inward from one triangle per grid cell to find walls thinner than `min_feature` (`preview_thin`, red, lifted 0.08 mm over the part).

Supports: `pc.overhang_mask` / `pc.overhang_area` use the slicer convention (faces flatter than `support_angle` degrees from horizontal, not on the bed). `orient=min_support` scores 400 up directions with `pc.support_scores` (projected overhang area) and builds the rotation with `pc.frame_with_up`, keeping the screen's x axis. `helix_core=1` adds `_helix_cores` rods (axis from CA bisectors, radius sized to fuse with the ribbon) to the pieces before placement. `filaments=auto` uses `pc.pick_filaments` with `pc.FILAMENT_LIBRARY`.

`join_strands=1` (default) adds `_strand_joins` tubes: at every strand/loop boundary, one tube swept by `pc.sweep_tube` along the cartoon's own center line (Catmull-Rom through `_cartoon_anchor_map` points, re-centered on the cartoon mesh by `_center_on_mesh`), from one residue inside the strand to just past the loop's start. Its cross-section is an ellipse: the loop's width across, and in the strand's face normal (flattest direction of the nearby cartoon) the strand's thickness, easing to round past the arrow tip and slimming at the far end so it stays buried in the loop tube. The result is a blunt arrow tip the width of the loop that flows into the loop, as in ChimeraX. `show_strand_joins` (run by `printstyle` with `show_joins=1`) draws the same tubes as the display-only CGO `strand_joins`; `_is_overlay()` keeps it and `preview_*` out of every export.

Fine PyMOL surfaces (`surface_quality` 1 and up) have a few seams where 3 or 4 triangles share an edge; `pc.repair_nonmanifold` cuts those spots out with a ring of neighbours and refills the holes, so the surface stays one closed solid.

Nozzle and layer height (`nozzle`, `layer_height`) drive every option left at 0 ("auto"): thinnest feature, fit and socket clearances, pillar size, base sizes rounded to whole layers and rod facets (`pc.auto_segments`). New print-related options should follow the same pattern instead of hard-coding millimetres.

Keep PyMOL calls in `print3d.py` and geometry in `print3d_core.py`, so the core stays testable without PyMOL.

## PyMOL quirks (learned the hard way)

- `run file.py` sets `__file__` to PyMOL's own package. To find the script's folder use `inspect.currentframe().f_code.co_filename`.
- `cmd.get_view()[:9].reshape(3, 3)` is `R` with `camera = xyz @ R`.
- `get_vrml()` writes sticks as flat-capped cylinders (half-bonds), and cartoon plus surfaces as one triangle list whose winding alternates; use the exported normals.
- Cartoon sizes: helix width/thickness = 2 x `cartoon_oval_length/width`; strands about 1.414 x `cartoon_rect_length/width`; loops = 2 x `cartoon_loop_radius`.
- Cartoon segments take the color of their CA atom unless `cartoon_color` is set.
- `cmd.create` copies atom colors and representations, but not object-level settings. Copy them with `cmd.get_object_settings`.
- Hex colors work as `0xRRGGBB`.
- Struts live in `print3d._struts_registry[name]` = `dict(struts=[dict(p1, p2, m1, m2, l1, l2)], radius, color, labels)`; `color=None` means each half takes its chain's color (`_chain_color_lookup`), recomputed whenever struts are drawn or exported. They are drawn by `_draw_struts` as a small molecule object (strut k = residue STR k of chain S, pseudo-atoms E1/E2 of element PS, one bond, flagged `ignore`), so they can be clicked, selected and labeled. Always change struts through the registry and redraw; never edit the object directly. `p1`/`p2` are always the atom positions (used to match atoms, colors and parts); optional `d1`/`d2` are where the ends are drawn, on the cartoon center line (`_cartoon_anchor_map`: centroid of the cartoon mesh in a thin slab across each CA shown as cartoon). Use `_ends(st)` for geometry. Per-object exports build temporary CGO copies with `_strut_cgo`.

## manifold3d quirks

- Build every Manifold from arrays with `pc.make_manifold` (never `mf.Mesh(...)` directly): manifold3d rejects read-only arrays with a TypeError, and numpy 1.x (as bundled with PyMOL) returns read-only arrays from `to_mesh()`. The helper copies, reshapes and returns None for anything that isn't a valid solid.
- `colorprint ..., use_manifold=0` (GUI: "Use merge library") turns the library off for a run; every manifold-dependent step must check `solid_ok`.

- Optional: always check `pc.have_manifold()` and fall back gracefully.
- Operators: `+` union, `-` difference, `^` intersection.
- When rebuilding a Manifold from `to_mesh()` (recoloring, moving vertices), keep `merge_from_vert`/`merge_to_vert` (see `_mesh_with`, `tint`, `dilate`). Without them the mesh can fail validation and silently become `None`.

## 3MF output

OrcaSlicer, Anycubic Slicer Next and Bambu Studio read: one assembly object whose components are mesh objects (one per filament part), `Metadata/model_settings.config` with `<part id=... subtype="normal_part">` and `<metadata key="extruder" value="N"/>`, and per-triangle `paint_color` codes for mixed colors (`pc.paint_code`: 1 -> "4", 2 -> "8", 3 -> "0C", 4 -> "1C").

## Testing headless

```
pip install pymol-open-source-whl PyQt5 manifold3d
QT_QPA_PLATFORM=offscreen python test_something.py
```

```python
import os
os.environ["QT_QPA_PLATFORM"] = "offscreen"
from pymol.Qt import QtWidgets
app = QtWidgets.QApplication([])
from pymol import cmd
cmd.run("/path/to/MOLtiTool/script_hub.py")
import sys; sys.path.insert(0, "/path/to/MOLtiTool")
import script_hub
cmd.fab("AWGKLLEAAWK", "pep", ss=1)            # or cmd.load(...) a structure
panel = script_hub.make_panel()
for b in panel.findChildren(QtWidgets.QPushButton):
    if b.text() == "Apply thickness":
        b.click()
print(panel.status.text())
panel.grab().save("panel.png")                    # look at the layout
```

Before you finish: every file compiles, the panel builds without errors, your command runs from both the command line and the hub, and README.md is updated.
