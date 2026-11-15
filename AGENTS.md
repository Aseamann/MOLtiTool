# AGENTS.md: working on MOLtiTool

Guide for AI coding agents (and humans) adding tools to this repo. Read this before editing.

## Layout

```
script_hub.py      Script Hub panel (PyMOL Qt dock): "3D Print" tab + "Scripts" tab
print3d.py         PyMOL commands for 3D printing (thin layer that talks to PyMOL)
print3d_core.py    geometry engine, pure numpy (+ optional manifold3d), no PyMOL imports
*.py (others)      standalone tools; each registers commands with cmd.extend
chimerax/          ChimeraX tools (separate from the PyMOL hub)
docs/images/       README images
```

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
- `_help(layout, key)`: fold-out "About" text from the `HELP` dict at the top of the file; add one to every new section

Steps: write the command in `print3d.py` (or a new module) and register it with `cmd.extend`; add a `QGroupBox` in `_print_tab` with fields registered under a new section name; connect a button to `self._run(...)`; if the command should remember options for presets, store them in `print3d._last['section']`.

Presets (`print3d_presets.json`) are `{name: {section: {key: value}}}` and are shared by the GUI and the `print3d_preset` command.

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
- Struts are CGO objects. `print3d._struts_registry[name]` holds `dict(struts=[dict(p1, p2, m1, m2)], radius, color)`; `color=None` means each half takes its chain's color (`_chain_color_lookup`), recomputed whenever struts are drawn or exported.

## manifold3d quirks

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
