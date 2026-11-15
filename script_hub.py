"""
script_hub.py  -  MOLtiTool: a PyMOL panel for the scripts in this folder.

Tabs
  3D Print : thicken cartoon/sticks, add struts, separate ligands/chains
             (each surfaced alone), export colored OBJ / per-color STL files,
             optionally one file per part (uses print3d.py).
  Scripts  : every .py in the folder, with the commands each one adds
             (found automatically). Pick a command, fill in its arguments,
             press Run. New scripts appear after "Refresh".

Install (once): add this line to ~/.pymolrc (or ~/.pymolrc.py as cmd.run(...)):

    run /path/to/MOLtiTool/script_hub.py

Then open it from  Plugin > Script Hub  or type  script_hub
"""

import ast
import inspect
import os
import sys
import traceback

from pymol import cmd

# Folder this file lives in (PyMOL's "run" sets __file__ to PyMOL's own
# package, but the code object knows the real path). Scripts in this folder
# are listed in the panel.
SCRIPTS_DIR = os.path.dirname(os.path.abspath(inspect.currentframe().f_code.co_filename))
if not os.path.exists(os.path.join(SCRIPTS_DIR, 'print3d.py')):
    SCRIPTS_DIR = os.path.expanduser('~/MOLtiTool')
if SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, SCRIPTS_DIR)

# files that are part of the hub itself, not user scripts
_HIDDEN = {'script_hub.py', 'print3d_core.py'}

_state = {'dock': None, 'loaded': set()}

HELP = {
    'thickness': (
        'Sets how thick the cartoon and sticks are, in Angstrom (full width and thickness, '
        'the same meaning as ChimeraX). Applies to every object. Once scaled, parts thinner '
        'than about 1.5 mm break easily on an FDM printer; the export report warns about '
        'them. Smoothness is PyMOL\'s cartoon_sampling (higher is smoother, files get bigger).'),
    'separate': (
        'Copies ligands or chains into their own objects, so each is surfaced on its own '
        '(the protein never cuts into a glycan surface) and can be exported as its own file.\n'
        'Each ligand: polymer plus one object per connected ligand. All ligands together: '
        'polymer plus one ligand object. Chains: one object per chain. Objects as they are: '
        'no copying. The original object is hidden, and struts made earlier follow the parts. '
        'Tick "One file per object" when exporting.'),
    'struts': (
        'Struts are rods that brace the model so it survives printing. They join atoms that '
        'are close in space (Max length, A) but far apart along the chain (Min loop, A). '
        'Leave Color blank and each strut takes the color of its chain (half and half if it '
        'joins two chains); type a color to use one color for all. "Add struts" replaces the '
        'previous set in that object.\n'
        'Single struts: type two atoms or residues (a residue uses its CA), or select two in '
        'the viewer and press "From sele". "Add strut" adds one, "Remove strut" removes the '
        'one between them; leave Atom 2 blank to remove every strut at Atom 1.'),
    'snapfit': (
        'For parts printed separately, such as a glycan that should click onto its protein. '
        'The guest gets a seat cut to the host\'s shape, plus a peg (the host gets a walled '
        'socket) or pockets for small magnets in both. Sizes are in mm after scaling. '
        'Clearance on "auto" is half the nozzle size (0.2 mm for a 0.4 mm nozzle). Applied '
        'when exporting with "One file per object"; needs the merge library.'),
    'export': (
        '3MF opens ready to slice in OrcaSlicer, Anycubic Slicer Next and Bambu Studio, with '
        'each part on its filament slot. OBJ keeps vertex colors; STL writes one file per '
        'filament.\n'
        'Nozzle and layer height tune the print: values left on "auto" (thinnest feature, '
        'fit and socket clearances, pillar thickness, rod smoothness) follow from them, and '
        'base sizes are rounded to whole layers. "Preview in PyMOL" shows the result over '
        'the molecule without writing files.\n'
        'Longest side sets the printed size (or use the scale). "As on screen" prints the '
        'model the way you see it, screen-down on the bed. Filament colors: every color snaps '
        'to the closest loaded filament. Single-color printer: one file per color, smaller '
        'pieces plug into sockets. Base: flat cut, or a stand with pillars. The report lists '
        'size, grams, thin parts, floating pieces and bed contact.'),
    'scripts': (
        'Every .py file in this folder, with the commands it adds. Pick a command, fill in '
        'the form and press Run; the exact command line appears in the PyMOL log. Scripts '
        'without commands get a "Run script" button. New files appear after Refresh. To '
        'write your own, see AGENTS.md in this folder.'),
}


# ---------------------------------------------------------------------------
# Script discovery (static, no code is executed until you press Load/Run)
# ---------------------------------------------------------------------------

def _doc_first_lines(doc, n=12):
    if not doc:
        return ''
    lines = inspect.cleandoc(doc).splitlines()
    return '\n'.join(lines[:n])


def scan_script(path):
    """Return dict(path, name, doc, commands=[dict(name, params, doc)]).
    params: list of (arg, default_text or None)."""
    info = dict(path=path, name=os.path.basename(path), doc='', commands=[], error=None)
    try:
        src = open(path, encoding='utf-8', errors='replace').read()
        tree = ast.parse(src)
    except Exception as e:
        info['error'] = str(e)
        return info
    info['doc'] = ast.get_docstring(tree) or ''
    if not info['doc']:
        comments = [l.lstrip('# ').rstrip() for l in src.splitlines()[:8] if l.startswith('#')]
        info['doc'] = '\n'.join(comments)
    funcs = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    found = []
    for node in ast.walk(tree):
        # cmd.extend("name", func)  or  cmd.extend(func)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) \
                and node.func.attr == 'extend' and getattr(node.func.value, 'id', '') == 'cmd':
            args = node.args
            if len(args) >= 2 and isinstance(args[0], ast.Constant) and isinstance(args[1], ast.Name):
                found.append((args[0].value, args[1].id))
            elif len(args) == 1 and isinstance(args[0], ast.Name):
                found.append((args[0].id, args[0].id))
        # @cmd.extend decorator
        if isinstance(node, ast.FunctionDef):
            for d in node.decorator_list:
                if isinstance(d, ast.Attribute) and d.attr == 'extend':
                    found.append((node.name, node.name))
    seen = set()
    for cname, fname in found:
        if cname in seen:
            continue
        seen.add(cname)
        fn = funcs.get(fname)
        params, doc = [], ''
        if fn is not None:
            doc = ast.get_docstring(fn) or ''
            a = fn.args
            pos = a.args
            defaults = [None] * (len(pos) - len(a.defaults)) + list(a.defaults)
            for arg, d in zip(pos, defaults):
                if arg.arg.startswith('_'):
                    continue
                params.append((arg.arg, None if d is None else _literal_text(d)))
            for arg, d in zip(a.kwonlyargs, a.kw_defaults):
                if not arg.arg.startswith('_'):
                    params.append((arg.arg, None if d is None else _literal_text(d)))
        info['commands'].append(dict(name=cname, params=params, doc=doc))
    return info


def _literal_text(node):
    try:
        v = ast.literal_eval(node)
        return '' if v is None else str(v)
    except Exception:
        try:
            return ast.unparse(node)
        except Exception:
            return ''


def runtime_params(func):
    params = []
    try:
        sig = inspect.signature(func)
    except (TypeError, ValueError):
        return params
    for p in sig.parameters.values():
        if p.name.startswith('_') or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD):
            continue
        d = None if p.default is p.empty else ('' if p.default is None else str(p.default))
        params.append((p.name, d))
    return params


def list_scripts(folder=None):
    folder = folder or SCRIPTS_DIR
    out = []
    try:
        names = sorted(os.listdir(folder), key=lambda s: s.lstrip('.').lower())
    except OSError:
        return out
    for n in names:
        if n.endswith('.py') and n not in _HIDDEN:
            out.append(scan_script(os.path.join(folder, n)))
    return out


def load_script(path):
    """cmd.run the file; return names of commands it newly registered."""
    before = set(cmd.keyword)
    cmd.run(path)
    _state['loaded'].add(path)
    return sorted(set(cmd.keyword) - before)


def build_command(name, values):
    """values: list of (arg, text, default). Positional for leading filled
    args, keywords after any gap."""
    parts, use_kw = [], False
    for arg, text, default in values:
        text = (text or '').strip()
        if text == '' or (default is not None and text == default):
            use_kw = True
            continue
        parts.append(('%s=%s' % (arg, text)) if use_kw else text)
    return name + (' ' + ', '.join(parts) if parts else '')


# ---------------------------------------------------------------------------
# GUI
# ---------------------------------------------------------------------------

def _qt():
    from pymol.Qt import QtWidgets, QtCore, QtGui
    return QtWidgets, QtCore, QtGui


def _enum(obj, *names):
    """Qt5/Qt6-compatible enum lookup: _enum(QtCore.Qt, 'RightDockWidgetArea',
    'DockWidgetArea.RightDockWidgetArea')."""
    for n in names:
        cur = obj
        try:
            for part in n.split('.'):
                cur = getattr(cur, part)
            return cur
        except AttributeError:
            continue
    raise AttributeError(names[0])


def make_panel(parent=None):
    QtWidgets, QtCore, QtGui = _qt()
    import print3d  # noqa: registers commands
    import print3d_core as pc

    W = QtWidgets

    class Hub(W.QWidget):
        def __init__(self):
            super().__init__(parent)
            lay = W.QVBoxLayout(self)
            lay.setContentsMargins(6, 6, 6, 6)
            self.tabs = W.QTabWidget()
            lay.addWidget(self.tabs)
            self.status = W.QLabel('')
            self.status.setWordWrap(True)
            lay.addWidget(self.status)
            self.tabs.addTab(self._print_tab(), '3D Print')
            self.tabs.addTab(self._scripts_tab(), 'Scripts')

        # ---------- helpers ----------
        def _say(self, text, error=False):
            self.status.setText(text)
            self.status.setStyleSheet('color: #c0392b;' if error else '')

        def _run(self, label, fn, *a, **k):
            try:
                r = fn(*a, **k)
                self._say(label + ' done' + ('' if r in (None, True) else ': %s' % r))
                return r
            except Exception as e:
                traceback.print_exc()
                self._say('%s failed: %s' % (label, e), error=True)

        @staticmethod
        def _shrink(w):
            w.setMinimumWidth(60)
            w.setSizePolicy(W.QSizePolicy.Policy.Expanding if hasattr(W.QSizePolicy, 'Policy')
                            else W.QSizePolicy.Expanding,
                            W.QSizePolicy.Policy.Fixed if hasattr(W.QSizePolicy, 'Policy')
                            else W.QSizePolicy.Fixed)
            return w

        def _help(self, layout, key):
            """Collapsible 'About' text (hidden until the arrow is clicked)."""
            box = W.QWidget()
            vb = W.QVBoxLayout(box)
            vb.setContentsMargins(0, 0, 0, 4)
            vb.setSpacing(2)
            btn = W.QToolButton()
            btn.setText('About')
            btn.setCheckable(True)
            btn.setAutoRaise(True)
            btn.setToolButtonStyle(_enum(QtCore.Qt, 'ToolButtonTextBesideIcon',
                                         'ToolButtonStyle.ToolButtonTextBesideIcon'))
            right = _enum(QtCore.Qt, 'RightArrow', 'ArrowType.RightArrow')
            down = _enum(QtCore.Qt, 'DownArrow', 'ArrowType.DownArrow')
            btn.setArrowType(right)
            lab = W.QLabel(HELP.get(key, ''))
            lab.setWordWrap(True)
            lab.setStyleSheet('padding: 2px 4px;')
            fnt = lab.font()
            if fnt.pointSizeF() > 0:
                fnt.setPointSizeF(fnt.pointSizeF() * 0.92)
            lab.setFont(fnt)
            lab.setVisible(False)

            def toggled(on):
                lab.setVisible(on)
                btn.setArrowType(down if on else right)
            btn.toggled.connect(toggled)
            vb.addWidget(btn)
            vb.addWidget(lab)
            if isinstance(layout, W.QFormLayout):
                layout.addRow(box)
            else:
                layout.addWidget(box)
            return box

        @staticmethod
        def _auto(spin):
            """0 shows as 'auto' (value derived from nozzle / layer height)."""
            spin.setSpecialValueText('auto')
            return spin

        @staticmethod
        def _dspin(val, lo=0.0, hi=50.0, step=0.05, dec=2):
            s = W.QDoubleSpinBox()
            s.setMinimumWidth(96)
            s.setRange(lo, hi)
            s.setSingleStep(step)
            s.setDecimals(dec)
            s.setValue(val)
            return s

        @staticmethod
        def _ispin(val, lo=0, hi=1000):
            s = W.QSpinBox()
            s.setMinimumWidth(96)
            s.setRange(lo, hi)
            s.setValue(val)
            return s

        def _color_row(self, default):
            row = W.QWidget()
            h = W.QHBoxLayout(row)
            h.setContentsMargins(0, 0, 0, 0)
            edit = W.QLineEdit(default)
            btn = W.QPushButton('Pick')
            btn.setFixedWidth(48)

            def pick():
                c = W.QColorDialog.getColor(parent=self)
                if c.isValid():
                    edit.setText('0x%02x%02x%02x' % (c.red(), c.green(), c.blue()))
            btn.clicked.connect(pick)
            h.addWidget(edit)
            h.addWidget(btn)
            return row, edit

        # ---------- 3D print tab ----------
        def _combo(self, items):
            """items: list of (value, label). Value read/written by presets."""
            c = W.QComboBox()
            c._values = [v for v, _ in items]
            c.addItems([lab for _, lab in items])
            return c

        def _reg(self, section, key, widget):
            self._fields.setdefault(section, {})[key] = widget
            return widget

        @staticmethod
        def _get(w):
            if isinstance(w, W.QComboBox):
                return w._values[w.currentIndex()] if hasattr(w, '_values') else w.currentText()
            if isinstance(w, W.QCheckBox):
                return int(w.isChecked())
            if isinstance(w, (W.QDoubleSpinBox, W.QSpinBox)):
                return w.value()
            return w.text()

        @staticmethod
        def _set(w, v):
            if isinstance(w, W.QComboBox):
                vals = getattr(w, '_values', None)
                if vals and v in vals:
                    w.setCurrentIndex(vals.index(v))
                elif not vals:
                    w.setEditText(str(v))
            elif isinstance(w, W.QCheckBox):
                w.setChecked(bool(int(float(v))) if not isinstance(v, bool) else v)
            elif isinstance(w, W.QDoubleSpinBox):
                w.setValue(float(v))
            elif isinstance(w, W.QSpinBox):
                w.setValue(int(float(v)))
            else:
                w.setText(str(v))

        def _values(self, section):
            return {k: self._get(w) for k, w in self._fields.get(section, {}).items()}

        def _print_tab(self):
            self._fields = {}
            page = W.QWidget()
            outer = W.QVBoxLayout(page)

            # presets bar
            pb = W.QHBoxLayout()
            self.preset = W.QComboBox()
            self.preset.setEditable(True)
            self.preset.lineEdit().setPlaceholderText('preset name')
            self._refresh_presets()
            for label, fn in (('Load', self._preset_load), ('Save', self._preset_save),
                              ('Delete', self._preset_delete)):
                b_ = W.QPushButton(label)
                b_.setFixedWidth(58)
                b_.clicked.connect(fn)
                pb.addWidget(b_) if label != 'Load' else None
                if label == 'Load':
                    pb.addWidget(self.preset, 1)
                    pb.addWidget(b_)
            outer.addLayout(pb)

            scroll = W.QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setHorizontalScrollBarPolicy(_enum(QtCore.Qt, 'ScrollBarAlwaysOff',
                                                      'ScrollBarPolicy.ScrollBarAlwaysOff'))
            inner = W.QWidget()
            v = W.QVBoxLayout(inner)
            scroll.setWidget(inner)
            outer.addWidget(scroll)

            # 1. thickness
            g = W.QGroupBox('1. Thickness (Angstrom)')
            f = W.QFormLayout(g)
            self._help(f, 'thickness')
            R = lambda k, w: self._reg('style', k, w)
            for lab, key, val in (('Helix width', 'helix_width', 3.0),
                                  ('Helix thickness', 'helix_thickness', 1.5),
                                  ('Strand width', 'strand_width', 3.0),
                                  ('Strand thickness', 'strand_thickness', 1.0),
                                  ('Loop thickness', 'loop_thickness', 1.25),
                                  ('Stick radius', 'stick_radius', 0.7)):
                f.addRow(lab, R(key, self._dspin(val)))
            f.addRow('Smoothness', R('sampling', self._ispin(14, 1, 40)))
            hh = R('hide_h', W.QCheckBox('Hide hydrogens'))
            hh.setChecked(True)
            f.addRow(hh)
            b = W.QPushButton('Apply thickness')
            b.clicked.connect(lambda: self._run('Thickness', print3d.printstyle,
                                                **self._values('style')))
            f.addRow(b)
            v.addWidget(g)

            # 2. separate parts
            g = W.QGroupBox('2. Separate parts (optional)')
            f = W.QFormLayout(g)
            self._help(f, 'separate')
            R = lambda k, w: self._reg('split', k, w)
            f.addRow('Split by', R('by', self._combo([
                ('ligands', 'Each ligand'), ('ligands_together', 'All ligands together'),
                ('chains', 'Chains'), ('objects', 'Objects as they are')])))
            f.addRow('Show as', R('show', self._combo([
                ('surface_ligands', 'Surface ligands'), ('surface', 'Surface every part'),
                ('keep', 'Keep look')])))
            f.addRow('Atoms', R('selection', W.QLineEdit('all')))
            prow, pcol = self._color_row('')
            pcol.setPlaceholderText('blank = keep colors')
            R('color', pcol)
            f.addRow('Surface color', prow)
            f.addRow('Surface quality', R('quality', self._ispin(1, -1, 4)))
            b = W.QPushButton('Separate')
            b.clicked.connect(self._split)
            f.addRow(b)
            v.addWidget(g)

            # 3. struts
            g = W.QGroupBox('3. Struts')
            f = W.QFormLayout(g)
            self._help(f, 'struts')
            R = lambda k, w: self._reg('struts', k, w)
            f.addRow('Atoms', R('selection', W.QLineEdit('name CA')))
            f.addRow('Max length', R('length', self._dspin(8.0, 1, 100, 0.5, 1)))
            f.addRow('Min loop', R('loop', self._dspin(60.0, 1, 500, 5, 0)))
            f.addRow('Radius', R('radius', self._dspin(0.5, 0.05, 5)))
            crow, scol = self._color_row('')
            scol.setPlaceholderText('blank = chain color')
            R('color', scol)
            f.addRow('Color', crow)
            f.addRow('Object name', R('name', W.QLineEdit('struts')))
            h_ = W.QHBoxLayout()
            b1, b2 = W.QPushButton('Add struts'), W.QPushButton('Delete struts')
            b1.clicked.connect(lambda: self._run('Struts', print3d.struts, replace=1,
                                                 **self._values('struts')))
            b2.clicked.connect(lambda: self._run('Delete struts', print3d.struts_delete,
                                                 self._values('struts')['name'] or 'struts'))
            h_.addWidget(b1)
            h_.addWidget(b2)
            f.addRow(h_)
            single = W.QLabel('Single struts')
            single.setStyleSheet('font-weight: bold; margin-top: 6px;')
            f.addRow(single)
            self.sa1, self.sa2 = W.QLineEdit(), W.QLineEdit()
            self.sa1.setPlaceholderText('e.g. chain A and resi 45')
            self.sa2.setPlaceholderText('e.g. chain A and resi 112')
            f.addRow('Atom 1', self.sa1)
            f.addRow('Atom 2', self.sa2)
            h2 = W.QHBoxLayout()
            bs = W.QPushButton('From sele')
            bs.setToolTip('Fill Atom 1 and 2 from the two atoms or residues in (sele)')
            ba, br = W.QPushButton('Add strut'), W.QPushButton('Remove strut')
            bs.clicked.connect(self._struts_from_sele)
            ba.clicked.connect(lambda: self._run(
                'Add strut', print3d.strut_add, self.sa1.text(), self.sa2.text(),
                self._values('struts')['name'] or 'struts', self._values('struts')['radius'],
                self._values('struts')['color']))
            br.clicked.connect(lambda: self._run(
                'Remove strut', print3d.strut_remove, self.sa1.text(), self.sa2.text(), ''))
            for b_ in (bs, ba, br):
                h2.addWidget(b_)
            f.addRow(h2)
            v.addWidget(g)

            # 4. snap-fit
            g = W.QGroupBox('4. Snap-fit (parts printed separately)')
            f = W.QFormLayout(g)
            self._help(f, 'snapfit')
            R = lambda k, w: self._reg('snapfit', k, w)
            self.sf_host, self.sf_guest = W.QComboBox(), W.QComboBox()
            for c in (self.sf_host, self.sf_guest):
                c.setEditable(True)
            rf = W.QPushButton('Objects')
            rf.setFixedWidth(70)
            rf.clicked.connect(self._refresh_objects)
            hr = W.QHBoxLayout()
            hr.addWidget(self.sf_host, 1)
            hr.addWidget(rf)
            f.addRow('Host (e.g. protein)', hr)
            f.addRow('Guest (e.g. glycan)', self.sf_guest)
            f.addRow('Type', R('kind', self._combo([('peg', 'Peg + socket'),
                                                     ('magnet', 'Magnet pockets')])))
            f.addRow('Peg diameter (mm)', R('diameter', self._dspin(4.0, 1, 20, 0.5, 1)))
            f.addRow('Peg length (mm)', R('length', self._dspin(6.0, 1, 40, 0.5, 1)))
            f.addRow('Magnet diameter (mm)', R('magnet_d', self._dspin(6.0, 1, 30, 0.5, 1)))
            f.addRow('Magnet height (mm)', R('magnet_h', self._dspin(3.0, 0.5, 20, 0.5, 1)))
            f.addRow('Clearance (mm)', R('clearance', self._auto(self._dspin(0.0, 0, 2, 0.05, 2))))
            f.addRow('Wall (mm)', R('wall', self._dspin(1.6, 0.4, 10, 0.2, 1)))
            h_ = W.QHBoxLayout()
            b1, b2 = W.QPushButton('Add connector'), W.QPushButton('Clear')
            b1.clicked.connect(self._add_snapfit)
            b2.clicked.connect(lambda: (print3d.snapfit_clear(), self._show_snapfits()))
            h_.addWidget(b1)
            h_.addWidget(b2)
            f.addRow(h_)
            self.sf_list = W.QLabel('')
            self.sf_list.setWordWrap(True)
            f.addRow(self.sf_list)
            v.addWidget(g)

            # 5. export
            g = W.QGroupBox('5. Export')
            f = W.QFormLayout(g)
            self._help(f, 'export')
            R = lambda k, w: self._reg('export', k, w)
            orow = W.QWidget()
            oh = W.QHBoxLayout(orow)
            oh.setContentsMargins(0, 0, 0, 0)
            self.eobj = R('objects', W.QLineEdit(''))
            self.eobj.setPlaceholderText('blank = everything shown')
            ob = W.QPushButton('Shown')
            ob.setFixedWidth(60)
            ob.clicked.connect(lambda: self.eobj.setText(
                ' '.join(cmd.get_names('objects', enabled_only=1))))
            oh.addWidget(self.eobj)
            oh.addWidget(ob)
            frow = W.QWidget()
            fh = W.QHBoxLayout(frow)
            fh.setContentsMargins(0, 0, 0, 0)
            self.efile = W.QLineEdit(os.path.expanduser('~/Desktop/model.3mf'))
            fb = W.QPushButton('...')
            fb.setFixedWidth(32)

            def browse():
                p, _ = W.QFileDialog.getSaveFileName(self, 'Export to', self.efile.text(),
                                                     '3MF (*.3mf);;OBJ (*.obj);;STL (*.stl)')
                if p:
                    self.efile.setText(p)
            fb.clicked.connect(browse)
            fh.addWidget(self.efile)
            fh.addWidget(fb)
            f.addRow('Objects', orow)
            f.addRow('File', frow)
            f.addRow('Format', R('format', self._combo([
                ('3mf', '3MF: ready to slice'), ('obj', 'OBJ with colors'),
                ('stl', 'STL per color')])))
            f.addRow('Nozzle (mm)', R('nozzle', self._dspin(0.4, 0.1, 2.0, 0.05, 2)))
            f.addRow('Layer height (mm)', R('layer_height', self._dspin(0.2, 0.04, 1.0, 0.02, 2)))
            self.eeach = R('each', W.QCheckBox('One file per object (parts line up)'))
            f.addRow(self.eeach)
            sc = R('single_color', W.QCheckBox('Single-color printer: one file per color'))
            sc.setToolTip('Each color becomes separate pieces to print in that color and glue '
                          'together; bigger pieces get sockets for the smaller ones.')
            f.addRow(sc)
            sz = R('size_mm', self._dspin(0.0, 0, 1000, 5, 0))
            sz.setSpecialValueText('off (use scale)')
            f.addRow('Longest side (mm)', sz)
            f.addRow('Scale (mm per A)', R('mm_per_A', self._dspin(1.0, 0.05, 50, 0.1, 2)))
            f.addRow('Orientation', R('orient', self._combo([
                ('view', 'As on screen (down = bed)'), ('model', 'Molecule frame')])))
            # filament slots
            self.slots = []
            for i in range(4):
                srow, sedit = self._color_row('')
                sedit.setPlaceholderText('empty')
                self.slots.append(sedit)
                f.addRow('Filament %d' % (i + 1), srow)
            f.addRow('Base', R('base', self._combo([
                ('none', 'None'), ('cut', 'Flat cut bottom'), ('stand', 'Stand with pillars')])))
            f.addRow('Cut depth (mm)', R('base_depth', self._dspin(1.0, 0.1, 20, 0.5, 1)))
            f.addRow('Stand pillars', R('pillars', self._ispin(3, 1, 12)))
            f.addRow('Pillar thickness (mm)', R('pillar_d', self._auto(self._dspin(0.0, 0, 20, 0.5, 1))))

            adv_toggle = W.QCheckBox('Show advanced options')
            adv = W.QWidget()
            af = W.QFormLayout(adv)
            af.setContentsMargins(0, 0, 0, 0)
            for key, lab, val in (('union', 'Merge each color into one solid', 1),
                                  ('attach', 'Join loose pieces with short rods', 1),
                                  ('round_ends', 'Round stick ends', 1),
                                  ('close_ends', 'Cap open ends', 1)):
                cb = R(key, W.QCheckBox(lab))
                cb.setChecked(bool(val))
                af.addRow(cb)
            af.addRow('Thinnest feature (mm)', R('min_feature', self._auto(self._dspin(0.0, 0, 10, 0.1, 2))))
            af.addRow('Join reach (mm)', R('attach_mm', self._dspin(4.0, 0.5, 30, 0.5, 1)))
            af.addRow('Socket clearance (mm)', R('socket_clearance', self._auto(self._dspin(0.0, 0, 1, 0.02, 2))))
            af.addRow('Facets per rod', R('segments', self._auto(self._ispin(0, 0, 64))))
            af.addRow('Color tolerance', R('color_tolerance', self._ispin(0, 0, 128)))
            adv.setVisible(False)
            adv_toggle.toggled.connect(adv.setVisible)
            f.addRow(adv_toggle)
            f.addRow(adv)
            hb = W.QHBoxLayout()
            bp, bc, b = (W.QPushButton('Preview in PyMOL'), W.QPushButton('Clear preview'),
                         W.QPushButton('Export'))
            bp.setToolTip('Build the print (base, connectors, rods, colors) and show it in '
                          'PyMOL over the molecule, without writing files')
            bp.clicked.connect(lambda: self._export(preview=True))
            bc.clicked.connect(lambda: (cmd.delete('preview_*'), self._say('Preview removed')))
            b.clicked.connect(self._export)
            hb.addWidget(bp)
            hb.addWidget(bc)
            f.addRow(hb)
            f.addRow(b)
            self.report = W.QPlainTextEdit()
            self.report.setReadOnly(True)
            self.report.setPlaceholderText('The printability report appears here after export.')
            self.report.setMinimumHeight(140)
            f.addRow(self.report)
            self.merge_lbl = W.QLabel()
            self.merge_lbl.setWordWrap(True)
            self.merge_btn = W.QPushButton('Install merge library (manifold3d)')
            self.merge_btn.clicked.connect(self._install_merge)
            f.addRow(self.merge_lbl)
            f.addRow(self.merge_btn)
            self._update_merge_status()
            v.addWidget(g)
            v.addStretch(1)
            for e in inner.findChildren(W.QLineEdit):
                if not isinstance(e.parent(), (W.QAbstractSpinBox, W.QComboBox)):
                    self._shrink(e)
            for c in inner.findChildren(W.QComboBox):
                self._shrink(c)
            self._refresh_objects()
            return page

        # --- presets ---
        def _refresh_presets(self):
            cur = self.preset.currentText()
            self.preset.clear()
            self.preset.addItems(sorted(print3d._load_presets()))
            self.preset.setEditText(cur)

        def _gather(self):
            d = {sec: self._values(sec) for sec in self._fields}
            d['export']['filaments'] = ' '.join(s.text().strip() for s in self.slots if s.text().strip())
            d['export'].pop('objects', None)
            return d

        def _preset_save(self):
            name = self.preset.currentText().strip()
            if not name:
                self._say('Type a preset name first', error=True)
                return
            data = print3d._load_presets()
            data[name] = self._gather()
            print3d._save_presets(data)
            self._refresh_presets()
            self.preset.setEditText(name)
            self._say('Saved preset "%s"' % name)

        def _preset_load(self):
            name = self.preset.currentText().strip()
            p = print3d._load_presets().get(name)
            if not p:
                self._say('No preset "%s"' % name, error=True)
                return
            for sec, vals in p.items():
                for k, val in vals.items():
                    w = self._fields.get(sec, {}).get(k)
                    if w is not None and k != 'objects':
                        self._set(w, val)
            fil = str(p.get('export', {}).get('filaments', '')).split()
            for i, s in enumerate(self.slots):
                s.setText(fil[i] if i < len(fil) else '')
            self._say('Loaded preset "%s" (press Apply / Export to use it)' % name)

        def _preset_delete(self):
            name = self.preset.currentText().strip()
            data = print3d._load_presets()
            if data.pop(name, None) is not None:
                print3d._save_presets(data)
                self._refresh_presets()
                self.preset.setEditText('')
                self._say('Deleted preset "%s"' % name)

        def _struts_from_sele(self):
            if 'sele' not in cmd.get_names('selections'):
                self._say('Select two atoms or residues in the viewer first', error=True)
                return
            m = cmd.get_model('sele')
            if len(m.atom) == 2:
                specs = ['%s and index %d' % (a.model, a.index) for a in m.atom]
            else:
                res = []
                for a in m.atom:
                    k = (a.model, a.chain, a.resi)
                    if k not in res:
                        res.append(k)
                if len(res) != 2:
                    self._say('(sele) has %d residues; select exactly two' % len(res), error=True)
                    return
                specs = ['%s and chain "%s" and resi %s' % k for k in res]
            self.sa1.setText(specs[0])
            self.sa2.setText(specs[1])
            self._say('Atom 1 and 2 filled from (sele)')

        # --- snap-fit ---
        def _refresh_objects(self):
            names = [n for n in cmd.get_names('objects') if n not in print3d._struts_registry]
            for c in (self.sf_host, self.sf_guest):
                cur = c.currentText()
                c.clear()
                c.addItems(names)
                if cur:
                    c.setEditText(cur)

        def _add_snapfit(self):
            vals = self._values('snapfit')
            self._run('Connector', print3d.snapfit, self.sf_host.currentText(),
                      self.sf_guest.currentText(), **vals)
            self.eeach.setChecked(True)
            self._show_snapfits()

        def _show_snapfits(self):
            lst = print3d.snapfit_list()
            self.sf_list.setText('Connectors: ' + ('; '.join(lst) if lst else 'none'))

        def _update_merge_status(self):
            ok = pc.have_manifold()
            self.merge_lbl.setText('Merge library: installed' if ok else
                                   'Merge library not installed: merging, connectors, base, '
                                   'single-color sockets and the full report need it.')
            self.merge_btn.setVisible(not ok)

        def _install_merge(self):
            self._say('Installing manifold3d (this can take a minute)...')
            W.QApplication.processEvents()
            ok = print3d.print3d_install_merge()
            if ok:
                import importlib
                importlib.invalidate_caches()
            self._update_merge_status()
            self._say('manifold3d installed' if pc.have_manifold() else
                      'Install failed; see the PyMOL log for details', error=not pc.have_manifold())

        def _split(self):
            vals = self._values('split')
            sel = vals.pop('selection') or 'all'
            parts = self._run('Separate', print3d.split_parts, sel, **vals)
            if parts:
                self._say('Parts: %s. Tip: tick "One file per object" when exporting.'
                          % ', '.join(parts))
                self.eeach.setChecked(True)
                self._refresh_objects()

        def _export(self, preview=False):
            vals = self._values('export')
            vals['preview'] = int(bool(preview))
            vals['filaments'] = ' '.join(s.text().strip() for s in self.slots if s.text().strip())
            r = self._run('Export', print3d.colorprint, self.efile.text(), **vals)
            self.report.setPlainText('\n'.join(print3d._last.get('report', [])))
            if r:
                self._say(('Preview: %s (Clear preview to remove)' if preview else 'Exported: %s') % r)

        # ---------- scripts tab ----------
        def _scripts_tab(self):
            page = W.QWidget()
            v = W.QVBoxLayout(page)
            self._help(v, 'scripts')
            top = W.QHBoxLayout()
            self.filter = W.QLineEdit()
            self.filter.setPlaceholderText('Filter scripts / commands')
            self.filter.textChanged.connect(self._apply_filter)
            rb = W.QPushButton('Refresh')
            rb.clicked.connect(self._populate)
            fb = W.QPushButton('Folder')
            fb.clicked.connect(lambda: self._open_path(SCRIPTS_DIR))
            top.addWidget(self.filter)
            top.addWidget(rb)
            top.addWidget(fb)
            v.addLayout(top)

            split = W.QSplitter(_enum(QtCore.Qt, 'Vertical', 'Orientation.Vertical'))
            self.tree = W.QTreeWidget()
            self.tree.setHeaderHidden(True)
            self.tree.currentItemChanged.connect(self._select)
            split.addWidget(self.tree)

            detail = W.QWidget()
            dv = W.QVBoxLayout(detail)
            dv.setContentsMargins(0, 0, 0, 0)
            self.title = W.QLabel('')
            self.title.setStyleSheet('font-weight: bold;')
            self.doc = W.QPlainTextEdit()
            self.doc.setReadOnly(True)
            self.doc.setMaximumHeight(120)
            self.form_box = W.QWidget()
            self.form = W.QFormLayout(self.form_box)
            form_scroll = W.QScrollArea()
            form_scroll.setWidgetResizable(True)
            form_scroll.setWidget(self.form_box)
            form_scroll.setMinimumHeight(180)
            self.cmd_preview = W.QLineEdit()
            self.cmd_preview.setReadOnly(True)
            bh = W.QHBoxLayout()
            self.run_btn = W.QPushButton('Run')
            self.run_btn.clicked.connect(self._run_selected)
            self.load_btn = W.QPushButton('Load script')
            self.load_btn.clicked.connect(self._load_selected)
            self.edit_btn = W.QPushButton('Edit')
            self.edit_btn.clicked.connect(lambda: self._open_path(self._current_file()))
            for b in (self.run_btn, self.load_btn, self.edit_btn):
                bh.addWidget(b)
            dv.addWidget(self.title)
            dv.addWidget(self.doc)
            dv.addWidget(form_scroll, 1)
            dv.addWidget(self.cmd_preview)
            dv.addLayout(bh)
            split.addWidget(detail)
            split.setSizes([220, 520])
            v.addWidget(split)
            self.inputs = []
            self._populate()
            return page

        def _populate(self):
            self.tree.clear()
            self.scripts = list_scripts()
            role = _enum(QtCore.Qt, 'UserRole', 'ItemDataRole.UserRole')
            for info in self.scripts:
                it = W.QTreeWidgetItem([info['name']])
                it.setData(0, role, ('file', info))
                first = (info['doc'] or '').strip().splitlines()
                it.setToolTip(0, first[0] if first else info['path'])
                self.tree.addTopLevelItem(it)
                for c in info['commands']:
                    ci = W.QTreeWidgetItem([c['name']])
                    ci.setData(0, role, ('cmd', info, c))
                    it.addChild(ci)
                it.setExpanded(bool(info['commands']))
            self._role = role
            self._apply_filter(self.filter.text())
            self._say('%d scripts in %s' % (len(self.scripts), SCRIPTS_DIR))

        def _apply_filter(self, text):
            t = (text or '').lower()
            for i in range(self.tree.topLevelItemCount()):
                it = self.tree.topLevelItem(i)
                kids = [it.child(j) for j in range(it.childCount())]
                kid_match = False
                for k in kids:
                    m = t in k.text(0).lower()
                    k.setHidden(bool(t) and not m)
                    kid_match |= m
                it.setHidden(bool(t) and not (t in it.text(0).lower() or kid_match))

        def _current(self):
            it = self.tree.currentItem()
            return it.data(0, self._role) if it else None

        def _current_file(self):
            d = self._current()
            return d[1]['path'] if d else SCRIPTS_DIR

        def _clear_form(self):
            while self.form.rowCount():
                self.form.removeRow(0)
            self.inputs = []

        def _select(self, *_):
            d = self._current()
            self._clear_form()
            if not d:
                return
            info = d[1]
            loaded = info['path'] in _state['loaded']
            if d[0] == 'file':
                self.title.setText(info['name'] + ('  (loaded)' if loaded else ''))
                self.doc.setPlainText(info['error'] and ('Could not read: ' + info['error'])
                                      or _doc_first_lines(info['doc'], 40))
                self.cmd_preview.setText('run ' + info['path'])
                self.run_btn.setVisible(False)
                self.load_btn.setText('Run script' if not info['commands'] else 'Load script')
                if not info['commands']:
                    self.load_btn.setToolTip('This script has no commands; running it '
                                             'executes its code right away.')
                return
            c = d[2]
            self.title.setText(c['name'] + '   (from %s)' % info['name'])
            self.doc.setPlainText(_doc_first_lines(c['doc'], 40) or
                                  _doc_first_lines(info['doc'], 20))
            for arg, default in c['params']:
                e = W.QLineEdit('' if default is None else default)
                if default is None:
                    e.setPlaceholderText('required')
                e.textChanged.connect(self._update_preview)
                self.form.addRow(arg, e)
                self.inputs.append((arg, e, default))
            self.run_btn.setVisible(True)
            self.load_btn.setText('Load script')
            self._update_preview()

        def _update_preview(self, *_):
            d = self._current()
            if d and d[0] == 'cmd':
                vals = [(a, e.text(), dflt) for a, e, dflt in self.inputs]
                self.cmd_preview.setText(build_command(d[2]['name'], vals))

        def _ensure_loaded(self, info):
            if info['path'] in _state['loaded'] and \
                    all(c['name'] in cmd.keyword for c in info['commands']):
                return
            new = load_script(info['path'])
            # pick up commands registered in ways the scanner can't see
            known = {c['name'] for c in info['commands']}
            for n in new:
                if n not in known:
                    fn = cmd.keyword[n][0]
                    info['commands'].append(dict(name=n, params=runtime_params(fn),
                                                 doc=inspect.getdoc(fn) or ''))

        def _load_selected(self):
            d = self._current()
            if not d:
                return
            info = d[1]
            try:
                self._ensure_loaded(info)
                self._say('Loaded %s' % info['name'])
                # refresh tree children for runtime-discovered commands
                it = self.tree.currentItem()
                if d[0] == 'cmd':
                    it = it.parent()
                have = {it.child(j).text(0) for j in range(it.childCount())}
                for c in info['commands']:
                    if c['name'] not in have:
                        ci = W.QTreeWidgetItem([c['name']])
                        ci.setData(0, self._role, ('cmd', info, c))
                        it.addChild(ci)
                it.setExpanded(True)
            except Exception as e:
                traceback.print_exc()
                self._say('Loading %s failed: %s' % (info['name'], e), error=True)

        def _run_selected(self):
            d = self._current()
            if not d or d[0] != 'cmd':
                return
            missing = [a for a, e, dflt in self.inputs if dflt is None and not e.text().strip()]
            if missing:
                self._say('Fill in: ' + ', '.join(missing), error=True)
                return
            try:
                self._ensure_loaded(d[1])
                line = self.cmd_preview.text()
                cmd.do(line)
                self._say('Ran: ' + line)
            except Exception as e:
                traceback.print_exc()
                self._say('Failed: %s' % e, error=True)

        def _open_path(self, path):
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(path))

    return Hub()


def script_hub():
    """
DESCRIPTION

    Open the Script Hub panel (3D printing tools + the scripts in this folder).
    """
    QtWidgets, QtCore, QtGui = _qt()
    dock = _state.get('dock')
    if dock is not None:
        try:
            dock.show()
            dock.raise_()
            return
        except RuntimeError:          # widget was deleted
            _state['dock'] = None
    win = None
    try:
        from pymol import gui
        win = gui.get_qtwindow()
    except Exception:
        pass
    panel = make_panel()
    if win is not None:
        dock = QtWidgets.QDockWidget('Script Hub', win)
        dock.setWidget(panel)
        win.addDockWidget(_enum(QtCore.Qt, 'RightDockWidgetArea',
                                'DockWidgetArea.RightDockWidgetArea'), dock)
        dock.setMinimumWidth(400)
    else:
        dock = panel
        dock.setWindowTitle('Script Hub')
    dock.show()
    _state['dock'] = dock


def __init_plugin__(app=None, _tries=[0]):
    """Add Plugin > Script Hub. Also lets the PyMOL Plugin Manager install
    this file. If the main window is not ready yet (e.g. when run from
    .pymolrc), try again shortly."""
    if _state.get('menu_added'):
        return
    try:
        from pymol.plugins import addmenuitemqt, get_pmgapp
        if get_pmgapp() is None:
            if _tries[0] < 10:
                _tries[0] += 1
                from pymol.Qt import QtCore
                QtCore.QTimer.singleShot(1000, __init_plugin__)
            return
        addmenuitemqt('Script Hub', script_hub)
        _state['menu_added'] = True
    except Exception:
        pass


cmd.extend('script_hub', script_hub)
try:
    import print3d  # noqa: makes printstyle/struts/glycan_surface/colorprint available
except Exception:
    traceback.print_exc()
__init_plugin__()
