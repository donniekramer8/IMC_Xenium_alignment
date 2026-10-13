"""Landmark-based alignment of Xenium DAPI to IMC nuclear (Ir191 + Ir193) in napari.

Two side-by-side canvases in one napari window: IMC (left) and Xenium DAPI (right), each
with its own rotate/flip buttons (display only). Click a landmark on the left, then the
matching one on the right. From 3 pairs a least-squares affine (DAPI -> IMC) is fitted and the
two views are linked: panning/zooming one keeps the other on the corresponding location.

Matrices returned / saved are 2x3 in (x, y) order, i.e. cv2.warpAffine / FIJI convention.

After Validate, an optional elastic step (ECC + Farneback flow, see elastic.py) refines the
mapping; IMC pixels are never resampled, the field refines where Xenium lands in native IMC px.

Run the viewer with `python run_alignment.py` (settings in config.py). Later, in analysis code:
    from imc_xenium_align import load_alignment, map_xenium_um_to_imc
    res = load_alignment('outputs/<roi>_alignment.json')      # also loads the elastic field
    xy_imc = map_xenium_um_to_imc(cells[['x_centroid', 'y_centroid']].values, res)

Don't import cv2 in the viewer process: opencv-python's bundled Qt breaks PyQt5.
"""
import json
import os

import numpy as np
import tifffile

XENIUM_UM_PER_PX = 0.2125
_P = np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]], float)  # (row, col) <-> (x, y) swap


# ----------------------------------------------------------------------------- preprocessing
def _stretch(x, lo, hi):
    return (np.clip((x.astype(np.float32) - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)


def load_imc_nuclear(path, chans=(37, 38), cofactor=20, pct=(30, 99.5)):
    """IMC DNA channels summed, arcsinh, percentile-stretched to uint8."""
    imc = tifffile.imread(path)
    nuc = np.arcsinh(imc[list(chans)].sum(0) / cofactor)
    lo, hi = np.percentile(nuc, pct)
    return _stretch(nuc, lo, hi)


def load_dapi(morph_dir, level=2, n_levels=4, pct=(30, 99.5)):
    """Xenium DAPI pyramid levels [level, level+1, ...] as uint8, all with the same lo/hi."""
    f = os.path.join(morph_dir, 'ch0000_dapi.ome.tif')
    # is_ome=False: the OME metadata references the other channel files, we only want DAPI
    with tifffile.TiffFile(f, is_ome=False) as t:
        levels = t.series[0].levels
        stop = min(level + n_levels, len(levels))
        raw = [levels[k].asarray() for k in range(level, stop)]
    lo, hi = np.percentile(raw[0][::4, ::4], pct)
    return [_stretch(r, lo, hi) for r in raw]


# ----------------------------------------------------------------------------- transforms
def _homog(L, t):
    A = np.eye(3)
    A[:2, :2] = L
    A[:2, 2] = t
    return A


def fit_affine(src, dst):
    """Least-squares affine, (N, 2) src -> dst, N >= 3. Returns 3x3."""
    X = np.c_[src, np.ones(len(src))]
    sol, *_ = np.linalg.lstsq(X, dst, rcond=None)  # (3, 2)
    A = np.eye(3)
    A[:2] = sol.T
    return A


def apply(A, pts):
    """Apply 3x3 (or 2x3) A to (N, 2) points."""
    A = np.asarray(A, float)
    return np.c_[pts, np.ones(len(pts))] @ A[:2].T


def rc_to_xy(A):
    """3x3 (row, col) matrix -> 3x3 (x, y) matrix."""
    return _P @ A @ _P


def apply_to_points(xy, M):
    """(N, 2) (x, y) points through a saved 2x3 (x, y) matrix."""
    return apply(M, np.asarray(xy, float))


def load_alignment(path):
    """Load a saved alignment JSON; matrices come back as numpy arrays. If the elastic step
    was used, its field is loaded too: res['elastic_field'] = (ecc 2x3, flow H x W x 2)."""
    with open(path) as fh:
        d = json.load(fh)
    res = {k: (np.array(v) if k.startswith('M_') or k.startswith('landmarks') else v)
           for k, v in d.items()}
    res['elastic_field'] = None
    if res.get('elastic_used'):
        z = np.load(os.path.join(os.path.dirname(os.path.abspath(path)), res['elastic_field_file']))
        res['elastic_field'] = (z['ecc_2x3_xy'], z['flow_dxdy'])
    return res


def output_matrices(A, ds):
    """3x3 (row, col) DAPI-level px -> IMC px  ->  3x3 (x, y) matrices from DAPI ds px,
    full-res px and µm."""
    M_ds = rc_to_xy(A)
    M_full = M_ds @ np.diag([1 / ds, 1 / ds, 1])
    M_um = M_full @ np.diag([1 / XENIUM_UM_PER_PX, 1 / XENIUM_UM_PER_PX, 1])
    return M_ds, M_full, M_um


# ----------------------------------------------------------------------------- elastic field
# Field (ecc 2x3, flow H x W x 2) on the IMC grid, all (x, y) IMC px:
#   backward map  B(x) = ecc · [x + flow(x); 1]  :  IMC px -> affine-aligned frame
#   image:  warped(x) = img_affine_frame(B(x))
#   points: affine-frame point q -> IMC px x solving B(x) = q (fixed-point iteration)
def _sample_flow(flow, xy):
    from scipy.ndimage import map_coordinates

    rc = [xy[:, 1], xy[:, 0]]
    return np.stack([map_coordinates(flow[..., k], rc, order=1, mode='nearest')
                     for k in (0, 1)], axis=1)


def warp_with_field(img, ecc, flow, order=1):
    """Resample an image that is already in the affine-aligned IMC frame through the field."""
    from scipy.ndimage import map_coordinates

    H, W = flow.shape[:2]
    yy, xx = np.mgrid[0:H, 0:W].astype(np.float64)
    px, py = xx + flow[..., 0], yy + flow[..., 1]
    bx = ecc[0, 0] * px + ecc[0, 1] * py + ecc[0, 2]
    by = ecc[1, 0] * px + ecc[1, 1] * py + ecc[1, 2]
    return map_coordinates(img, [by, bx], order=order, mode='constant', cval=0).astype(img.dtype)


def invert_points(q, ecc, flow, n_iter=30):
    """Affine-frame points q (N, 2 x,y) -> IMC px x with B(x) = q."""
    p = apply(np.linalg.inv(np.vstack([ecc, [0, 0, 1]])), np.asarray(q, float))
    x = p.copy()
    for _ in range(n_iter):
        x = p - _sample_flow(flow, x)
    return x


def map_xenium_um_to_imc(xy_um, res):
    """Xenium µm (x, y) -> native IMC px (x, y): affine, then the elastic field if used."""
    q = apply(res['M_dapi_um_to_imc_px'], np.asarray(xy_um, float))
    if res.get('elastic_field') is not None:
        q = invert_points(q, *res['elastic_field'])
    return q


# ----------------------------------------------------------------------------- nuclei
class XeniumNuclei:
    """Xenium nucleus_boundaries.parquet (vertices in µm), rasterised as outlines in IMC px."""

    def __init__(self, path):
        import pandas as pd

        df = pd.read_parquet(path, columns=['cell_id', 'vertex_x', 'vertex_y'])
        ids = df['cell_id'].to_numpy()
        start = np.flatnonzero(np.r_[True, ids[1:] != ids[:-1]])  # rows are grouped per nucleus
        self.xy = df[['vertex_x', 'vertex_y']].to_numpy(np.float64)
        self.start = start
        self.bounds = np.r_[start, len(ids)]
        self.n = len(start)

    def outlines(self, M_um_xy, shape, upsample=2, field=None):
        """Nuclei whose polygon lies entirely inside the IMC frame `shape`, mapped by the affine
        M_um_xy and (optionally) the elastic `field` = (ecc, flow). Returns
        (outlines, n_nuclei, filled): boolean images of shape (H*upsample, W*upsample) whose
        pixel edges align with IMC pixel edges."""
        from skimage.draw import polygon
        from skimage.segmentation import find_boundaries

        H, W = shape
        v_aff = apply(M_um_xy, self.xy)
        # only nuclei near the frame go through the (slower) elastic mapping
        m = 2.0
        if field is not None:  # how far the field can move a point: flow + ECC shift/linear part
            ecc, flow = field
            m += (np.abs(flow).max() + np.abs(ecc[:, 2]).max()
                  + 2 * np.abs(ecc[:, :2] - np.eye(2)).max() * max(H, W))
        near = ((v_aff[:, 0] >= -m) & (v_aff[:, 0] <= W + m) &
                (v_aff[:, 1] >= -m) & (v_aff[:, 1] <= H + m))
        cand = np.flatnonzero(np.logical_or.reduceat(near, self.start))
        counts = np.diff(self.bounds)[cand]
        local_start = np.r_[0, np.cumsum(counts)[:-1]].astype(int)
        idx = np.repeat(self.bounds[cand] - local_start, counts) + np.arange(counts.sum())
        v = v_aff[idx]
        if field is not None:
            v = invert_points(v, *field)
        inside = ((v[:, 0] >= -0.5) & (v[:, 0] <= W - 0.5) &
                  (v[:, 1] >= -0.5) & (v[:, 1] <= H - 0.5))
        keep = np.flatnonzero(np.logical_and.reduceat(inside, local_start)) if len(cand) else []
        ends = np.r_[local_start[1:], len(v)]
        lab = np.zeros((H * upsample, W * upsample), np.int32)
        for k, i in enumerate(keep, 1):
            vv = (v[local_start[i]:ends[i]] + 0.5) * upsample - 0.5
            rr, cc = polygon(vv[:, 1], vv[:, 0], lab.shape)
            lab[rr, cc] = k
        return find_boundaries(lab, mode='inner'), len(keep), lab > 0


def nuclei_iou(imc, filled, upsample=2):
    """IoU between Otsu-thresholded IMC nuclear signal and the filled Xenium nucleus mask
    (Otsu-thresholded nuclei overlap, as used for H&E/IF tile QC)."""
    from skimage.filters import threshold_otsu

    H, W = imc.shape
    xen = filled.reshape(H, upsample, W, upsample).mean((1, 3)) >= 0.5
    imc_m = imc > threshold_otsu(imc)
    union = np.logical_or(imc_m, xen).sum()
    return float(np.logical_and(imc_m, xen).sum() / union) if union else 0.0


# ----------------------------------------------------------------------------- viewer
_ROT_CW = _homog([[0, 1], [-1, 0]], (0, 0))   # 90° clockwise on screen, in (row, col)
_FLIP_H = _homog([[1, 0], [0, -1]], (0, 0))
_FLIP_V = _homog([[-1, 0], [0, 1]], (0, 0))
_TEXT = {'string': '{idx}', 'size': 14, 'color': 'yellow', 'anchor': 'upper_left'}


class _Panel:
    """One canvas: greyscale image + its landmarks. Orientation (rot/flip) is display-only;
    landmarks are stored in native image pixels."""

    def __init__(self, model, image, name, multiscale=False):
        self.model = model
        self.O = np.eye(3)  # display orientation: native px -> world
        self.extra = []     # (layer, base affine) shown under the same orientation
        self.img = model.add_image(image, multiscale=multiscale, name=name, colormap='gray',
                                   contrast_limits=(0, 255))
        # canvas_size_limits caps the on-screen marker size so it stays small when zoomed in
        self.pts = model.add_points(np.empty((0, 2)), name=f'{name} landmarks', size=12,
                                    canvas_size_limits=(6, 12), face_color='transparent',
                                    border_color='yellow', border_width=0.15,
                                    features={'idx': np.empty(0, int)}, text=_TEXT)
        self.activate_points()

    def activate_points(self):
        self.model.layers.selection.active = self.pts
        self.pts.mode = 'add'  # click = add point, drag = pan

    def add_overlay(self, layer, base):
        """Put `layer` (already in model) just above the image, with native-px affine `base`."""
        self.extra.append((layer, base))
        layer.affine = self.O @ base
        self.model.layers.move(self.model.layers.index(layer), 1)
        self.activate_points()

    def center_data(self):
        """Native-pixel position at the centre of the view."""
        c = np.array(self.model.camera.center[-2:], float)
        return apply(np.linalg.inv(self.O), c[None])[0]

    def look_at(self, d):
        self.model.camera.center = tuple(apply(self.O, np.asarray(d, float)[None])[0])

    def orient(self, M):
        d = self.center_data()
        self.O = M @ self.O
        self.img.affine = self.O
        self.pts.affine = self.O
        for layer, base in self.extra:
            layer.affine = self.O @ base
        self.look_at(d)


class AlignmentViewer:
    UPSAMPLE = 2    # nucleus outlines drawn at 2x IMC resolution (thinner lines)
    MIN_PAIRS = 6   # pairs needed before Validate

    def __init__(self, imc_path, morph_dir, dapi_level=2, out_path=None,
                 nuclei_path=None, imc_kwargs=None, dapi_kwargs=None, elastic_defaults=None):
        import napari
        from napari._qt.qt_viewer import QtViewer
        from napari.components import ViewerModel

        from elastic import DEFAULTS as ELASTIC_DEFAULTS

        self.imc_path = imc_path
        self.dapi_level = dapi_level
        self.ds = 2 ** dapi_level
        roi = os.path.splitext(os.path.basename(imc_path))[0]
        self.out_path = out_path or f'{roi}_alignment.json'
        self.result = None
        self.A = None        # DAPI ds px (row, col) -> IMC px (row, col), from >=3 pairs
        self._validated_A = None  # fit that was validated; Save only allowed for this fit
        self._syncing = False
        self.el_defaults = {**ELASTIC_DEFAULTS, **(elastic_defaults or {})}
        self.elastic = None  # result of the last elastic run for the validated fit
        self._el_job = None  # running elastic subprocess
        self._overlays = {}  # key -> napari layer on the IMC panel

        print('loading IMC ...')
        self.imc = load_imc_nuclear(imc_path, **(imc_kwargs or {}))
        print('loading DAPI pyramid ...')
        self.dapi = load_dapi(morph_dir, level=dapi_level, **(dapi_kwargs or {}))
        self.nuclei = None
        if nuclei_path:
            print('loading Xenium nucleus boundaries ...')
            self.nuclei = XeniumNuclei(nuclei_path)

        # left = napari's own viewer, right = second canvas embedded in the same window
        self.viewer = napari.Viewer(title=f'IMC ← Xenium DAPI alignment: {roi}')
        self._dapi_model = ViewerModel(title='Xenium DAPI')
        self._qt_dapi = QtViewer(self._dapi_model)  # create the canvas before adding layers
        self.L = _Panel(self.viewer, self.imc, 'IMC')
        self.R = _Panel(self._dapi_model, self.dapi, 'DAPI', multiscale=True)
        self._build_window()

        for p in (self.L, self.R):
            p.pts.events.data.connect(self._on_points)
        for src, dst in ((self.L, self.R), (self.R, self.L)):
            src.model.camera.events.center.connect(lambda e, s=src, d=dst: self._sync(s, d))
            src.model.camera.events.zoom.connect(lambda e, s=src, d=dst: self._sync(s, d))

        self._update()
        self.viewer.reset_view()
        self._dapi_model.reset_view()

    # ---- window --------------------------------------------------------------------
    def _build_window(self):
        from qtpy.QtCore import Qt, QTimer
        from qtpy.QtWidgets import (QApplication, QCheckBox, QDoubleSpinBox, QHBoxLayout, QLabel,
                                    QPushButton, QSpinBox, QSplitter, QVBoxLayout, QWidget)

        qt_main = self.viewer.window._qt_window
        qt_left = self.viewer.window._qt_viewer
        for dock in (qt_left.dockLayerList, qt_left.dockLayerControls):
            qt_main.removeDockWidget(dock)
            dock.hide()

        def column(title, panel, canvas):
            w = QWidget()
            v = QVBoxLayout(w)
            v.setContentsMargins(2, 2, 2, 2)
            bar = QHBoxLayout()
            lbl = QLabel(f'<b>{title}</b>')
            bar.addWidget(lbl)
            bar.addStretch()
            for txt, M in (('⟲ 90°', np.linalg.inv(_ROT_CW)), ('⟳ 90°', _ROT_CW),
                           ('Flip ↔', _FLIP_H), ('Flip ↕', _FLIP_V)):
                b = QPushButton(txt)
                b.clicked.connect(lambda _=False, M=M: panel.orient(M))
                bar.addWidget(b)
            v.addLayout(bar)
            v.addWidget(canvas, 1)
            return w

        split = QSplitter(Qt.Orientation.Horizontal)
        split.addWidget(column('IMC nuclear', self.L, qt_left))   # reparents napari's canvas
        split.addWidget(column('Xenium DAPI', self.R, self._qt_dapi))
        split.setSizes([1, 1])

        # row 1: landmarks / affine
        bottom = QHBoxLayout()
        self._status = QLabel()
        bottom.addWidget(self._status, 1)
        self._link = QCheckBox('Link views')
        self._link.setToolTip('Available from 3 pairs: panning/zooming one view follows in the other')
        self._link.toggled.connect(lambda on: on and self._sync(self.L, self.R))
        bottom.addWidget(self._link)
        self._del_idx = QSpinBox(minimum=1, maximum=1)
        b_del = QPushButton('Delete pair #')
        b_del.clicked.connect(lambda: self.delete_pair(self._del_idx.value() - 1))
        b_last = QPushButton('Delete last pair')
        b_last.clicked.connect(lambda: self.delete_pair(-1))
        self._btn_validate = QPushButton('Validate')
        self._btn_validate.setToolTip(f'Available from {self.MIN_PAIRS} pairs: transform the Xenium '
                                      'nuclei (fully inside IMC) onto the IMC image, in cyan')
        self._btn_validate.clicked.connect(self.validate)
        self._btn_save = QPushButton('Save')
        self._btn_save.setToolTip('Available after Validate (overwrites the previous result); '
                                  'includes the elastic field if "use in Save" is ticked')
        self._btn_save.clicked.connect(self.save)
        for wdg in (b_del, self._del_idx, b_last, self._btn_validate, self._btn_save):
            bottom.addWidget(wdg)

        # row 2: optional elastic refinement (after Validate)
        self._el_row = QWidget()
        er = QHBoxLayout(self._el_row)
        er.setContentsMargins(0, 0, 0, 0)
        er.addWidget(QLabel('<b>Elastic (optional):</b>'))
        d = self.el_defaults
        self._el_ecc = QCheckBox('ECC')
        self._el_ecc.setToolTip('ECC affine refinement before the optical flow')
        self._el_ecc.setChecked(bool(d['ecc']))
        er.addWidget(self._el_ecc)
        self._el_w = {}
        for key, lo, hi, step, dec, tip in (
                ('levels', 1, 8, 1, None, 'Farneback pyramid levels'),
                ('winsize', 3, 101, 2, None, 'Farneback averaging window (px): larger = smoother, '
                                             'more robust, less local detail'),
                ('iterations', 1, 20, 1, None, 'Farneback iterations per level'),
                ('poly_n', 5, 7, 2, None, 'Farneback polynomial neighbourhood (5 or 7)'),
                ('poly_sigma', 0.5, 3.0, 0.1, 1, 'Farneback polynomial sigma (1.1 for poly_n 5, '
                                                 '1.5 for 7)'),
                ('smooth_sigma', 0.0, 50.0, 1.0, 1, 'Gaussian smoothing of the flow field (IMC px); '
                                                    '0 = off'),
                ('strength', 0.0, 2.0, 0.1, 2, 'Scale applied to the flow (1.0 = full correction)')):
            w = QSpinBox() if dec is None else QDoubleSpinBox(decimals=dec)
            w.setRange(lo, hi)
            w.setSingleStep(step)
            w.setValue(d[key])
            w.setPrefix(f'{key} ')
            w.setToolTip(tip)
            self._el_w[key] = w
            er.addWidget(w)
        self._btn_elastic = QPushButton('Run elastic')
        self._btn_elastic.setToolTip('Compute (or re-compute with the current parameters)')
        self._btn_elastic.clicked.connect(self.run_elastic)
        er.addWidget(self._btn_elastic)
        er.addWidget(QLabel('  show:'))
        self._show_aff = QCheckBox('affine (cyan)')
        self._show_el = QCheckBox('elastic (orange)')
        self._show_dapi = QCheckBox('DAPI overlay')
        self._show_dapi.setToolTip('Xenium DAPI warped into IMC (green); elastic version if '
                                   '"use in Save" is ticked, else affine')
        self._use_el = QCheckBox('use in Save')
        self._use_el.setToolTip('Save the affine + elastic registration (untick to discard elastic)')
        self._show_aff.toggled.connect(lambda on: self._set_visible('nuc_aff', on))
        self._show_el.toggled.connect(lambda on: self._set_visible('nuc_el', on))
        self._show_dapi.toggled.connect(lambda _: self._update_dapi_overlay())
        self._use_el.toggled.connect(lambda _: (self._update_dapi_overlay(), self._refresh_status()))
        for wdg in (self._show_aff, self._show_el, self._show_dapi, self._use_el):
            er.addWidget(wdg)
        self._el_status = QLabel()
        er.addWidget(self._el_status, 1)

        root = QWidget()
        v = QVBoxLayout(root)
        v.setContentsMargins(4, 0, 4, 4)
        v.addWidget(split, 1)
        v.addLayout(bottom)
        v.addWidget(self._el_row)
        qt_main.setCentralWidget(root)
        root.show()  # window is already shown, so the new central widget must be shown explicitly
        qt_main.resize(1800, 1000)

        self._el_timer = QTimer()
        self._el_timer.setInterval(150)
        self._el_timer.timeout.connect(self._poll_elastic)
        # close the embedded canvas with the window (otherwise vispy crashes at interpreter exit)
        self.viewer.window._qt_viewer.destroyed.connect(self._qt_dapi.close)
        QApplication.instance().aboutToQuit.connect(self._qt_dapi.close)
        QApplication.instance().aboutToQuit.connect(self._kill_elastic)

    # ---- overlays on the IMC panel -------------------------------------------------
    def _set_overlay(self, key, data, colormap, upsampled, opacity=1.0):
        layer = self._overlays.get(key)
        if layer is None:
            layer = self.viewer.add_image(data, name=key, colormap=colormap, blending='additive',
                                          contrast_limits=(0, 1 if data.max() <= 1 else 255),
                                          opacity=opacity)
            u = self.UPSAMPLE if upsampled else 1
            # upsampled px j -> IMC px (j + 0.5) / u - 0.5
            self.L.add_overlay(layer, _homog(np.eye(2) / u, (0.5 / u - 0.5,) * 2))
            self._overlays[key] = layer
        else:
            layer.data = data
        return layer

    def _set_visible(self, key, on):
        if key in self._overlays:
            self._overlays[key].visible = on

    def _update_dapi_overlay(self):
        if not self.validated:
            return
        use_el = self.elastic is not None and self._use_el.isChecked()
        data = self.elastic['dapi'] if use_el else self._dapi_aff
        self._set_overlay('dapi', data, 'green', upsampled=False, opacity=0.7)
        self._set_visible('dapi', self._show_dapi.isChecked())

    # ---- landmarks / fit -----------------------------------------------------------
    def pairs(self):
        n = min(len(self.L.pts.data), len(self.R.pts.data))
        return np.asarray(self.R.pts.data[:n], float), np.asarray(self.L.pts.data[:n], float)

    def residuals(self):
        src, dst = self.pairs()
        return np.linalg.norm(apply(self.A, src) - dst, axis=1)

    def _update(self):
        src, dst = self.pairs()
        had = self.A is not None
        self.A = fit_affine(src, dst) if len(src) >= 3 else None
        self._link.setEnabled(self.A is not None)
        if self.A is None:
            self._link.setChecked(False)
        elif not had:
            self._link.setChecked(True)  # turn linking on the first time it becomes available
        if self._validated_A is not None and not (
                self.A is not None and np.array_equal(self.A, self._validated_A)):
            self._invalidate()
        self._refresh_status()

    @property
    def validated(self):
        return self._validated_A is not None

    def validate(self):
        """Overlay the Xenium nuclei transformed with the current fit; unlocks Save and the
        optional elastic step."""
        from scipy.ndimage import affine_transform

        if self.A is None or len(self.pairs()[0]) < self.MIN_PAIRS:
            return
        Ai = np.linalg.inv(self.A)  # IMC (row, col) -> DAPI ds (row, col)
        self._dapi_aff = affine_transform(self.dapi[0], Ai[:2, :2], offset=Ai[:2, 2],
                                          output_shape=self.imc.shape, order=1)
        self._aff_nuc = None
        if self.nuclei is not None:
            edges, n, filled = self.nuclei_outlines()
            self._aff_nuc = dict(n=n, iou=nuclei_iou(self.imc, filled, self.UPSAMPLE))
            self._set_overlay('nuc_aff', edges.astype(np.uint8), 'cyan', upsampled=True)
            self._show_aff.setChecked(True)
            self._set_visible('nuc_aff', True)
        self._validated_A = self.A.copy()
        self.elastic = None
        self._el_row.setEnabled(True)
        self._show_el.setEnabled(False)
        self._use_el.setEnabled(False)
        self._update_dapi_overlay()
        self._refresh_status()

    def _invalidate(self):
        self._validated_A = None
        self.elastic = None
        self._kill_elastic()
        for layer in self._overlays.values():
            layer.visible = False
        for cb in (self._show_aff, self._show_el, self._use_el):
            cb.setChecked(False)
        self._el_row.setEnabled(False)
        self._el_status.setText('')

    def nuclei_outlines(self, field=None):
        M_um = output_matrices(self._validated_A if self.validated else self.A, self.ds)[2]
        return self.nuclei.outlines(M_um, self.imc.shape, self.UPSAMPLE, field=field)

    # ---- elastic -------------------------------------------------------------------
    def elastic_params(self):
        p = {k: w.value() for k, w in self._el_w.items()}
        p['ecc'] = self._el_ecc.isChecked()
        return {**self.el_defaults, **p}

    def run_elastic(self):
        """Start the elastic registration (IMC fixed, affine-warped DAPI moving) in a separate
        process (OpenCV must not be loaded into the viewer process)."""
        import subprocess
        import sys
        import tempfile

        if not self.validated or self._el_job is not None:
            return
        params = self.elastic_params()
        tmp = tempfile.mkdtemp(prefix='imc_elastic_')
        src, dst = os.path.join(tmp, 'in.npz'), os.path.join(tmp, 'out.npz')
        np.savez(src, fixed=self.imc, moving=self._dapi_aff, params=json.dumps(params))
        script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'elastic.py')
        log = open(os.path.join(tmp, 'log.txt'), 'w')
        proc = subprocess.Popen([sys.executable, script, src, dst], stdout=log,
                                stderr=subprocess.STDOUT)
        self._el_job = dict(proc=proc, tmp=tmp, out=dst, log=log, params=params,
                            A=self._validated_A.copy())
        self._btn_elastic.setEnabled(False)
        self._btn_elastic.setText('Running …')
        self._el_status.setText('running elastic …')
        self._el_timer.start()

    def _kill_elastic(self):
        import shutil

        job, self._el_job = self._el_job, None
        if job is None:
            return
        self._el_timer.stop()
        if job['proc'].poll() is None:
            job['proc'].kill()
            job['proc'].wait()
        job['log'].close()
        shutil.rmtree(job['tmp'], ignore_errors=True)
        self._btn_elastic.setEnabled(True)
        self._btn_elastic.setText('Run elastic')

    def _poll_elastic(self):
        job = self._el_job
        if job is None or job['proc'].poll() is None:
            return
        rc = job['proc'].returncode
        job['log'].close()
        msg = open(os.path.join(job['tmp'], 'log.txt')).read().strip().splitlines()
        out = None
        if rc == 0 and self.validated and np.array_equal(job['A'], self._validated_A):
            z = np.load(job['out'])
            out = dict(ecc=z['ecc'], flow=z['flow'], ecc_ok=bool(z['ecc_ok']))
        self._kill_elastic()
        if rc != 0:
            self._el_status.setText(f'elastic failed: {msg[-1] if msg else rc}')
            return
        if out is None:  # landmarks changed while it was running
            return
        self._accept_elastic_run(out, job['params'])

    def _accept_elastic_run(self, out, params):
        field = (out['ecc'], out['flow'])
        el = dict(params=params, ecc=out['ecc'], flow=out['flow'], ecc_ok=out['ecc_ok'],
                  dapi=warp_with_field(self._dapi_aff, *field),
                  max_disp=float(np.linalg.norm(out['flow'], axis=-1).max()))
        if self.nuclei is not None:
            edges, n, filled = self.nuclei_outlines(field=field)
            el.update(n=n, iou=nuclei_iou(self.imc, filled, self.UPSAMPLE))
            self._set_overlay('nuc_el', edges.astype(np.uint8), 'bop orange', upsampled=True)
        self.elastic = el
        self._show_el.setEnabled(True)
        self._show_el.setChecked(True)
        self._set_visible('nuc_el', True)
        self._use_el.setEnabled(True)
        self._use_el.setChecked(True)
        self._update_dapi_overlay()
        self._refresh_status()

    # ---- navigation ----------------------------------------------------------------
    def _on_points(self, event):
        layer = event.source
        layer.features = {'idx': np.arange(1, len(layer.data) + 1)}
        layer.refresh_text()
        action = getattr(event, 'action', None)
        if action in ('changing', 'adding', 'removing'):
            return
        self._update()
        if action == 'added' and self.A is not None:
            # jump the other view to where the current fit predicts the matching point
            src, dst = (self.L, self.R) if layer is self.L.pts else (self.R, self.L)
            if len(dst.pts.data) < len(src.pts.data):
                M = np.linalg.inv(self.A) if src is self.L else self.A
                self._syncing = True
                try:
                    dst.look_at(apply(M, layer.data[-1:])[0])
                finally:
                    self._syncing = False

    def _sync(self, src, dst):
        if self._syncing or self.A is None or not self._link.isChecked():
            return
        M = np.linalg.inv(self.A) if src is self.L else self.A  # src native -> dst native
        self._syncing = True
        try:
            dst.look_at(apply(M, src.center_data()[None])[0])
            s = np.sqrt(abs(np.linalg.det(M[:2, :2])))           # dst px per src px
            dst.model.camera.zoom = src.model.camera.zoom / s
        finally:
            self._syncing = False

    def delete_pair(self, i):
        for p in (self.L, self.R):
            n = len(p.pts.data)
            j = i if i >= 0 else n + i
            if 0 <= j < n:
                p.pts.data = np.delete(p.pts.data, j, axis=0)
        self._update()

    def _refresh_status(self):
        n = len(self.pairs()[0])
        nl, nr = len(self.L.pts.data), len(self.R.pts.data)
        s = f'pairs: {n}'
        if nl != nr:
            s += f'   ·   waiting for {"DAPI (right)" if nl > nr else "IMC (left)"} point #{max(nl, nr)}'
        if n < 3:
            s += f'   ·   {3 - n} more to fit, {self.MIN_PAIRS} to validate (spread them out)'
        tip = ''
        if self.A is not None:
            r = self.residuals()
            s += f'   ·   affine RMSE {np.sqrt((r ** 2).mean()):.2f} IMC px'
            if n >= 4:
                k = int(r.argmax())
                s += f'   ·   worst #{k + 1}: {r[k]:.2f} px'
            tip = '\n'.join(f'#{i + 1}: {e:.2f} px' for i, e in enumerate(r))
        if self.validated:
            if self._aff_nuc is not None:
                s += (f'   ·   VALIDATED: {self._aff_nuc["n"]} Xenium nuclei (cyan), '
                      f'IoU {self._aff_nuc["iou"]:.3f} – Save, or try elastic below')
            else:
                s += '   ·   VALIDATED – Save, or try elastic below'
        elif n >= self.MIN_PAIRS:
            s += '   ·   ready → Validate'
        elif n >= 3:
            s += f'   ·   {self.MIN_PAIRS - n} more for Validate'
        self._status.setText(s)
        self._status.setToolTip(tip)
        self._del_idx.setMaximum(max(1, max(nl, nr)))
        self._btn_validate.setEnabled(self.A is not None and n >= self.MIN_PAIRS)
        self._btn_save.setEnabled(self.validated)
        el = self.elastic
        if el is not None:
            e = f'max |flow| {el["max_disp"]:.1f} px' + ('' if el['ecc_ok'] or not el['params']['ecc']
                                                       else ' · ECC did not converge (identity)')
            if 'iou' in el:
                e = f'IoU affine {self._aff_nuc["iou"]:.3f} → elastic {el["iou"]:.3f} · ' + e
            e += '  →  Save: ' + ('affine + elastic' if self._use_el.isChecked() else 'affine only')
            self._el_status.setText(e)
        if not self.validated:
            self._el_row.setEnabled(False)

    # ---- output --------------------------------------------------------------------
    def use_elastic(self):
        return self.elastic is not None and self._use_el.isChecked()

    def field(self):
        return (self.elastic['ecc'], self.elastic['flow']) if self.use_elastic() else None

    def compute_result(self):
        src, dst = self.pairs()
        if len(src) < 3:
            raise RuntimeError('need at least 3 landmark pairs')
        A = fit_affine(src, dst)
        r = np.linalg.norm(apply(A, src) - dst, axis=1)
        M_ds, M_full, M_um = output_matrices(A, self.ds)  # DAPI ds px / full-res px / µm -> IMC px
        xy = lambda p: p[:, ::-1]  # (row, col) -> (x, y)
        res = {
            'imc_path': self.imc_path,
            'dapi_level': self.dapi_level,
            'xenium_um_per_px': XENIUM_UM_PER_PX,
            'n_pairs': int(len(src)),
            'rmse_imc_px': float(np.sqrt((r ** 2).mean())),
            'residuals_imc_px': r.tolist(),
            'landmarks_imc_px_xy': xy(dst).tolist(),
            'landmarks_dapi_ds_px_xy': xy(src).tolist(),
            'landmarks_dapi_fullres_px_xy': (xy(src) * self.ds).tolist(),
            'M_dapi_fullres_px_to_imc_px': M_full[:2].tolist(),
            'M_dapi_um_to_imc_px': M_um[:2].tolist(),
            'M_dapi_ds_px_to_imc_px': M_ds[:2].tolist(),
            'M_imc_px_to_dapi_fullres_px': np.linalg.inv(M_full)[:2].tolist(),
            'M_imc_px_to_dapi_um': np.linalg.inv(M_um)[:2].tolist(),
            'iou_affine': self._aff_nuc['iou'] if self._aff_nuc else None,
            'elastic_used': self.use_elastic(),
        }
        if self.use_elastic():
            el = self.elastic
            res.update({
                'elastic_params': el['params'],
                'elastic_ecc_converged': el['ecc_ok'],
                'elastic_max_flow_px': el['max_disp'],
                'iou_elastic': el.get('iou'),
                'elastic_field_file': os.path.basename(self._field_path()),
                'elastic_definition': ELASTIC_DEFINITION,
            })
        return res

    def _field_path(self, path=None):
        return os.path.splitext(path or self.out_path)[0].removesuffix('_alignment') + \
            '_elastic_field.npz'

    def save_qc(self, path):
        """QC PNG: IMC, DAPI warped into IMC, overlay (IMC magenta / DAPI green, landmarks
        yellow), Xenium nuclei with the affine (cyan) and, if used, the elastic (orange) fit."""
        from matplotlib.colors import ListedColormap
        from matplotlib.figure import Figure

        H, W = self.imc.shape
        use_el = self.use_elastic()
        warped = self.elastic['dapi'] if use_el else self._dapi_aff
        dst = self.pairs()[1]
        r = self.residuals()
        panels = ['imc', 'dapi', 'overlay']
        if self.nuclei is not None:
            panels += ['nuc_aff'] + (['nuc_el'] if use_el else [])
        ncol = 3
        nrow = int(np.ceil(len(panels) / ncol))
        fig = Figure(figsize=(6 * ncol, 6.3 * nrow))
        axs = np.atleast_1d(fig.subplots(nrow, ncol)).ravel()
        tag = 'affine + elastic' if use_el else 'affine'

        def outlines(a, edges, color):
            a.imshow(self.imc, cmap='gray')
            a.imshow(np.ma.masked_where(~edges, edges), cmap=ListedColormap([color]),
                     extent=(-0.5, W - 0.5, H - 0.5, -0.5), interpolation='nearest')

        for a, p in zip(axs, panels):
            if p == 'imc':
                a.imshow(self.imc, cmap='gray')
                a.set_title('IMC nuclear')
            elif p == 'dapi':
                a.imshow(warped, cmap='gray')
                a.set_title(f'Xenium DAPI → IMC ({tag})')
            elif p == 'overlay':
                a.imshow(np.dstack([self.imc, warped, self.imc]))
                a.set_title(f'overlay ({tag}), {len(r)} pairs, RMSE {np.sqrt((r ** 2).mean()):.2f} px')
                a.scatter(dst[:, 1], dst[:, 0], s=40, facecolors='none', edgecolors='yellow')
                for i, (rr, cc) in enumerate(dst):
                    a.text(cc + 8, rr, str(i + 1), color='yellow', fontsize=9, clip_on=True)
                a.set_xlim(-0.5, W - 0.5)
                a.set_ylim(H - 0.5, -0.5)
            elif p == 'nuc_aff':
                edges, n, _ = self.nuclei_outlines()
                outlines(a, edges, 'cyan')
                a.set_title(f'affine: {n} Xenium nuclei, IoU {self._aff_nuc["iou"]:.3f}')
            elif p == 'nuc_el':
                edges, n, _ = self.nuclei_outlines(field=self.field())
                outlines(a, edges, 'orange')
                a.set_title(f'affine + elastic: {n} Xenium nuclei, IoU {self.elastic["iou"]:.3f}')
        for a in axs:
            a.axis('off')
        fig.tight_layout()
        fig.savefig(path, dpi=150)

    def save(self, path=None):
        """Write JSON (+ elastic field npz if used) + QC PNG for the validated fit; overwrites
        previous results for this ROI."""
        if not self.validated:
            raise RuntimeError('Validate the fit before saving')
        path = path or self.out_path
        res = self.compute_result()
        fpath = self._field_path(path)
        if res['elastic_used']:
            np.savez_compressed(fpath, ecc_2x3_xy=self.elastic['ecc'],
                                flow_dxdy=self.elastic['flow'].astype(np.float32),
                                definition=ELASTIC_DEFINITION)
        elif os.path.exists(fpath):
            os.remove(fpath)  # stale field from an earlier save
        with open(path, 'w') as fh:
            json.dump(res, fh, indent=2)
        self.result = load_alignment(path)
        qc = os.path.splitext(path)[0] + '_qc.png'
        self.save_qc(qc)
        msg = (f'saved {res["n_pairs"]} pairs, RMSE {res["rmse_imc_px"]:.2f} px, '
               f'{"affine + elastic" if res["elastic_used"] else "affine only"} → {path} (+ _qc.png'
               + (', _elastic_field.npz)' if res['elastic_used'] else ')'))
        self.viewer.status = msg
        self._status.setText(msg)
        print(msg)
        return self.result


ELASTIC_DEFINITION = (
    'Field on the native IMC pixel grid, (x, y) = (col, row) IMC px. Backward map '
    'B(x) = ecc_2x3_xy @ [x + flow_dxdy[y, x]; 1] sends an IMC pixel to the affine-aligned frame '
    '(Xenium mapped by M_dapi_um_to_imc_px). Images: warped(x) = img_affine(B(x)). Points: solve '
    'B(x) = q, see invert_points / map_xenium_um_to_imc in imc_xenium_align.py.')


def launch(imc_path, morph_dir, **kwargs):
    """Open the alignment viewer; call napari.run() afterwards when running as a script."""
    return AlignmentViewer(imc_path, morph_dir, **kwargs)
