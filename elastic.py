"""Elastic refinement (ECC affine + Farneback optical flow), adapted from an
elastic_registration_opencv routine used in an H&E/IF training-tile pipeline.

Runs as a separate process (python elastic.py in.npz out.npz) because opencv-python
bundles its own Qt, which must not be loaded into the napari viewer process.

Convention (all in IMC pixels, x = col, y = row):
    moving(W · [x + flow(x); 1]) ≈ fixed(x)
i.e. ecc (2x3) and flow (H, W, 2 = dx, dy) form a backward map from the fixed (IMC) grid
into the moving image (DAPI already affine-warped into the IMC grid).
"""
import json
import sys

import numpy as np

DEFAULTS = dict(ecc=True, ecc_iters=100, pyr_scale=0.5, levels=3, winsize=15, iterations=3,
                poly_n=5, poly_sigma=1.2, smooth_sigma=3.0, strength=1.0)


def elastic_register(fixed, moving, params=None):
    """fixed, moving: uint8 (H, W), same grid. Returns ecc (2x3 float64), flow (H, W, 2 float32),
    ecc_ok (bool)."""
    import cv2

    p = {**DEFAULTS, **(params or {})}
    fixed = np.ascontiguousarray(fixed, np.uint8)
    moving = np.ascontiguousarray(moving, np.uint8)
    h, w = fixed.shape

    # 1. ECC affine refine (moving -> fixed); identity if it fails to converge
    ecc = np.eye(2, 3, dtype=np.float32)
    ecc_ok = False
    if p['ecc']:
        criteria = (cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT, int(p['ecc_iters']), 1e-3)
        try:
            _, ecc = cv2.findTransformECC(fixed.astype(np.float32), moving.astype(np.float32),
                                          ecc, cv2.MOTION_AFFINE, criteria, None, 1)
            ecc_ok = True
        except cv2.error:
            ecc = np.eye(2, 3, dtype=np.float32)
    moving_ecc = cv2.warpAffine(moving, ecc, (w, h), flags=cv2.INTER_LINEAR + cv2.WARP_INVERSE_MAP)

    # 2. dense optical flow on the ECC-warped moving image: fixed(x) ≈ moving_ecc(x + flow(x))
    flow = cv2.calcOpticalFlowFarneback(fixed, moving_ecc, None, float(p['pyr_scale']),
                                        int(p['levels']), int(p['winsize']), int(p['iterations']),
                                        int(p['poly_n']), float(p['poly_sigma']), 0)
    # 3. regularise + scale
    if p['smooth_sigma'] > 0:
        flow = cv2.GaussianBlur(flow, (0, 0), float(p['smooth_sigma']))
    flow = (flow * float(p['strength'])).astype(np.float32)
    return ecc.astype(np.float64), flow, ecc_ok


if __name__ == '__main__':
    src, dst = sys.argv[1], sys.argv[2]
    d = np.load(src)
    ecc, flow, ecc_ok = elastic_register(d['fixed'], d['moving'], json.loads(str(d['params'])))
    np.savez(dst, ecc=ecc, flow=flow, ecc_ok=ecc_ok)
