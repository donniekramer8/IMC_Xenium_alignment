"""Open the napari landmark viewer using the settings in config.py.

    conda activate imc_xenium_align
    python run_alignment.py              # or: python run_alignment.py my_other_config.py

Place >=6 landmark pairs, click Validate (Xenium nuclei on IMC), optionally Run elastic
(re-run with other parameters until it looks right), then Save -> JSON with the DAPI -> IMC
affine (+ elastic field .npz if used) + QC PNG (overwrites previous results).
"""
import os
import sys

# opencv-python bundles its own Qt and exports this, which breaks PyQt5 -> drop it before Qt loads
os.environ.pop('QT_QPA_PLATFORM_PLUGIN_PATH', None)
os.environ.setdefault('QT_API', 'pyqt5')

import importlib.util  # noqa: E402

import napari  # noqa: E402

from imc_xenium_align import launch  # noqa: E402


def load_config(path):
    spec = importlib.util.spec_from_file_location('config', path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)
    return cfg


def check_inputs(cfg, cfg_path):
    """Stop with a readable message if config.py still has placeholder / wrong paths."""
    checks = [('IMC_PATH', cfg.IMC_PATH, os.path.isfile),
              ('MORPH_DIR', cfg.MORPH_DIR,
               lambda p: os.path.isfile(os.path.join(p, 'ch0000_dapi.ome.tif')))]
    if getattr(cfg, 'NUCLEI_PATH', None):
        checks.append(('NUCLEI_PATH', cfg.NUCLEI_PATH, os.path.isfile))
    bad = [(name, path) for name, path, ok in checks if not ok(path)]
    if bad:
        print(f'\nPlease edit the input paths in {cfg_path}:')
        for name, path in bad:
            extra = ' (folder must contain ch0000_dapi.ome.tif)' if name == 'MORPH_DIR' else ''
            print(f'  {name} = {path!r}  -> not found{extra}')
        sys.exit(1)


if __name__ == '__main__':
    here = os.path.dirname(os.path.abspath(__file__))
    cfg_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'config.py')
    cfg = load_config(cfg_path)
    check_inputs(cfg, cfg_path)

    os.makedirs(cfg.OUT_DIR, exist_ok=True)
    roi = os.path.splitext(os.path.basename(cfg.IMC_PATH))[0]
    al = launch(
        cfg.IMC_PATH, cfg.MORPH_DIR,
        out_path=os.path.join(cfg.OUT_DIR, f'{roi}_alignment.json'),
        dapi_level=cfg.DAPI_LEVEL,
        nuclei_path=getattr(cfg, 'NUCLEI_PATH', None),
        imc_kwargs=dict(chans=cfg.IMC_NUC_CHANNELS, cofactor=cfg.IMC_COFACTOR, pct=cfg.IMC_PCT),
        dapi_kwargs=dict(n_levels=cfg.DAPI_N_LEVELS, pct=cfg.DAPI_PCT),
        elastic_defaults=getattr(cfg, 'ELASTIC', None),
    )
    napari.run()
