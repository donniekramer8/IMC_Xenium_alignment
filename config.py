"""Settings for run_alignment.py (Xenium DAPI -> IMC landmark alignment).

Edit the three input paths below, then run:  python run_alignment.py
On Windows write paths as r'C:\\data\\...' (note the r before the quote).
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

# inputs: CHANGE THESE
IMC_PATH = '/path/to/imc_roi.tiff'                                       # IMC ROI, all channels stacked
MORPH_DIR = '/path/to/xenium_output/morphology_focus'                    # contains ch0000_dapi.ome.tif
NUCLEI_PATH = '/path/to/xenium_output/nucleus_boundaries.parquet'        # Xenium nucleus outlines (µm); None to skip

# outputs: <OUT_DIR>/<IMC file name>_alignment.json, _alignment_qc.png, _elastic_field.npz
OUT_DIR = os.path.join(HERE, 'outputs')

# IMC nuclear image: arcsinh(sum(channels) / cofactor), percentile stretch -> uint8
IMC_NUC_CHANNELS = (37, 38)   # Ir191 / Ir193 DNA (0-based); adjust to your panel
IMC_COFACTOR = 20
IMC_PCT = (30, 99.5)

# Xenium DAPI: pyramid level (2 = ds x4, 0.85 µm/px), percentile stretch -> uint8
DAPI_LEVEL = 2
DAPI_N_LEVELS = 4             # extra coarser levels for smooth zooming
DAPI_PCT = (30, 99.5)

# Optional elastic step after Validate (ECC affine + Farneback optical flow, see elastic.py);
# these are the starting values of the controls in the viewer
ELASTIC = dict(
    ecc=True,           # ECC affine refinement before the flow
    ecc_iters=100,
    pyr_scale=0.5,      # Farneback parameters
    levels=3,
    winsize=15,
    iterations=3,
    poly_n=5,
    poly_sigma=1.2,
    smooth_sigma=3.0,   # Gaussian smoothing of the flow field (IMC px), 0 = off
    strength=1.0,       # scale applied to the flow (1.0 = full correction)
)
