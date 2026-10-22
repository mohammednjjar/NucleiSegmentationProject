# Nuclei Segmentation and Spatial Transcriptomics: Research Notebooks

Exploratory notebooks from the M.Sc. in Computer Engineering at Bahçeşehir
University, Istanbul (supervisor: Dr. Lavdie Rada). They record the work that
led to the final benchmark, which lives in its own repository:
**[moffitt-merfish-benchmark](https://github.com/mohammednjjar/moffitt-merfish-benchmark)**.

These notebooks are kept as a record. They are not maintained, and their file
paths point to the Kaggle or Colab environments they were run in.

## Contents

| Folder | What it holds |
|---|---|
| `UNET_SEGMENTATION/` | Four U-Net experiments for nuclei segmentation on the 2018 Data Science Bowl images |
| `VISIUM_HD/` | StarDist nuclei segmentation on Visium HD mouse small-intestine H&E images |
| `EROSION_DILATION/` | Morphological erosion and dilation with OpenCV on a sample image |
| `MANUAL_MORPHOLOGICAL_OPERATIONS/` | Thresholding, adaptive thresholding, blurring, histogram equalisation and dilation on Allen MERFISH, osmFISH and Visium HD nuclei images |
| `MOFFITT_PROJECT/` | Earlier stages of the Moffitt (2018) MERFISH work: data merging, boundary unification, Voronoi cell reconstruction, visualisation, the v4 Random Forest / LightGBM benchmark and the v8 Random Forest runs |
| `COMPARISON_METHODS/` | SpaGCN and DeepST run on a 50,000-cell subset of the Moffitt data |
| `MERFISH_ANALYSIS_PYTHON_TRANSLATION/` | Python translation of parts of the ZhuangLab MERFISH_analysis MATLAB code (see below) |

## Third-party code

`MERFISH_ANALYSIS_PYTHON_TRANSLATION/` translates parts of
[ZhuangLab/MERFISH_analysis](https://github.com/ZhuangLab/MERFISH_analysis)
from MATLAB to Python, and includes the original MATLAB files for reference.
The original code is © Presidents and Fellows of Harvard College and remains
under the terms of its original repository.
