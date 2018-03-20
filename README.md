# Predict carbonate sediment thickness

Generate carbonate sediment thickness grids from age, mean distance and bathymetry grids over the time range 0-230Ma (in 1My increments).

## Prerequisites

You'll need to install GMT (and make sure the 'gmt' executable is in the PATH).
And also install the Python module 'scipy'.

The source code is compatible with Python 2.7.

## Usage

To generate the carbonate thickness grids you can either:

- load the Jupyter notebook `carbonate_sediment_thickness.ipynb` and run all cells, or
- type `python run_carbonate_sediment_thickness.py` in a console/terminal window.

In either case there are a bunch of top-level parameters that you can change/configure.
By default `use_all_cpu_cores` is set to `True` to run on all CPU cores
(otherwise it takes too long; up to 25 hours at 0.5 degree resolution using just a single core).
Note that you can increase the `grid_spacing` parameter to reduce the running time.

**Note:** If you choose the Jupyter notebook *and* you edit a parameter *outside* the notebook
(such as inside the imported module *carbonate_sediment_thickness*) then you'll need to restart the notebook kernel
after each modification (or insert `reload(carbonate_sediment_thickness)` after `import carbonate_sediment_thickness`).

The location of the age, bathymetry and distance grids will need to be changed to point to your local grids.

The mean-distance-to-passive-margins grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Dutkiewicz_etal_2017_G3/passive_margin_mean_distance_grids/

The age grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_Dutkiewicz_2018_SciAdv/Supplementary_grids/age_grids_AREPS2016/
