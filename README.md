# Predict carbonate sediment thickness

Generate carbonate sediment thickness grids from age, mean distance and bathymetry grids over the time range 0-230Ma (in 1My increments).

## Prerequisites

You'll need to install GMT (and make sure the 'gmt' executable is in the PATH).
And also install the Python module 'scipy'.

The source code is compatible with Python 2.7.

## Usage

You can either run the Jupyter notebook `carbonate_sediment_thickness.ipynb` or run `python run_carbonate_sediment_thickness.py` in a console/terminal window.

In either case there are a bunch of top-level parameters that you can change/configure. You'll need to at least change the location of the age, bathymetry and distance grids.

The mean-distance-to-passive-margins grids can be downloaded from ftp://ftp.earthbyte.org/Data_Collections/Dutkiewicz_etal_2017_G3/passive_margin_mean_distance_grids/

The age grids can be downloaded from ftp://ftp.earthbyte.org/Data_Collections/Muller_Dutkiewicz_2018_SciAdv/Supplementary_grids/age_grids_AREPS2016/
