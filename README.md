# Predict carbonate sediment thickness

Generate carbonate sediment thickness grids from age, mean distance and bathymetry grids over the time range 0-230Ma (in 1My increments).

## Prerequisites

You'll need to install GMT (and make sure the 'gmt' executable is in the PATH).
And also install the Python module 'scipy'. And on Windows platforms also install the Python module 'psutil'.

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

The mean-distance-to-passive-margins grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Dutkiewicz_etal_2017_G3/passive_margin_mean_distance_grids/. Or you can run the workflow https://github.com/EarthByte/predicting-sediment-thickness to generate the grids.

The age grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF.zip .

The bathymetry grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Wright_etal_2020_ESR/Grids/Paleobathymetry_RHCW18/ .

## Reference

Dutkiewicz, A., Müller, R.D., Cannon, J., Vaughan, S. and Zahirovic, S., 2019, Sequestration and subduction of deep-sea carbonate in the global ocean since the Early Cretaceous. Geology, 47(1), pp.91-94. DOI:  https://doi.org/10.1130/G45424.1
