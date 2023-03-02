# Calculate carbonate sediment thickness

Generate carbonate sediment thickness grids from age and bathymetry grids over the time range 0-230Ma (in 1My increments).

## Dependencies

- [GMT](https://www.generic-mapping-tools.org/download/) (and make sure the 'gmt' executable is in the PATH).
- [PyGPlates](https://www.gplates.org/docs/pygplates/pygplates_getting_started.html#installation).
- SciPy.
- And on Windows platforms also install [psutil](https://pypi.org/project/psutil/).

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

The location of the age and bathymetry will need to be changed to point to your local grids.

The age grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_Agegrids/Muller_etal_2019_Tectonics_v2.0_netCDF.zip .

The bathymetry grids can be downloaded from https://www.earthbyte.org/webdav/ftp/Data_Collections/Wright_etal_2020_ESR/Grids/Paleobathymetry_RHCW18/ .

You can either use the supplied topological model (`2019_v2`) or provide your own. In either case you'll set the `topology_model_name` variable (in `carbonate_sediment_thickness.ipynb` or `run_carbonate_sediment_thickness.py`) to the name of the topological model. If you're using the supplied model you only need to specify `topology_model_name = '2019_v2'`. However if you're providing your own model then please follow these steps to use your model:

- Create a new sub-directory of `input_data/topology_model/` that is the name of your model (eg, `my_model`).
  - For example: `input_data/topology_model/my_model/`
- Copy your model files into that directory.
  - For example, copy the model GPML(Z) and ROT files into `input_data/topology_model/my_model/`.
- Create a new text file called `rotation_files.txt` that lists the model's *rotation* files.
  - For example: `input_data/topology_model/my_model/rotation_files.txt`
  - The paths of the rotation files (listed in `rotation_files.txt`) should be relative to the base directory. For example:
    ```
    input_data/topology_model/my_model/rotations_1.rot
    input_data/topology_model/my_model/rotations_2.rot
    ...
    ```
  - A convenient way to generate this file is (using the `my_model` example) is to run the following command from the base directory:
    ```
    ls -A1 input_data/topology_model/my_model/*.rot > input_data/topology_model/my_model/rotation_files.txt
    ```
- Create a new text file called `topology_files.txt` that lists the model's *topology* files.
  - For example: `input_data/topology_model/my_model/topology_files.txt`
  - The paths of the topology files (listed in `topology_files.txt`) should be relative to the base directory. For example:
    ```
    input_data/topology_model/my_model/topologies_1.gpmlz
    input_data/topology_model/my_model/topologies_2.gpmlz
    ...
    ```
  - A convenient way to generate this file is (using the `my_model` example) is to run the following command from the base directory. This assumes a mix of `.gpml` and `.gpmlz` files (also note the `>>` on the second command to append):
    ```
    ls -A1 input_data/topology_model/my_model/*.gpml > input_data/topology_model/my_model/topology_files.txt
    ls -A1 input_data/topology_model/my_model/*.gpmlz >> input_data/topology_model/my_model/topology_files.txt
    ```
- Then set the `topology_model_name` variable  to the name of your model (ie, name of the sub-directory).
  - For example: `topology_model_name = 'my_model'`
- Note: It may be helpful to look at the supplied `2019_v2` model as an example of how to do this.
  - That is, have a look in the `input_data/topology_model/2019_v2/` directory.


## Reference

Dutkiewicz, A., Müller, R.D., Cannon, J., Vaughan, S. and Zahirovic, S., 2019, Sequestration and subduction of deep-sea carbonate in the global ocean since the Early Cretaceous. Geology, 47(1), pp.91-94. DOI:  https://doi.org/10.1130/G45424.1
