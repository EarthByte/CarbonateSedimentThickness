#
# Generate carbonate sediment thickness grids from paleobathymetry grids.
#
# See "carbonate_sediment_thickness.py" for details of the algorithm.
#

import carbonate_sediment_thickness
import os.path
import os


#
# Input Parameters.
#

# Grid spacing of ocean basin points (in degrees).
grid_spacing = 0.5

# Regional latitude and longitude range.
# Keep latitudes in the range [-90, 90] and longitudes the in range [-180, 180].
# Set to [-90, 90] and [-180, 180] to cover the entire globe.
min_lat, max_lat = -90, 90
min_lon, max_lon = -180, 180

# Times to generate sediment thickness grids.
# Must also have age and bathymetry grids at these times.
times = list(range(0, 181))

# Whether to use all CPU cores (parallel) or just one (serial).
# Note: Each process is set to a low priority so as not to interfere with your regular tasks.
use_all_cpu_cores = True

# The topological model used to assign plate IDs to ocean crust at paleo times (including crust subducted at present day).
#
# This is the name of a sub-directory in 'input_data/topology_model/'.
# Currently the only builtin model is '2019_v2'.
# However you can provide your own topological model by following the instructions in the main README.
topology_model_name = '2019_v2'

# The anchor plate ID used to reconstruct grid points to sample paleobathymetry grids.
anchor_plate_id = 0

# CCD (calcite compensation depth) curve filename.
# This file maps time to CCD depth (negative).
ccd_curve_filename = 'input_data/Boss_Wilkinson_1991_global_CCD.txt'

# Maximum carbonate decompacted sediment rate curve filename (in cm/ky).
# This file maps time to the maximum carbonate rate (at mid-ocean ridge depth; reduces to zero at CCD).
max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename = 'input_data/sed_rate_v6.txt'

# Location of age grids.
# The full path to grids including filename, but with time and filename extension removed.
age_grid_filename_prefix = '/home/michael/workspace/CarbonateSedimentThickness/Muller_etal_2019_Tectonics_v2.0_netCDF/Muller_etal_2019_Tectonics_v2.0_AgeGrid-'
# The number of decimal places in the time when it is added to the filename prefix (and following by the filename extension).
# Typically either 0 or 1.
age_grid_filename_decimal_places_in_time = 0
# Filename extension (typically 'nc' or 'grd').
age_grid_filename_extension = 'nc'

# Location of bathymetry grids.
# The full path to grids including filename, but with time and filename extension removed.
bathymetry_filename_prefix = '/home/michael/workspace/CarbonateSedimentThickness/Paleobathymetry_RHCW18/paleobathymetry_'
# The number of decimal places in the time when it is added to the filename prefix (and following by the filename extension).
# Typically either 0 or 1.
bathymetry_filename_decimal_places_in_time = 0
# Filename extension (typically 'nc' or 'grd').
#bathymetry_filename_extension = 'grd'
bathymetry_filename_extension = 'nc'

# Time of the oldest bathymetry grid.
bathymetry_filename_oldest_time = 230

# Location of output carbonate decompacted and compacted sediment thickness grids.
# The full path to thickness grids including the base filename
# (grid spacing, time and filename extension will get added later).
output_data_dir = 'sediment_thickness'
carbonate_decompacted_sediment_thickness_filename_prefix = os.path.join(output_data_dir, 'decompacted_sediment_thickness')
carbonate_compacted_sediment_thickness_filename_prefix = os.path.join(output_data_dir, 'compacted_sediment_thickness')
carbonate_deposition_mask_filename_prefix = os.path.join(output_data_dir, 'deposition_mask')


if __name__ == '__main__':
    
    # Create output directory if it doesn't exist.
    if not os.path.exists(output_data_dir):
        os.makedirs(output_data_dir)
    
    # All functionality is delegated to "carbonate_sediment_thickness.py" to enable parallel processing.
    carbonate_sediment_thickness.calc_sedimentation_and_write_data_for_times(
        times,
        (min_lat, max_lat),
        (min_lon, max_lon),
        grid_spacing,
        topology_model_name,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        (age_grid_filename_prefix, age_grid_filename_decimal_places_in_time, age_grid_filename_extension),
        (bathymetry_filename_prefix, bathymetry_filename_decimal_places_in_time, bathymetry_filename_extension),
        bathymetry_filename_oldest_time,
        carbonate_decompacted_sediment_thickness_filename_prefix,
        carbonate_compacted_sediment_thickness_filename_prefix,
        carbonate_deposition_mask_filename_prefix,
        anchor_plate_id,
        use_all_cpu_cores)
