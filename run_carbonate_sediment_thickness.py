#
# Generate carbonate sediment thickness grids from age, mean distance and bathymetry grids
# over the time range 0-230Ma (in 1My increments).
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

# Times to generate sediment thickness grids.
# Must also have age and bathymetry grids at these times.
times = range(0, 231)

# Whether to use all CPU cores (parallel) or just one (serial).
# Note: Each process is set to a low priority so as not to interfere with your regular tasks.
use_all_cpu_cores = True

# CCD (calcite compensation depth) curve filename.
# This file maps time to CCD depth (negative).
ccd_curve_filename = 'input_data/Boss_Wilkinson_1991_global_CCD.txt'

# Maximum carbonate decompacted sediment rate curve filename (in cm/ky).
# This file maps time to the maximum carbonate rate (at mid-ocean ridge depth; reduces to zero at CCD).
max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename = 'input_data/sed_rate_400_to_0_Ma.txt'

# Location of age grids.
# The full path to grids including filename, but with time and filename extension removed.
age_grid_filename_prefix = 'E:/Users/John/Downloads/GPlates/data/rasters/Muller_etal_2016_AREPS_Agegrids/netCDF_0-230Ma/EarthByte_AREPS_Muller_etal_2016_AgeGrid-'
# Filename extension (typically 'nc' or 'grd').
age_grid_filename_extension = 'nc'

# Location of bathymetry grids.
# The full path to grids including filename, but with time and filename extension removed.
bathymetry_filename_prefix = 'E:/Users/John/Downloads/GPlates/data/rasters/Sioned/02_Paleobathymetry_with_seds_LIPs_Sioned/pbathy-with-seds-and-LIPs-AREPS-muller-etal-'
# Filename extension (typically 'nc' or 'grd').
bathymetry_filename_extension = 'grd'

# Location of mean-distance-to-passive-margins grids.
# The full path to grids including filename, but with time and filename extension removed.
distance_filename_prefix = 'E:/Users/John/Downloads/GPlates/data/PythonWorkflows/SedimentationRate/shortest_path_distance_grids_clamped_1d/grid_reg_mean_distance_1.0d_'
# Filename extension (typically 'nc' or 'grd').
distance_filename_extension = 'nc'

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
    carbonate_sediment_thickness.predict_sedimentation_and_write_data_for_times(
        times,
        grid_spacing,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        age_grid_filename_prefix,
        age_grid_filename_extension,
        distance_filename_prefix,
        distance_filename_extension,
        bathymetry_filename_prefix,
        bathymetry_filename_extension,
        carbonate_decompacted_sediment_thickness_filename_prefix,
        carbonate_compacted_sediment_thickness_filename_prefix,
        carbonate_deposition_mask_filename_prefix,
        use_all_cpu_cores)
