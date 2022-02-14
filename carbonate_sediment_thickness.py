
"""
    Copyright (C) 2018 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""


############################################################################################
# Generate carbonate sediment thickness grids from age, mean distance and bathymetry grids #
# over the time range 0-230Ma (in 1 Myr increments).                                       #
############################################################################################



from call_system_command import call_system_command
import csv
import math
import multiprocessing
import os.path
import pybacktrack
from scipy.interpolate import interp1d
import sys


#
# Input Parameters.
#


# Average ocean floor sediment density (kg/m^3), porosity (%) and decay (metres).
AVERAGE_OCEAN_FLOOR_SEDIMENT_DENSITY = 2647.0
AVERAGE_OCEAN_FLOOR_SEDIMENT_POROSITY = 0.66
AVERAGE_OCEAN_FLOOR_SEDIMENT_DECAY = 1333.0

# Density mantle and water (kg/m^3).
DENSITY_MANTLE = 3330.0
DENSITY_WATER = 1030.0

# Dissolution distance (in metres) above the CCD.
# This distance above the CCD is where carbonate dissolution is 100%.
# Below this distance dissolution decreases linearly to 0% at the CCD.
CCD_DISSOLUTION_DISTANCE = 300.0

# Polynomial coefficients for predicting compacted sediment thickness and decompacted sediment rate from
# age and mean distance to passive margins in the paper...
#
#   "Predicting sediment thickness on vanished ocean crust since 200 Ma".
#
SEDIMENT_THICKNESS_POLYNOMIAL_COEFFICIENTS = [
    5.3732890044120172, 0.44092176, -0.15401756,
    -0.23843168, -0.06208386, 0.00957594,
    0.07160925, 0.00344379, 0.0, -0.32534525]
SEDIMENT_RATE_POLYNOMIAL_COEFFICIENTS = [
    -1.0051420927669603, -0.30916322, -0.19923406,
    0.38827883, -0.12533169, 0.0,
    -0.11374372, 0.0297582, -0.02391933, -0.36943835]
# Parameters used to standardize age and distance for both the sediment thickness and rate polynomials.
SEDIMENT_POLYNOMIAL_MEAN_AGE = 60.1842831
SEDIMENT_POLYNOMIAL_MEAN_DISTANCE = 1878.23124959
SEDIMENT_POLYNOMIAL_VARIANCE_AGE = 1893.88287649
SEDIMENT_POLYNOMIAL_STD_DEVIATION_AGE = math.sqrt(SEDIMENT_POLYNOMIAL_VARIANCE_AGE)
SEDIMENT_POLYNOMIAL_VARIANCE_DISTANCE = 1159561.12717194
SEDIMENT_POLYNOMIAL_STD_DEVIATION_DISTANCE = math.sqrt(SEDIMENT_POLYNOMIAL_VARIANCE_DISTANCE)
SEDIMENT_POLYNOMIAL_MAX_AGE = 196.88598633
SEDIMENT_POLYNOMIAL_MAX_DISTANCE = 3000.


#
# Functions to read/generate input data.
#


# Generate uniformly spaced lat/lon grid of points in the specified lat/lon range.
# Points use GMT grid registration (start and end exactly on the dateline when using global lat/lon range).
def generate_input_points_grid(
        grid_spacing_degrees,
        latitude_range,   # (min, max) tuple
        longitude_range): # (min, max) tuple
    
    min_lat, max_lat = latitude_range
    min_lon, max_lon = longitude_range
    
    if (grid_spacing_degrees == 0 or
        max_lat <= min_lat or
        max_lon <= min_lon):
        return
    
    input_points = []
    
    # Data points start *on* 'min_lat' and 'min_lon'.
    # If '(max_lat - min_lat)' is an integer multiple of grid spacing then final latitude also lands on 'max_lat'.
    # If '(max_lon - min_lon)' is an integer multiple of grid spacing then final longitude also lands on 'max_lon'.
    num_latitudes = int(math.floor((max_lat - min_lat) / grid_spacing_degrees)) + 1
    num_longitudes = int(math.floor((max_lon - min_lon) / grid_spacing_degrees)) + 1
    for lat_index in range(num_latitudes):
        lat = min_lat + lat_index * grid_spacing_degrees
        
        for lon_index in range(num_longitudes):
            lon = min_lon + lon_index * grid_spacing_degrees
            
            input_points.append((lon, lat))
    
    return input_points


# Returns a list of scalars (one per (lon, lat) point in the 'input_points' list).
# For input points outside the scalar grid then scalars will be Nan (ie, 'math.isnan(scalar)' will return True).
def get_positions_and_scalars(input_points, scalar_grid_filename, max_scalar=None):
    
    input_points_data = ''.join('{0} {1}\n'.format(lon, lat) for lon, lat in input_points)

    # The command-line strings to execute GMT 'grdtrack'.
    grdtrack_command_line = ["gmt", "grdtrack", "-nl", "-G{0}".format(scalar_grid_filename)]
    stdout_data = call_system_command(grdtrack_command_line, stdin=input_points_data.encode('utf-8'), return_stdout=True)
    
    lon_lat_scalar_list = []
    
    # Read lon, lat and scalar values from the output of 'grdtrack'.
    for line in stdout_data.splitlines():
        if (line.strip().startswith(b'#') or
            line.strip().startswith(b'>')):
            continue
        
        line_data = line.split()
        num_values = len(line_data)
        
        # If just a line containing white-space then skip to next line.
        if num_values == 0:
            continue
        
        if num_values < 3:
            print('Ignoring line "{0}" - has fewer than 3 white-space separated numbers.'.format(line), file=sys.stderr)
            continue
            
        try:
            # Convert strings to numbers.
            lon = float(line_data[0])
            lat = float(line_data[1])
            
            # The scalar got appended to the last column by 'grdtrack'.
            scalar = float(line_data[-1])
            
            # If the point is outside the grid then the scalar grid will return 'NaN'.
            if math.isnan(scalar):
                # print('Ignoring line "{0}" - point is outside scalar grid.'.format(line), file=sys.stderr)
                continue
            
            # Clamp to max value if requested.
            if (max_scalar is not None and
                scalar > max_scalar):
                scalar = max_scalar
            
        except ValueError:
            print('Ignoring line "{0}" - cannot read floating-point lon, lat and scalar values.'.format(line), file=sys.stderr)
            continue
        
        lon_lat_scalar_list.append((lon, lat, scalar))
    
    return lon_lat_scalar_list


def read_curve(curve_filename):
    
    x = []
    y = []
    with open(curve_filename, 'r') as curve_file:
        curve_reader = csv.reader(curve_file, delimiter='\t',)
        for row in curve_reader:
            # Skip comments.
            if (row[0].startswith("#") or
                row[0].startswith(">")):
                continue

            if len(row) < 2:
                raise ValueError('Curve file "{0}" does not have at least 2 columns at line {1}.'.format(
                                 curve_filename,
                                 curve_reader.line_num - 1))

            x.append(float(row[0]))
            y.append(float(row[1]))
    
    # Convert curve to a piecewise linear interpolating 1D function that accepts x and returns y.
    # Return NaN on bounds error instead of raising ValueError...
    return interp1d(x, y, bounds_error=False, fill_value='extrapolate')


#
# Functions to write output data.
#


def write_xyz_file(output_filename, output_data):
    with open(output_filename, 'w') as output_file:
        for output_line in output_data:
            output_file.write(' '.join(str(item) for item in output_line) + '\n')


def write_grid_file_from_xyz(
        grid_filename,
        xyz_filename,
        grid_spacing,
        latitude_range,
        longitude_range):
    
    # The command-line strings to execute GMT 'nearneighbor'.
    # For example "nearneighbor thickness.xy -R-180/180/-90/90 -I1 -N4 -S1d -Gthickness.nc".
    gmt_command_line = [
        "gmt",
        "nearneighbor",
        xyz_filename,
        "-N4",
        "-S{0}d".format(1.5 * grid_spacing),
        "-I{0}".format(grid_spacing),
        "-R{0}/{1}/{2}/{3}".format(longitude_range[0], longitude_range[1], latitude_range[0], latitude_range[1]),
        # Use GMT gridline registration since our input point grid has data points on the grid lines.
        # Gridline registration is the default so we don't need to force pixel registration...
        # "-r",  # Force pixel registration since data points are at centre of cells.
        "-G{0}".format(grid_filename)]
    call_system_command(gmt_command_line)


def write_data(
        data,
        output_filename_prefix,
        grid_spacing,
        latitude_range,
        longitude_range):
    
    # Write XYZ file.
    xyz_filename = '{0}.xy'.format(output_filename_prefix)
    write_xyz_file(xyz_filename, data)
    
    # Write grid file.
    grid_filename = '{0}.nc'.format(output_filename_prefix)
    write_grid_file_from_xyz(grid_filename, xyz_filename, grid_spacing, latitude_range, longitude_range)


#########################################################################################
# GDH1 model (Stein and Stein 1992)                                                     #
# "Model for the global variation in oceanic depth and heat flow with lithospheric age" #
#########################################################################################


def age_to_depth_GDH1(age):
    """Convert age to depth (negative)."""
    
    if age < 0:
        raise ValueError('Age must be non-negative')
    elif age < 20:
        return -(2600.0 + 365.0 * math.sqrt(age))
    else:
        return -(5651.0 - 2473.0 * math.exp(-0.0278 * age))


def depth_to_age_GDH1(depth):
    """
    Convert depth (negative) to age.
    
    This is the reverse of GDH1 age to depth:
        -depth = 2600 + 365 * age                       ; age <  20
        -depth = 5651 - 2473 * e^(-0.0278 * age)        ; age >= 20
    ...which is:
        age = 0                                         ; -depth <  2600
        age = ((-depth - 2600) / 365)^2                 ; -depth <  4232.5
        age = ln((5651.0 + depth) / 2473.0) / -0.0278   ; -depth >= 4232.5
    """
    
    if depth > 0:
        raise ValueError('Depth must be negative')
    elif depth > -2600.0:
        # Minimum depth is 2600.
        return 0.0
    elif depth > -4232.5:  # depth at age 20
        return math.pow((-depth - 2600.0) / 365.0, 2.0)
    else:
        return math.log((5651.0 + depth) / 2473.0) / -0.0278


#############################################################################
# Richards et al. (2020)                                                    #
# "Structure and dynamics of the oceanic lithosphere-asthenosphere system". #
#############################################################################

_RHCW18_age_to_depth_function = None

def age_to_depth_RHCW18(age):

    # Create the model function the first time we're called.
    global _RHCW18_age_to_depth_function
    if _RHCW18_age_to_depth_function is None:
        RHCW18_age_to_depth_file = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'input_data', 'RHCW18', 'depth-1333-130-2500.dat')
        # Read the age-to-depth curve depth=function(age) from age-to-depth data file.
        _RHCW18_age_to_depth_function = read_curve(RHCW18_age_to_depth_file)

    return _RHCW18_age_to_depth_function(age)


#
# Functions to predict, compact and isostatically correct sediment thickness.
#


def predict_total_compacted_sediment_thickness(
        age,
        distance):
    """
    Predict total compacted sediment thickness from age and mean distance to passive margins.
    """
    
    # We need to remove the mean and scale to unit variance (for the age and distance values)
    # based on the machine learning training scaler.
    # See http://scikit-learn.org/stable/modules/preprocessing.html#standardization-or-mean-removal-and-variance-scaling
    age = (age - SEDIMENT_POLYNOMIAL_MEAN_AGE) / SEDIMENT_POLYNOMIAL_STD_DEVIATION_AGE
    distance = (distance - SEDIMENT_POLYNOMIAL_MEAN_DISTANCE) / SEDIMENT_POLYNOMIAL_STD_DEVIATION_DISTANCE
    
    polynomial_features = [
        1, age, distance,
        age * age, age * distance, distance * distance,
        age * age * age, age * age * distance, age * distance * distance, distance * distance * distance]
    
    # Evaluate the polynomial to get the log of the predicated sediment thickness.
    log_sediment_thickness = sum(
        SEDIMENT_THICKNESS_POLYNOMIAL_COEFFICIENTS[i] * polynomial_features[i]
        for i in range(10))
    
    # Return the predicted sediment thickness (not as a logarithm).
    return math.exp(log_sediment_thickness)


def predict_total_decompacted_sediment_rate(
        age,
        distance):
    """
    Predict total decompacted sediment rate from age and mean distance to passive margins.
    """
    
    # We need to remove the mean and scale to unit variance (for the age and distance values)
    # based on the machine learning training scaler.
    # See http://scikit-learn.org/stable/modules/preprocessing.html#standardization-or-mean-removal-and-variance-scaling
    age = (age - SEDIMENT_POLYNOMIAL_MEAN_AGE) / SEDIMENT_POLYNOMIAL_STD_DEVIATION_AGE
    distance = (distance - SEDIMENT_POLYNOMIAL_MEAN_DISTANCE) / SEDIMENT_POLYNOMIAL_STD_DEVIATION_DISTANCE
    
    polynomial_features = [
        1, age, distance, age * age,
        age * distance, distance * distance, age * age * age,
        age * age * distance, age * distance * distance, distance * distance * distance]
    
    # Evaluate the polynomial to get the log of the predicated sedimentation rate.
    log_sedimentation_rate = sum(
        SEDIMENT_RATE_POLYNOMIAL_COEFFICIENTS[i] * polynomial_features[i]
        for i in range(10))
    
    # Return the predicted sedimentation rate (not as a logarithm).
    #
    # NOTE: The factor of 10 converts cm/Ky to m/My (the polynomial generates units of cm/Ky).
    return 10.0 * math.exp(log_sedimentation_rate)


# Decompact total sediment thickness assuming a single lithology of the specified surface porosity and exponential decay.
#
# The decompacted (ie, surface) thickness d(z) of a compacted layer of height t(z) at depth z is:
#
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# ...due to the fact that the grain volume doesn't change with depth:
#
#   (1 - porosity(z1)) * t(z1) = (1 - porosity(z2)) * t(z2)
#
# ...and substituting z1=z and z2=0, and using d(z) = t(0) gives:
#
#   (1 - porosity(z)) * t(z) = (1 - porosity(0)) * d(z)
#   d(z) = t(z) * (1 - porosity(z)) / (1 - porosity(0))
#
# Where the porosity is:
#
#   porosity(z) = porosity(0) * exp(-c * z)
#
# The total decompacted thickness is:
#
#   D = Integral(d(z), z = 0 -> T)
#
# ...where T is the total compacted thickness.
#
#   D = Integral(dz * (1 - porosity(z)) / (1 - porosity(0)), z = 0 -> T)
#   D = Integral(dz * (1 - porosity(0) * exp(-c * z)) / (1 - porosity(0)), z = 0 -> T)
#
#     = (T + (porosity(0) / c) * (exp(-c * T) - 1)) / (1 - porosity(0))
#
# Rearranging gives...
#
#   T = D * (1 - porosity(0)) - (porosity(0) / c) * (exp(-c * T) - 1)
#     = a * exp(-c * T) + b
#
# ...where a and b are:
#
#    a = -porosity(0) / c
#
#    b = D * (1 - porosity(0)) - a
#
# So the compacted thickness T:
#
#    T = a * exp(-c * T) + b
#
# ...can be solved iteratively by repeatedly substituting the left hand side back into the right hand side
# until T converges on a solution. The initial T is chosen to be D (the decompacted thickness).
#
def compact_sediment_thickness(thickness, surface_porosity, porosity_exp_decay):
    
    # Constants 'a' and 'b' are calculated outside the iteration loop for efficiency.
    a = -surface_porosity / porosity_exp_decay
    b = thickness * (1 - surface_porosity) - a

    # Start out with initial estimate - choose the decompacted thickness.
    compacted_thickness = thickness

    # Limit the number of iterations in case we never converge.
    # Although should converge within around 20 iterations (for 1e-6 accuracy).
    for iteration in range(1000):
        new_compacted_thickness = a * math.exp(-porosity_exp_decay * compacted_thickness) + b
        
        # If we've converged within a tolerance then we're done.
        if math.fabs(new_compacted_thickness - compacted_thickness) < 1e-6:
            compacted_thickness = new_compacted_thickness
            break
        
        # The new thickness becomes the old thickness for the next loop iteration.
        compacted_thickness = new_compacted_thickness
    
    return compacted_thickness


def sediment_isostatic_correction(sediment_thickness, average_sediment_density):
    """
    Returns the isostatic correction of sediment.

    The returned correction can be added to a known water depth to obtain the deeper isostatically compensated,
    sediment-free water depth (tectonic subsidence). Or the correction could be subtracted from a
    known tectonic subsidence (unloaded water depth) to get the depth at sediment/water interface.
    """

    return sediment_thickness * (DENSITY_MANTLE - average_sediment_density) / (DENSITY_MANTLE - DENSITY_WATER)


def bathymetry_from_tectonic_subsidence_and_sedimentation(
        age,
        distance,
        age_to_depth_curve,
        dynamic_topography=0.0,
        bathymetry_model_adjustment=0.0):
    """
    Calculates tectonic subsidence of ocean floor (from age/depth curve) and predicts total compacted
    sediment thickness. Both are then combined to get bathymetry (with an offset correction).
    """
    
    # Tectonic subsidence is just the age-to-depth model of ocean basement (ie, sediment free), plus optional dynamic topography.
    tectonic_subsidence = age_to_depth_curve(age) + dynamic_topography
    
    # Predict total compacted sediment thickness for the current age and mean distance to passive margins.
    total_compacted_sediment_thickness = predict_total_compacted_sediment_thickness(age, distance)
    
    # Load sediment on top of the sediment-free depth.
    # Note we add (instead of subtract) because bathymetry and subsidence are negative here (instead of positive).
    bathymetry = tectonic_subsidence + sediment_isostatic_correction(
        total_compacted_sediment_thickness, AVERAGE_OCEAN_FLOOR_SEDIMENT_DENSITY)
    
    # print('ts, ct, it, b: ',
    #         tectonic_subsidence,
    #         total_compacted_sediment_thickness,
    #         bathymetry - tectonic_subsidence,
    #         bathymetry)
    
    # Add in the constant bathymetry offset (difference between actual bathymetry and bathymetry model).
    return bathymetry + bathymetry_model_adjustment


def predict_carbonate_decompacted_sediment_thickness(
        lon_lat_age_distance_bathymetry_list,
        age_to_depth_curve,
        dynamic_topography_model_name,  # Can be None
        ccd_curve,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        time):
    """
    Predicts carbonate decompacted sediment thickness.
    
    Returns a list of accumulated carbonate sediment thicknesses (one for each input point).
    """
    
    # The number of ocean points at 'time'.
    num_points = len(lon_lat_age_distance_bathymetry_list)

    # Extract each input data scalar into its own list.
    locations = []  # (lon, lat) tuples
    ages = []
    distances = []
    bathymetrys = []
    for lon, lat, age, distance, bathymetry in lon_lat_age_distance_bathymetry_list:
        locations.append((lon, lat))
        ages.append(age)
        distances.append(distance)
        bathymetrys.append(bathymetry)

    # Initialise all carbonate deposition to zero.
    # We will accumulate carbonate deposition only above the CCD.
    carbonate_decompacted_sediment_thicknesses = [0.0] * num_points

    # Create the dynamic topography model (from model name) if requested.
    if dynamic_topography_model_name:
        dynamic_topography_model = pybacktrack.InterpolateDynamicTopography.create_from_bundled_model(dynamic_topography_model_name)
    else:
        dynamic_topography_model = None
    
    # Initialise dynamic topography for all points at 'time'.
    if dynamic_topography_model:
        # Sample dynamic topography for all points in one call (for efficiency).
        dynamic_topography_at_time = dynamic_topography_model.sample(time, locations)
    else:
        # Use zero values for dynamic topography when no dynamic topography model is present.
        dynamic_topography_at_time = [0.0] * num_points

    # Initialise the offsets between actual bathymetry and modelled bathymetry for all ocean points at 'time'.
    bathymetry_model_adjustments = []
    for point_index in range(num_points):

        # Model bathymetry from tectonic subsidence and predicted total compacted sediment thickness.
        bathymetry_from_model = bathymetry_from_tectonic_subsidence_and_sedimentation(
                ages[point_index], distances[point_index], age_to_depth_curve, dynamic_topography_at_time[point_index])
        
        # There will be a difference between the actual bathymetry and the bathymetry modelled on
        # age-to-depth subsidence. Assume this offset is constant over the lifetime of this parcel of
        # ocean crust and use it to adjust the bathymetry obtained from age-to-depth model for younger ages.
        bathymetry_model_adjustment = bathymetrys[point_index] - bathymetry_from_model

        bathymetry_model_adjustments.append(bathymetry_model_adjustment)
    
    # Time interval to reconstruct ocean points (from 'time' to each ocean point's time of appearance).
    # Currently set to 1 Myr.
    time_interval = 1.0
    
    # The '0.5 * time_interval' is so we sample at interval mid-points for better integration accuracy.
    reconstruction_time = time + 0.5 * time_interval
    reconstructed_ages = [(age - 0.5 * time_interval) for age in ages]
    reconstructed_point_indices = list(range(num_points))
    while True:

        # Remove any ocean points that don't exist at the current reconstruction time
        # (ie, they appeared after current reconstruction time).
        index = 0
        while index < len(reconstructed_point_indices):
            reconstructed_point_index = reconstructed_point_indices[index]
            if reconstructed_ages[reconstructed_point_index] <= 0:
                del reconstructed_point_indices[index]
                continue
            index += 1
        
        # Finished if all ocean points appeared after current reconstruction time.
        if not reconstructed_point_indices:
            break
        
        # Get dynamic topography values for ocean points existing at the current reconstruction time.
        if dynamic_topography_model:
            reconstructed_locations = [locations[reconstructed_point_index]
                    for reconstructed_point_index in reconstructed_point_indices]
            
            # Sample dynamic topography for all reconstructed points in one call (for efficiency).
            dynamic_topography_at_reconstruction_time = dynamic_topography_model.sample(
                    reconstruction_time, reconstructed_locations)
        
        else:
            # Use zero values for dynamic topography when no dynamic topography model is present.
            dynamic_topography_at_reconstruction_time = [0.0] * len(reconstructed_point_indices)
            
        # The CCD depth associated with the reconstructed time.
        ccd_depth_at_reconstruction_time = ccd_curve(reconstruction_time)

        # Reconstruct bathymetry for all ocean points that exist at the current reconstruction time.
        for index, reconstructed_point_index in enumerate(reconstructed_point_indices):
            
            # Calculate bathymetry from tectonic subsidence and predicted total compacted sediment thickness.
            #
            # This is modelled from tectonic subsidence and predicted total compacted sediment thickness.
            #
            # Note: We predict total compacted sediment thickness for the younger age but assume the mean distance
            # to passive margins remains the same for younger ages (ie, is constant over the lifetime).
            # This means we can just sample the pre-generated mean-distance grids and not have to worry about
            # redoing those very lengthy calculations here.
            bathymetry_at_reconstruction_time = bathymetry_from_tectonic_subsidence_and_sedimentation(
                    reconstructed_ages[reconstructed_point_index],
                    distances[reconstructed_point_index],
                    age_to_depth_curve,
                    dynamic_topography_at_reconstruction_time[index],  # Note the use of 'index' instead of 'reconstructed_point_index'
                    bathymetry_model_adjustments[reconstructed_point_index])
            
            # If modelled bathymetry is above the CCD depth then add a 1 Myr time interval to the time spent above CDD.
            if bathymetry_at_reconstruction_time > ccd_depth_at_reconstruction_time:
                
                # Linearly interpolate the carbonate sedimentation rate with depth between the maximum rate
                # at CCD_DISSOLUTION_DISTANCE metres above the CDD and zero (at the CCD).
                #
                # Because the dissolution of carbonate in the ocean is pressure and temperature dependent,
                # more and more carbonate gets dissolved as you go down. This effect ultimately causes the
                # CCD to occur, but is actually reduces the sedimentation rates as we go down from the MOR
                # towards deeper seafloor, and this phenomenon can be clearly seen in the sedimentation rates
                # in the graph computed in "Predicting sediment thickness on vanished ocean crust since 200 Ma".
                # You can see (in the previously mentioned sed rate graph) that the sed rate declines substantially
                # from the MOR depth at 0 Ma to the depth of the recent CCD (around 45-50 Ma) where the sed rate
                # flattens because below the CCD sed rate is no longer depth/age dependent. From published papers
                # it appears that it is a good assumption that the sed rate decreases linearly between the MOR depth
                # and the CCD. We know that it is zero for carbonates at the CCD, and the (maximum) value for the
                # MOR is obtained from a time-dependent carbonate sed rate file.
                # Of course if we have an oceanic plateau that is substantially shallower than the normal MOR this
                # rate might be even higher than what the graph implies for the MOR but we cannot do anything about that.
                # So what we want to do is make the sed rate linearly dependent on the depth range from the MOR to the CCD.
                #
                # UPDATE: It was found that starting at the MOR ridge (2.5km) for 100% dissolution was too conservative.
                #         Instead the dissolution zone is now always a fixed distance 'CCD_DISSOLUTION_DISTANCE' above the CCD.
                #
                # Note that previously we used the sed rate far away from cont margins, ie at 3km distance,
                # as a proxy for carbonate sed rate. This used the (maximum) value for the MOR from the
                # polynomial relationship for sed rate, ie we just plugged in 0 Ma and 3000 km. This assumed
                # the carbonate sedimentation rate equals the total sedimentation rate because the maximum mean
                # distance to margins (mentioned above) essentially removes sediment contributions from continents.
                # However now we have a fixed carbonate sed rate curve loaded from a file instead.
                #
                # Note that this is a decompacted thickness rate (ie, deposited at surface where no compaction takes place,
                # but porosity effects still apply, albeit only the *surface* porosity).
                #
                # NOTE: The factor of 10 converts cm/Ky to m/My (the file has units of cm/Ky).
                max_carbonate_decompacted_sediment_rate = 10 * max_carbonate_decomp_sed_rate_cm_per_ky_curve(reconstruction_time)
                carbonate_decompacted_sediment_rate = (
                    max_carbonate_decompacted_sediment_rate *
                    min(bathymetry_at_reconstruction_time - ccd_depth_at_reconstruction_time, CCD_DISSOLUTION_DISTANCE) /
                    CCD_DISSOLUTION_DISTANCE
                )
                carbonate_decompacted_sediment_thicknesses[reconstructed_point_index] += time_interval * carbonate_decompacted_sediment_rate

            # Decrease the age of crust as we step backward in time by one time interval.
            reconstructed_ages[reconstructed_point_index] -= time_interval

        # Increase the reconstruction time as we step backward in time by one time interval.
        reconstruction_time += time_interval
    
    return carbonate_decompacted_sediment_thicknesses


def predict_sedimentation(
        input_points,  # List of (lon, lat) tuples,
        age_grid_filename,
        distance_filename,
        bathymetry_filename,
        age_to_depth_curve,
        dynamic_topography_model_name,  # Can be None
        ccd_curve,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        time):
    """
    Predicts decompacted and compacted carbonate thickness for each ocean basin grid point.
    Also determines when carbonate is currently being deposited (ie, only at the current time).
    
    Returns: A 3-tuple of lists each containing 3-tuples of (lon, lat, value).
             The values in the first list contain decompacted thickness.
             The values in the second list contain compacted thickness.
             The values in the third list contain a mask (value 1.0) for regions where
             carbonate is currently being deposited (ie, only at the current time).
    """
    
    # Get the input point ages, distances and bathymetries.
    lon_lat_age_list = get_positions_and_scalars(input_points, age_grid_filename, SEDIMENT_POLYNOMIAL_MAX_AGE)
    lon_lat_distance_list = get_positions_and_scalars(input_points, distance_filename, SEDIMENT_POLYNOMIAL_MAX_DISTANCE)
    lon_lat_bathymetry_list = get_positions_and_scalars(input_points, bathymetry_filename)
    if not lon_lat_age_list or not lon_lat_distance_list or not lon_lat_bathymetry_list:
        return
    
    # Merge the age, distance and bathymetry lists.
    # Only keep points where there are age *and* distance *and* bathymetry values.
    lon_lat_age_distance_bathymetry_list = []
    age_dict = dict(((lon, lat), age) for lon, lat, age in lon_lat_age_list)
    distance_dict = dict(((lon, lat), distance) for lon, lat, distance in lon_lat_distance_list)
    for lon, lat, bathymetry in lon_lat_bathymetry_list:
        if (lon, lat) in age_dict and (lon, lat) in distance_dict:
            age = age_dict[(lon, lat)]
            distance = distance_dict[(lon, lat)]
            lon_lat_age_distance_bathymetry_list.append((lon, lat, age, distance, bathymetry))
    
    # Predict *decompacted* thickness.
    # We do this for all ocean basin points together since it's more efficient to sample dynamic topography for all points at once.
    carbonate_decompacted_sediment_thickness_list = predict_carbonate_decompacted_sediment_thickness(
        lon_lat_age_distance_bathymetry_list,
        age_to_depth_curve, dynamic_topography_model_name, ccd_curve, max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        time)
    
    # The CCD at the current time.
    ccd_at_current_time = ccd_curve(time)
    
    # For each ocean basin point also generate *compacted* thickness.
    lon_lat_carbonate_decompacted_sediment_thickness_list = []
    lon_lat_carbonate_compacted_sediment_thickness_list = []
    lon_lat_carbonate_deposition_mask_list = []
    for grid_sample_index, (lon, lat, age, distance, bathymetry) in enumerate(lon_lat_age_distance_bathymetry_list):
        
        # Compact thickness.
        carbonate_decompacted_sediment_thickness = carbonate_decompacted_sediment_thickness_list[grid_sample_index]
        carbonate_compacted_sediment_thickness = compact_sediment_thickness(
            carbonate_decompacted_sediment_thickness,
            AVERAGE_OCEAN_FLOOR_SEDIMENT_POROSITY,
            AVERAGE_OCEAN_FLOOR_SEDIMENT_DECAY)
        
        lon_lat_carbonate_decompacted_sediment_thickness_list.append((lon, lat, carbonate_decompacted_sediment_thickness))
        lon_lat_carbonate_compacted_sediment_thickness_list.append((lon, lat, carbonate_compacted_sediment_thickness))
        
        # See if carbonate sediment is being generated at the *current* time (if bathymetry above CCD).
        # This mask determines where on the current ocean floor carbonate is currently being deposited and
        # hence does not include the deposition history over the lifetime of the current parcel of ocean floor.
        if bathymetry > ccd_at_current_time:
            carbonate_deposition_mask = 1.0
            lon_lat_carbonate_deposition_mask_list.append((lon, lat, carbonate_deposition_mask))
    
    return (lon_lat_carbonate_decompacted_sediment_thickness_list,
            lon_lat_carbonate_compacted_sediment_thickness_list,
            lon_lat_carbonate_deposition_mask_list)


def predict_sedimentation_and_write_data(
        input_points,
        time,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
        grid_spacing,
        age_to_depth_curve,
        dynamic_topography_model_name,  # Can be None
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
        carbonate_deposition_mask_filename_prefix):
    
    print('Time: ', time)
    
    # Read the CCD curve from the CCD file.
    # Returned curve is a function that accepts time and return depth.
    ccd_curve = read_curve(ccd_curve_filename)
    
    # Read the maximum carbonate decompacted sedimentation rate curve from the file.
    # Returned curve is a function that accepts time and return sedimentation rate (in cm/ky).
    max_carbonate_decomp_sed_rate_cm_per_ky_curve = read_curve(max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename)
    
    # Age, distance and bathymetry filenames for the current time.
    age_grid_filename = age_grid_filename_prefix + '{0}.{1}'.format(time, age_grid_filename_extension)
    distance_filename = distance_filename_prefix + '{0}.{1}'.format(time, distance_filename_extension)
    bathymetry_filename = bathymetry_filename_prefix + '{0}.{1}'.format(time, bathymetry_filename_extension)
    
    # Predict carbonate decompacted and compacted sediment thickness at each input point that is in
    # the age, distance and bathmetry grids (in unmasked regions of all three grids).
    sediment_thickness_data = predict_sedimentation(
        input_points,
        age_grid_filename,
        distance_filename,
        bathymetry_filename,
        age_to_depth_curve,
        dynamic_topography_model_name,
        ccd_curve,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        time)
    if sediment_thickness_data is None:
        print('Ignoring time "{0}" - no grid points inside age, distance and bathymetry grids.'.format(time), file=sys.stderr)
        return
    (carbonate_decompacted_sediment_thickness_data,
     carbonate_compacted_sediment_thickness_data,
     carbonate_deposition_mask_data) = sediment_thickness_data

    # Carbonate decompacted and compacted sediment thickness filenames without the filename extension.
    carbonate_decompacted_sediment_thickness_base_filename = (
        carbonate_decompacted_sediment_thickness_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    carbonate_compacted_sediment_thickness_base_filename = (
        carbonate_compacted_sediment_thickness_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    carbonate_deposition_mask_base_filename = (
        carbonate_deposition_mask_filename_prefix + '_{0}_{1}'.format(grid_spacing, time))
    
    # Output Carbonate decompacted and compacted sediment thicknesses to '.xy' files and '.nc' files.
    write_data(
        carbonate_decompacted_sediment_thickness_data,
        carbonate_decompacted_sediment_thickness_base_filename,
        grid_spacing, latitude_range, longitude_range)
    write_data(
        carbonate_compacted_sediment_thickness_data,
        carbonate_compacted_sediment_thickness_base_filename,
        grid_spacing, latitude_range, longitude_range)
    write_data(
        carbonate_deposition_mask_data,
        carbonate_deposition_mask_base_filename,
        grid_spacing, latitude_range, longitude_range)


#
# Parallel optimization
#


# Wraps around 'predict_sedimentation_and_write_data()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def predict_sedimentation_and_write_data_parallel_pool_function(args):
    try:
        return predict_sedimentation_and_write_data(*args)
    except KeyboardInterrupt:
        pass


def low_priority():
    """ Set the priority of the process to below-normal."""

    import sys
    try:
        sys.getwindowsversion()
    except AttributeError:
        isWindows = False
    else:
        isWindows = True

    if isWindows:
        import psutil
        
        p = psutil.Process()
        p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    else:
        import os

        os.nice(1)


def predict_sedimentation_and_write_data_for_times(
        times,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
        grid_spacing,
        age_to_depth_curve,
        dynamic_topography_model_name,  # Can be None
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
        use_all_cpu_cores):
    
    print('Starting...')
    
    # Generate a uniform global grid of lat/lon points.
    input_points = generate_input_points_grid(grid_spacing, latitude_range, longitude_range)

    # Either distribute each time iteration (over array of times) across all CPU cores, or
    # run time iteration loop serially (ie, use only one CPU core).
    if use_all_cpu_cores:
        
        # Split the workload across the CPUs.
        try:
            pool = multiprocessing.Pool(initializer=low_priority)
            pool_map_async_result = pool.map_async(
                predict_sedimentation_and_write_data_parallel_pool_function,
                (
                    (
                        input_points,
                        time,
                        latitude_range,
                        longitude_range,
                        grid_spacing,
                        age_to_depth_curve,
                        dynamic_topography_model_name,
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
                        carbonate_deposition_mask_filename_prefix
                    ) for time in times
                ),
                1)  # chunksize

            # Apparently if we use pool.map_async instead of pool.map and then get the results
            # using a timeout, then we avoid a bug in Python where a keyboard interrupt does not work properly.
            # See http://stackoverflow.com/questions/1408356/keyboard-interrupts-with-pythons-multiprocessing-pool
            try:
                pool_map_async_result.get(99999)
            except KeyboardInterrupt:
                # Note: 'finally' block below gets executed before returning.
                pass
        
        finally:
            pool.close()
            pool.join()

    else:  # process in serial...

        # Iterate over times and generate grids for each time.
        for time in times:
            # Predict carbonate decompacted and compacted sediment thickness, and write results to output grids.
            predict_sedimentation_and_write_data(
                input_points,
                time,
                latitude_range,
                longitude_range,
                grid_spacing,
                age_to_depth_curve,
                dynamic_topography_model_name,
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
                carbonate_deposition_mask_filename_prefix)
    
    print('...finished.')
