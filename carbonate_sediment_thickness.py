
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


##############################################################################
# Generate carbonate sediment thickness grids from age and bathymetry grids. #
##############################################################################



from call_system_command import call_system_command
import csv
import math
import multiprocessing
import os.path
# We'll use the slower path instead of using the faster PlateTectonicTools to assign plate IDs.
# This removes our dependency on PlateTectonicTools.
# The increased time does not affect the total running time very much
# (most of the time is spent calculating bathymetry from subsidence model).
USE_PTT = False
if USE_PTT:
    from ptt.utils import points_in_polygons
import pygplates
from scipy.interpolate import interp1d
import sys


#
# Input Parameters.
#


# Average ocean floor sediment density (kg/m^3), porosity (%) and decay (metres).
AVERAGE_OCEAN_FLOOR_SEDIMENT_DENSITY = 2647.0
AVERAGE_OCEAN_FLOOR_SEDIMENT_POROSITY = 0.66
AVERAGE_OCEAN_FLOOR_SEDIMENT_DECAY = 1333.0

# Dissolution distance (in metres) above the CCD.
# This distance above the CCD is where carbonate dissolution is 100%.
# Below this distance dissolution decreases linearly to 0% at the CCD.
CCD_DISSOLUTION_DISTANCE = 300.0


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


def gmt_grdtrack(
        input,
        grid_filename):
    """
    Samples a grid file at the specified locations.
    
    'input' is a list of (longitude, latitude) tuples where latitude and longitude are in degrees.
    
    Returns a list of float values.
    """
    
    # Create a multiline string (one line per lon/lat row).
    location_data = ''.join(
            ' '.join(str(item) for item in row) + '\n' for row in input)

    # The command-line strings to execute GMT 'grdtrack'.
    grdtrack_command_line = ["gmt", "grdtrack",
        "-G{0}".format(grid_filename),
        # Geographic input/output coordinates...
        "-fg",
        # Use linear interpolation, and avoid anti-aliasing...
        "-nl+a+bg+t0.5"]
    
    # Call the system command.
    stdout_data = call_system_command(grdtrack_command_line, stdin=location_data, return_stdout=True)

    # Extract the sampled values.
    output_values = []
    for line in stdout_data.splitlines():
        # Each line returned by GMT grdtrack contains "longitude latitude grid_value".
        # Note that if GMT returns "NaN" then we'll return float('nan').
        output_value = float(line.split()[2])
        output_values.append(output_value)
    
    return output_values


def read_curve(curve_filename):
    
    x_array = []
    y_array = []
    xmin, xmax = None, None
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

            x = float(row[0])
            y = float(row[1])
            
            x_array.append(x)
            y_array.append(y)
            
            # Track the y values of the two endpoints (ie, min/max x).
            if xmin is None or x < xmin:
                xmin = x
                y_at_xmin = y
            if xmax is None or x > xmax:
                xmax = x
                y_at_xmax = y
    
    # Convert curve to a piecewise linear interpolating 1D function that accepts x and returns y.
    # Return the y value at the boundary on bounds error instead of raising ValueError...
    return interp1d(x_array, y_array, bounds_error=False, fill_value=(y_at_xmin, y_at_xmax))


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


#
# Function to compact sediment thickness.
#

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


def calc_carbonate_decompacted_sediment_thickness(
        input_points,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        topology_filenames,
        rotation_model,
        ccd_curve,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        carbonate_anchor_plate_id,
        time):
    """
    Calculates carbonate decompacted sediment thickness.
    
    Returns a list of tuples (lon, lat, age, bathymetry, carbonate_decompacted_sediment_thickness),
    Each tuple corresponds to an input point that has a valid age and a valid bathymetry.
    """

    # Convert points from (lon, lat) to pygplates.PointOnSphere.
    points = [pygplates.PointOnSphere(lat, lon) for lon, lat in input_points]
    
    # Resolve the topologies for the current 'time'.
    # We do this in the carbonate reference frame because the input points are in the carbonate reference frame.
    resolved_topologies = []
    pygplates.resolve_topologies(topology_filenames, rotation_model, resolved_topologies, time, anchor_plate_id=carbonate_anchor_plate_id)

    # Assign a plate ID to each point using the resolved topologies.
    # We'll need plate IDs to reconstruct each point prior to sampling reconstructed bathymetry.
    resolved_topology_boundaries = [resolved_topology.get_resolved_boundary()
            for resolved_topology in resolved_topologies]
    resolved_topology_plate_ids = [resolved_topology.get_feature().get_reconstruction_plate_id()
            for resolved_topology in resolved_topologies]
    if USE_PTT:
        point_plate_ids = points_in_polygons.find_polygons(points, resolved_topology_boundaries, resolved_topology_plate_ids)
        # Default to zero plate ID for any points that happen to fall outside all topologies (eg, in a sliver crack).
        for point_index, plate_id in enumerate(point_plate_ids):
            if plate_id is None:
                point_plate_ids[point_index] = 0
    else:
        # Default to zero plate ID for any points that happen to fall outside all topologies (eg, in a sliver crack).
        point_plate_ids = [0] * len(points)
        for point_index, point in enumerate(points):
            for resolved_topology_index, resolved_topology_boundary in enumerate(resolved_topology_boundaries):
                if resolved_topology_boundary.is_point_in_polygon(point):
                    point_plate_ids[point_index] = resolved_topology_plate_ids[resolved_topology_index]
                    break

    # Convert the points from the carbonate reference frame to the default reference frame so
    # that we can sample the age and bathymetry grids at 'time' (using the default reference frame).
    input_points_in_default_reference_frame = []
    for point_index, (input_longitude, input_latitude) in enumerate(input_points):
        plate_id = point_plate_ids[point_index]
        # Get the rotation from 'time' to present day (using anchor plate 'carbonate_anchor_plate_id') and then
        # from present day back to 'time' (using anchor plate zero).
        rotation = (
            rotation_model.get_rotation(time, plate_id, 0) *
            rotation_model.get_rotation(0, plate_id, time, anchor_plate_id=carbonate_anchor_plate_id)
        )
        input_point = rotation * pygplates.PointOnSphere(input_latitude, input_longitude)
        input_lat, input_lon = input_point.to_lat_lon()
        input_points_in_default_reference_frame.append((input_lon, input_lat))
    
    age_grid_filename_prefix, age_grid_filename_decimal_places_in_time, age_grid_filename_extension = age_grid_filename_components
    bathymetry_filename_prefix, bathymetry_filename_decimal_places_in_time, bathymetry_filename_extension = bathymetry_filename_components
    
    # Age, distance and bathymetry filenames for the current time.
    age_grid_filename = age_grid_filename_prefix + '{0:.{1}f}.{2}'.format(time, age_grid_filename_decimal_places_in_time, age_grid_filename_extension)
    bathymetry_filename = bathymetry_filename_prefix + '{0:.{1}f}.{2}'.format(time, bathymetry_filename_decimal_places_in_time, bathymetry_filename_extension)
    
    # Get the ages and bathymetries at the input point locations (in default reference frame).
    ages = gmt_grdtrack(input_points_in_default_reference_frame, age_grid_filename)
    bathymetries = gmt_grdtrack(input_points_in_default_reference_frame, bathymetry_filename)

    # Initialise all carbonate deposition to zero.
    # We will accumulate carbonate deposition only above the CCD.
    carbonate_decompacted_sediment_thicknesses = [0.0] * len(input_points)

    # Only reference points where there are both age *and* bathymetry values.
    point_indices = [point_index
                     for point_index in range(len(input_points))
                     if not math.isnan(ages[point_index]) and not math.isnan(bathymetries[point_index])]

    # Return early if no points have both age *and* bathymetry values.
    if len(point_indices) == 0:
        return

    # Time interval to reconstruct ocean points (from 'time' to each ocean point's time of appearance).
    # Currently set to 1 Myr.
    time_interval = 1.0
    
    reconstruction_time = time + time_interval
    reconstructed_ages = [(age - time_interval) for age in ages]
    reconstructed_point_indices = point_indices[:]
    while reconstruction_time <= bathymetry_filename_oldest_time:
        # print(reconstruction_time); sys.stdout.flush()

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
        
        # Reconstruct each point from 'time' to 'reconstruction_time'.
        #
        # And convert them from the carbonate reference frame to the default reference frame so that
        # we can sample the bathymetry grids at 'reconstruction_time' (using the default reference frame).
        reconstructed_locations = []
        for reconstructed_point_index in reconstructed_point_indices:
            plate_id = point_plate_ids[reconstructed_point_index]
            # Get the rotation from 'time' to present day (using anchor plate 'carbonate_anchor_plate_id') and then
            # from present day to 'reconstruction_time' (using anchor plate zero).
            rotation = (
                rotation_model.get_rotation(reconstruction_time, plate_id, 0) *
                rotation_model.get_rotation(0, plate_id, time, anchor_plate_id=carbonate_anchor_plate_id)
            )
            reconstructed_point = rotation * points[reconstructed_point_index]
            reconstructed_lat, reconstructed_lon = reconstructed_point.to_lat_lon()
            reconstructed_locations.append((reconstructed_lon, reconstructed_lat))
        
        # Reconstructed bathymetry filename for the current reconstruction time.
        bathymetry_filename_prefix, bathymetry_filename_decimal_places_in_time, bathymetry_filename_extension = bathymetry_filename_components
        reconstructed_bathymetry_filename = bathymetry_filename_prefix + '{0:.{1}f}.{2}'.format(
                reconstruction_time, bathymetry_filename_decimal_places_in_time, bathymetry_filename_extension)

        # Sample reconstructed bathymetry at the reconstructed locations.
        reconstructed_bathymetrys = gmt_grdtrack(reconstructed_locations, reconstructed_bathymetry_filename)
        
        # The CCD depth associated with the reconstructed time.
        ccd_depth_at_reconstruction_time = ccd_curve(reconstruction_time)

        # Reconstruct bathymetry for all ocean points that exist at the current reconstruction time.
        for index, reconstructed_point_index in enumerate(reconstructed_point_indices):
            
            # Sample bathymetry at current reconstruction time.
            #
            # Note the use of 'index' instead of 'reconstructed_point_index'
            bathymetry_at_reconstruction_time = reconstructed_bathymetrys[index]

            # print('  ', time, reconstruction_time, bathymetry_at_reconstruction_time); sys.stdout.flush()
            
            # If modelled bathymetry is above the CCD depth then add a 1 Myr time interval to the time spent above CDD.
            if (not math.isnan(bathymetry_at_reconstruction_time) and
                bathymetry_at_reconstruction_time > ccd_depth_at_reconstruction_time):
                
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

    # Combine carbonate decompacted sediment thickness with the lons, lats, ages and bathymetries.
    # And only include those points that have both age *and* bathymetry values (at 'time').
    return [input_points[point_index] + (ages[point_index], bathymetries[point_index], carbonate_decompacted_sediment_thicknesses[point_index])
            for point_index in point_indices]


def calc_sedimentation(
        input_points,  # List of (lon, lat) tuples,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        topology_filenames,
        rotation_model,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        carbonate_anchor_plate_id,
        time):
    """
    Calculates decompacted and compacted carbonate thickness for each ocean basin grid point.
    Also determines when carbonate is currently being deposited (ie, only at the current time).
    
    Returns: A 3-tuple of lists each containing 3-tuples of (lon, lat, value).
             The values in the first list contain decompacted thickness.
             The values in the second list contain compacted thickness.
             The values in the third list contain a mask (value 1.0) for regions where
             carbonate is currently being deposited (ie, only at the current time).
    """
    
    # Read the CCD curve from the CCD file.
    # Returned curve is a function that accepts time and return depth.
    ccd_curve = read_curve(ccd_curve_filename)
    
    # Read the maximum carbonate decompacted sedimentation rate curve from the file.
    # Returned curve is a function that accepts time and return sedimentation rate (in cm/ky).
    max_carbonate_decomp_sed_rate_cm_per_ky_curve = read_curve(max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename)
    
    # Predict *decompacted* thickness.
    # We do this for all ocean basin points together since it's more efficient to sample reconstructed bathymetry for all points at once.
    lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list = calc_carbonate_decompacted_sediment_thickness(
        input_points,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        topology_filenames, rotation_model,
        ccd_curve, max_carbonate_decomp_sed_rate_cm_per_ky_curve,
        carbonate_anchor_plate_id,
        time)
    if not lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list:
        return
    
    # The CCD at the current time.
    ccd_at_current_time = ccd_curve(time)
    
    # For each ocean basin point also generate *compacted* thickness.
    lon_lat_carbonate_decompacted_sediment_thickness_list = []
    lon_lat_carbonate_compacted_sediment_thickness_list = []
    lon_lat_carbonate_deposition_mask_list = []
    for lon, lat, age, bathymetry, carbonate_decompacted_sediment_thickness in lon_lat_age_bathymetry_carbonate_decompacted_sediment_thickness_list:
        
        # Compact thickness.
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


def calc_sedimentation_and_write_data(
        input_points,
        time,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
        grid_spacing,
        topology_model_name,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        carbonate_decompacted_sediment_thickness_filename_prefix,
        carbonate_compacted_sediment_thickness_filename_prefix,
        carbonate_deposition_mask_filename_prefix,
        carbonate_anchor_plate_id):
    
    print('Time: ', time)

    # Read filenames listed in topology list file.
    topology_list_filename = os.path.join('input_data', 'topology_model', topology_model_name, 'topology_files.txt')
    with open(topology_list_filename, 'r') as topology_list_file:
        topology_filenames = topology_list_file.read().splitlines()

    # Read filenames listed in rotation list file.
    rotation_list_filename = os.path.join('input_data', 'topology_model', topology_model_name, 'rotation_files.txt')
    with open(rotation_list_filename, 'r') as rotation_list_file:
        rotation_filenames = rotation_list_file.read().splitlines()

    # Create the rotation model from the rotation files.
    rotation_model = pygplates.RotationModel(rotation_filenames)
    
    # Calculate carbonate decompacted and compacted sediment thickness at each input point that is in
    # the age and bathmetry grids (in unmasked regions of both grids).
    sediment_thickness_data = calc_sedimentation(
        input_points,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        topology_filenames,
        rotation_model,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        carbonate_anchor_plate_id,
        time)
    if sediment_thickness_data is None:
        print('Ignoring time "{0}" - no grid points inside age and bathymetry grids.'.format(time), file=sys.stderr)
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


# Wraps around 'calc_sedimentation_and_write_data()' so can be used by multiprocessing.Pool.map()
# which requires a single-argument function.
def calc_sedimentation_and_write_data_parallel_pool_function(args):
    try:
        return calc_sedimentation_and_write_data(*args)
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
        try:
            import psutil
        except ImportError:
            pass
        else:
            p = psutil.Process()
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
    else:
        import os

        os.nice(1)


def calc_sedimentation_and_write_data_for_times(
        times,
        latitude_range,  # (min, max) tuple
        longitude_range, # (min, max) tuple
        grid_spacing,
        topology_model_name,
        ccd_curve_filename,
        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
        age_grid_filename_components,
        bathymetry_filename_components,
        bathymetry_filename_oldest_time,
        carbonate_decompacted_sediment_thickness_filename_prefix,
        carbonate_compacted_sediment_thickness_filename_prefix,
        carbonate_deposition_mask_filename_prefix,
        carbonate_anchor_plate_id,
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
                calc_sedimentation_and_write_data_parallel_pool_function,
                (
                    (
                        input_points,
                        time,
                        latitude_range,
                        longitude_range,
                        grid_spacing,
                        topology_model_name,
                        ccd_curve_filename,
                        max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
                        age_grid_filename_components,
                        bathymetry_filename_components,
                        bathymetry_filename_oldest_time,
                        carbonate_decompacted_sediment_thickness_filename_prefix,
                        carbonate_compacted_sediment_thickness_filename_prefix,
                        carbonate_deposition_mask_filename_prefix,
                        carbonate_anchor_plate_id
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
            calc_sedimentation_and_write_data(
                input_points,
                time,
                latitude_range,
                longitude_range,
                grid_spacing,
                topology_model_name,
                ccd_curve_filename,
                max_carbonate_decomp_sed_rate_cm_per_ky_curve_filename,
                age_grid_filename_components,
                bathymetry_filename_components,
                bathymetry_filename_oldest_time,
                carbonate_decompacted_sediment_thickness_filename_prefix,
                carbonate_compacted_sediment_thickness_filename_prefix,
                carbonate_deposition_mask_filename_prefix,
                carbonate_anchor_plate_id)
    
    print('...finished.')
