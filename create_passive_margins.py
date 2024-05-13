# To run this script: python create_passive_margins.py <path to yaml file>


from functools import partial
import gplately
import math
import multiprocessing
import numpy as np
import os
import os.path
import pygplates
from gplately.ptt import continent_contours
import sys

from datetime import datetime
import yaml
try:
    from yaml import Cloader as Loader
except ImportError:
    from yaml import Loader

# 2024-05-09: Modified this script to be compatible with a yaml file and output to user specified location (NW)

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ------------------------------------------
# --- Set paths and various parameters
# ------------------------------------------
today = datetime.today().strftime('%Y%m%d')

try:
    config_file = sys.argv[1]
    print("*** Parameters set from %s ***" % config_file)
    with open(config_file) as f:
        PARAMS = yaml.load(f, Loader=Loader)

    # ------------------------------------------
    # --- Set directories and input files ------
    model_name = PARAMS["InputFiles"]["model_name"]
    paleobathymetry_main_output_dir = PARAMS["OutputFiles"]["paleobathymetry_main_output_dir"]
    include_date_in_output_dir = PARAMS["OutputFiles"]["include_date_in_output_dir"]
    date = PARAMS["OutputFiles"]["date"]

    sediment_thickness_output_dir = PARAMS["OutputFiles"]["sediment_thickness_output_dir"]
    sediment_thickness_within_main_output_dir = PARAMS["OutputFiles"]["sediment_thickness_within_main_output_dir"]

    # --- agegrids
    agegrid_dir = PARAMS["InputFiles"]["agegrid_dir"]
    agegrid_filename = PARAMS["InputFiles"]["agegrid_filename"]
    agegrid_filename_ext = PARAMS["InputFiles"]["agegrid_filename_ext"]

    # --- input file
    proximity_features_files = [PARAMS["InputFiles"]["sediment_thickness_features"]]

    # --- Plate model files
    plate_model_dir = PARAMS["InputFiles"]["model_dir"]
    rotation_filenames = [os.path.join(plate_model_dir, i) for i in PARAMS['InputFiles']['rotation_files']]
    topology_filenames = [os.path.join(plate_model_dir, i) for i in PARAMS['InputFiles']['topology_files']]
    coastline_filename = '%s/%s' % (plate_model_dir, PARAMS['InputFiles']['coastline_file'])
    
    continents_file = '%s/%s' % (plate_model_dir, PARAMS['InputFiles']['continents_file'])
    cratons_file = '%s/%s' % (plate_model_dir, PARAMS['InputFiles']['cratons_file'])

    anchor_plate_id = PARAMS["InputFiles"]["anchor_plate_id"]

    # --- grid spacing
    grid_spacing = PARAMS["GridParameters"]["grid_spacing"]
    lon_min = PARAMS["GridParameters"]["lon_min"]
    lon_max = PARAMS["GridParameters"]["lon_max"]
    lat_min = PARAMS["GridParameters"]["lat_min"]
    lat_max = PARAMS["GridParameters"]["lat_max"]

    # --- time parameters
    min_time = int(PARAMS["TimeParameters"]["time_min"])           # Not truly a min_time, - parts 1 and 3 require a 0 Ma shapefile
    # oldest time to reconstruct to (will default to 0 Ma for min time)
    max_time = int(PARAMS["TimeParameters"]["time_max"])
    time_step = int(PARAMS["TimeParameters"]["time_step"])      # Myrs to increment age by in loop
    
    # running parameters
    num_cpus = PARAMS["Parameters"]["number_of_cpus"] # number of cpus to use. Reduce if required!
    max_memory_usage_in_gb = PARAMS["Parameters"]["max_memory_usage_in_gb_sedthickness"]

except IndexError:
    print('*** No yaml file given. Make sure you specify it ***')

# ------------------------------
# --- base output directory
if include_date_in_output_dir.lower() in ['true', '1', 't', 'y', 'yes']:
    if date == 'today':
        date = today
        paleobathymetry_main_output_dir = '%s/%s' % (paleobathymetry_main_output_dir, date)
    else:
        paleobathymetry_main_output_dir = '%s/%s' % (paleobathymetry_main_output_dir, date)
else:
    paleobathymetry_main_output_dir = '%s' % (paleobathymetry_main_output_dir)

output_dir = '%s/continent_contouring' % (paleobathymetry_main_output_dir)

os.makedirs(output_dir, exist_ok=True)

#############################
# Start of input parameters #
#############################

# Location of plate model.
model_dir = plate_model_dir

# Rotation files (relative to input directory).
rotation_features = rotation_filenames

# The reference frame to generate the output files relative to.
anchor_plate_id = anchor_plate_id

# Topology features (absolute file paths).
#
# Only include those GPML files that are used for topologies.
topology_features = topology_filenames

# Continent polygon features (absolute file paths).
continent_features = [continents_file, cratons_file]
# print(continent_features)
# Time range.
start_time = min_time
end_time = max_time
time_interval = time_step
times = np.arange(start_time, end_time + 0.5 * time_interval, time_interval)

# Use all CPUs.
#
# If False then use a single CPU.
# If True then use all CPUs (cores) - and make sure you don't interrupt the process.
# If a positive integer then use that specific number of CPUs (cores).
#
#use_all_cpus = False
#use_all_cpus = 4
use_all_cpus = num_cpus

# Maximum distance of an active margin from a subduction zone.
#
# If a contoured continent segment is near any subduction zone (ie, within this distance) then it's an active margin (ie, not a passive margin).
max_distance_of_subduction_zone_from_active_margin_kms = 500
max_distance_of_subduction_zone_from_active_margin_radians = max_distance_of_subduction_zone_from_active_margin_kms / pygplates.Earth.mean_radius_in_kms

# The grid spacing (in degrees) between points in the grid used for contouring/aggregrating blocks of continental polygons.
continent_contouring_point_spacing_degrees = 0.25

# Optional parameter specifying a minimum area threshold (in square radians) for including contoured continents.
#
# Contoured continents with area smaller than this threshold will be excluded.
# If this parameter is not specified then no area threshold is applied.
#
# Can also be a function (accepting time in Ma) and returning the area threshold.
#
# Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
#       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
#       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
continent_contouring_area_threshold_square_kms = 0
continent_contouring_area_threshold_steradians = continent_contouring_area_threshold_square_kms / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)

# Optional parameter specifying a minimum area threshold (in square radians) for contours that exclude continental crust.
#
# Polygon contours that exclude continental crust and have an area smaller than this threshold will be excluded
# (meaning they will now *include* continental crust, thus removing the contour).
# This is useful for removing small holes inside continents.
# If this parameter is not specified then no area threshold is applied.
#
# Can also be a function (accepting time in Ma) and returning the area threshold.
#
# Note: Units here are for normalised sphere (ie, steradians or square radians) so full Earth area is 4*pi.
#       So 0.1 covers an area of approximately 4,000,000 km^2 (ie, 0.1 * 6371^2, where Earth radius is 6371km).
#       Conversely 4,000,000 km^2 is equivalent to (4,000,000 / 6371^2) steradians.
continent_exclusion_area_threshold_square_kms = 800000
continent_exclusion_area_threshold_steradians = continent_exclusion_area_threshold_square_kms / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)

# Optional parameter specifying the distance threshold (in radians) above which continents are separated.
#
# Any continent polygons separated by a distance that is less than this threshold will become part of the same continent.
#
# Can also be a function (accepting time in Ma) and returning the distance threshold.
#
# Note: Units here are for normalised sphere (ie, radians).
#       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
#       Also 1.0 degree is approximately 110 km.
#
#continent_separation_distance_threshold_kms = 0
#continent_separation_distance_threshold_radians = continent_separation_distance_threshold_kms / pygplates.Earth.mean_radius_in_kms
continent_separation_distance_threshold_radians = continent_contours.DEFAULT_CONTINENT_SEPARATION_DISTANCE_THRESHOLD_RADIANS

# Optional parameter specifying a distance (in radians) to expand contours ocean-ward - this also
# ensures small gaps between continents are ignored during contouring.
#
# The continent(s) will be expanded by a buffer of this distance (in radians) when contouring/aggregrating blocks of continental polygons.
# If this parameter is not specified then buffer expansion is not applied.
#
# This parameter can also be a function (that returns the distance).
# The function can have a single function argument: (1) accepting time (in Ma).
# Or it can have two function arguments: (1) the first accepting time (in Ma) and
# (2) the second accepting the contoured continent (a 'ptt.continent_contours.ContouredContinent' object)
# of the (unexpanded) contoured continent that the buffer/gap distance will apply to.
# Hence a function with *two* arguments means a different buffer/gap distance can be specified for each contoured continent (eg, based on its area).
#
# Note: Units here are for normalised sphere (ie, radians).
#       So 1.0 radian is approximately 6371 km (where Earth radius is 6371 km).
#       Also 1.0 degree is approximately 110 km.
def continent_contouring_buffer_and_gap_distance_radians(time, contoured_continent):
    # One distance for time interval [1000, 300] and another for time interval [250, 0].
    # And linearly interpolate between them over the time interval [300, 250].
    pre_pangea_distance_radians = math.radians(2.5)  # convert degrees to radians
    post_pangea_distance_radians = math.radians(0.0)  # convert degrees to radians
    if time > 300:
        buffer_and_gap_distance_radians = pre_pangea_distance_radians
    elif time < 250:
        buffer_and_gap_distance_radians = post_pangea_distance_radians
    else:
        # Linearly interpolate between 250 and 300 Ma.
        interp = float(time - 250) / (300 - 250)
        buffer_and_gap_distance_radians = interp * pre_pangea_distance_radians + (1 - interp) * post_pangea_distance_radians
    
    # Area of the contoured continent.
    area_steradians = contoured_continent.get_area()

    # Linearly reduce the buffer/gap distance for contoured continents with area smaller than 1 million km^2.
    area_threshold_square_kms = 500000
    area_threshold_steradians = area_threshold_square_kms / (pygplates.Earth.mean_radius_in_kms * pygplates.Earth.mean_radius_in_kms)
    if area_steradians < area_threshold_steradians:
        buffer_and_gap_distance_radians *= area_steradians / area_threshold_steradians

    return buffer_and_gap_distance_radians

###########################
# End of input parameters #
###########################



# Rotation model.
rotation_model = pygplates.RotationModel(rotation_features, default_anchor_plate_id=anchor_plate_id)

# Create a ContinentContouring object.
continent_contouring = continent_contours.ContinentContouring(
        rotation_model,
        continent_features,
        continent_contouring_point_spacing_degrees,
        continent_contouring_area_threshold_steradians,
        continent_contouring_buffer_and_gap_distance_radians,
        continent_exclusion_area_threshold_steradians,
        continent_separation_distance_threshold_radians)

# Find passive margins at the specified time.
def find_passive_margins(
        time):
    print('time:', time)
    
    continent_contour_features = []
    passive_margin_features = []
    subduction_zone_features = []
    
    # Resolve the topological plate polygons for the current time.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time, shared_boundary_sections)
    
    # Get the resolved subduction zone segments.
    subduction_zone_lines = []
    for shared_boundary_section in shared_boundary_sections:
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
            shared_sub_segment_lines = [shared_sub_segment.get_resolved_geometry() for shared_sub_segment in shared_boundary_section.get_shared_sub_segments()]
            subduction_zone_lines.extend(shared_sub_segment_lines)
            # Also save subduction zone lines as features (so we can later save them to a file for debugging).
            subduction_zone_feature = pygplates.Feature()
            subduction_zone_feature.set_geometry(shared_sub_segment_lines)
            subduction_zone_feature.set_valid_time(time + 0.5 * time_interval - 1e-4,  # epsilon to avoid overlap at interval boundaries
                                                   time - 0.5 * time_interval)
            subduction_zone_features.append(subduction_zone_feature)

    # Get the continent mask and the continent contours at the current time.
    continent_mask, contoured_continents = continent_contouring.get_continent_mask_and_contoured_continents(time)

    # Write out the continent mask as NetCDF.
    
    continent_mask_filename = '%s/continent_mask_%s.nc' % (output_dir, time)
    # Note that we need to convert the boolean mask grid to a non-boolean number type for NetCDF (and it seems floating-point for gplately).
    continent_mask_grid = continent_mask.astype('float')
    gplately.grids.write_netcdf_grid(continent_mask_filename, continent_mask_grid)
    
    # Convert all continent contour geometries into features.
    for contoured_continent in contoured_continents:
        for continent_contour_geometry in contoured_continent.get_contours():
            continent_contour_feature = pygplates.Feature()
            continent_contour_feature.set_geometry(continent_contour_geometry)
            continent_contour_feature.set_valid_time(time + 0.5 * time_interval - 1e-4,  # epsilon to avoid overlap at interval boundaries
                                                     time - 0.5 * time_interval)
            continent_contour_features.append(continent_contour_feature)
    
    # Find passive margins along continent contours by removing active margins (contoured segments close to a subduction zone).
    passive_margin_geometries = []
    for contoured_continent in contoured_continents:
        for contour_polyline in contoured_continent.get_contours():
            # Add any passive margins found on the current contour.
            _find_passive_margin_geometries_on_contour(passive_margin_geometries,
                                                       subduction_zone_lines,
                                                       contour_polyline.get_segments())
    
    # Convert any passive margin geometries found into features.
    for passive_margin_geometry in passive_margin_geometries:
        passive_margin_feature = pygplates.Feature()
        passive_margin_feature.set_geometry(passive_margin_geometry)
        passive_margin_feature.set_valid_time(time + 0.5 * time_interval - 1e-4,  # epsilon to avoid overlap at interval boundaries
                                              time - 0.5 * time_interval)
        passive_margin_features.append(passive_margin_feature)

    # Save continent contours to GPML.
    pygplates.FeatureCollection(continent_contour_features).write('%s/continent_contour_features_%s.gpml' % (output_dir, time))

    # Save passive margins to GPML.
    pygplates.FeatureCollection(passive_margin_features).write('%s/passive_margin_features_%s.gpml' % (output_dir, time))

    # Save subducton zone segments to GPML (for debugging).
    pygplates.FeatureCollection(subduction_zone_features).write('%s/subduction_zone_features_%s.gpml' % (output_dir, time))


def _find_passive_margin_geometries_on_contour(
        passive_margin_geometries,
        subduction_zone_lines,
        contoured_continent_segments):
    
    # Points for the current passive margin (if one).
    passive_margin_adjacent_points = []
    
    # Iterate over great circle arc segments of the contour.
    for contoured_continent_segment in contoured_continent_segments:
        # Create a polyline from the current contoured segment (so we can do distance testing).
        contoured_continent_line_segment = pygplates.PolylineOnSphere((
            contoured_continent_segment.get_start_point(),
            contoured_continent_segment.get_end_point()))
        
        # If continent segment near any subduction zone then it's an active margin.
        is_active_margin = False
        for subduction_zone_line in subduction_zone_lines:
            # If distance less than threshold distance.
            if pygplates.GeometryOnSphere.distance(
                    contoured_continent_line_segment, subduction_zone_line, max_distance_of_subduction_zone_from_active_margin_radians) is not None:
                is_active_margin = True
                break  # skip remaining subduction zones
        
        # If it's not an active margin then it's a passive margin.
        is_passive_margin = not is_active_margin
        
        # If segment a passive margin then add segment to current passive margin.
        if is_passive_margin:
            if not passive_margin_adjacent_points:
                # Add segment start point for first segment.
                passive_margin_adjacent_points.append(contoured_continent_line_segment[0])
            # Add segment end point.
            passive_margin_adjacent_points.append(contoured_continent_line_segment[1])
        else:  # active margin
            if passive_margin_adjacent_points:
                # We have accumulated passive margin points but are now in an active margin, so submit a passive margin.
                passive_margin_geometries.append(
                    pygplates.PolylineOnSphere(passive_margin_adjacent_points))
                # Clear points for next passive margin geometry
                del passive_margin_adjacent_points[:]
    
    # If there's one last passive margin geometry in the current contour then submit it.
    if passive_margin_adjacent_points:
        passive_margin_geometries.append(
            pygplates.PolylineOnSphere(passive_margin_adjacent_points))


if __name__ == '__main__':
    
    print('... output to %s:' % output_dir)

    if use_all_cpus:
    
        # If 'use_all_cpus' is a bool (and therefore must be True) then use all available CPUs...
        if isinstance(use_all_cpus, bool):
            try:
                num_cpus = multiprocessing.cpu_count()
            except NotImplementedError:
                num_cpus = 1
        # else 'use_all_cpus' is a positive integer specifying the number of CPUs to use...
        elif isinstance(use_all_cpus, int) and use_all_cpus > 0:
            num_cpus = use_all_cpus
        else:
            raise TypeError('use_all_cpus: {} is neither a bool nor a positive integer'.format(use_all_cpus))
        
        # Distribute writing of each grid to a different CPU.
        with multiprocessing.Pool(num_cpus) as pool:
            pool.map(
                    partial(find_passive_margins),
                    times,
                    1) # chunksize
    
    else:
        for time in times:
            find_passive_margins(time)
    
    # Combine all features from each 'time'.
    continent_contour_features =  []
    passive_margin_features = []
    subduction_zone_features = []
    for time in times:
        continent_contour_filename = '%s/continent_contour_features_%s.gpml' % (output_dir, time)
        continent_contour_features.extend(pygplates.FeatureCollection(continent_contour_filename))
        # Remove temporary file at current 'time'.
        if os.access(continent_contour_filename, os.R_OK):
            os.remove(continent_contour_filename)
        
        passive_margin_filename = '%s/passive_margin_features_%s.gpml' % (output_dir, time)
        passive_margin_features.extend(pygplates.FeatureCollection(passive_margin_filename))
        # Remove temporary file at current 'time'.
        if os.access(passive_margin_filename, os.R_OK):
            os.remove(passive_margin_filename)
        
        subduction_zone_filename = '%s/subduction_zone_features_%s.gpml' % (output_dir, time)
        subduction_zone_features.extend(pygplates.FeatureCollection(subduction_zone_filename))
        # Remove temporary file at current 'time'.
        if os.access(subduction_zone_filename, os.R_OK):
            os.remove(subduction_zone_filename)
    
    # Save ALL continent contours to GPMLZ.
    pygplates.FeatureCollection(continent_contour_features).write('%s/continent_contour_features.gpmlz' % output_dir)
    
    # Save ALL passive margins to GPMLZ.
    pygplates.FeatureCollection(passive_margin_features).write('%s/passive_margin_features.gpmlz' % output_dir)

    # Save ALL subducton zone segments to GPMLZ.
    pygplates.FeatureCollection(subduction_zone_features).write('%s/subduction_zone_features.gpmlz' % output_dir)
