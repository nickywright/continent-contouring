# Continent contouring

Included are some workflows that use continent contouring from [`PlateTectonicTools`](https://github.com/EarthByte/PlateTectonicTools). This essentially amounts to creating a `ContinentContouring` object using rotation files, some continent/craton features (polygons), a contour resolution, a gap threshold, an area threshold and an age range. And then asking it to reconstruct the polygons to an `age` and contour them into continents:

```
from ptt.continent_contours import ContinentContouring

continent_contouring = ContinentContouring(...)
...
contoured_continents = continent_contouring.get_contoured_continents(age)
```

The workflows then use these contours in different ways.

## Installation

### Dependencies

The following Python packages are required:

- [`gplately`](https://github.com/GPlates/gplately)
- [`numpy`](http://numpy.org)
- [`PlateTectonicTools`](https://github.com/EarthByte/PlateTectonicTools)
- [`pygplates`](http://gplates.org/docs/pygplates/pygplates_getting_started.html#installation)

You can install these with conda:

```
conda create -n <conda-environment> -c conda-forge gplately numpy platetectonictools pygplates
conda activate <conda-environment>
```

...but until version 0.5 of PlateTectonicTools is available you'll need to install `platetectonictools` from Github:

```
conda create -n <conda-environment> -c conda-forge gplately numpy pygplates
conda activate <conda-environment>
conda install git pip
pip install git+https://github.com/EarthByte/PlateTectonicTools
```

...where `<conda-environment>` should be replaced with the name of your conda environment.

## Documentation

Current workflows include:

### create_passive_margins.py

Creates passive margins through time by excluding those parts of contours that are active considered to be margins (within a threshold distance from subduction zones).

To run the workflow you'll first need to copy the [Merdith et al (2021) 1Ga plate model files](https://github.sydney.edu.au/EarthByte/EarthBytePlateMotionModel-ARCHIVE/tree/master/Merdith_etal_2021_Published) into `models/Merdith_etal_2021_Published/`.

> _Note:_ You can use a different plate model by changing the parameters `rotation_features`, `topology_features` and `continent_features` to refer to files that you provide. And you can change the time range.

You can then run the workflow with `python create_passive_margins.py`. You should end up with a files `passive_margin_features.gpml` containing the passive margin features, `contoured_features.gpml` containing the continental-oceanic boundary (COB) features and `subduction_zone_features.gpml` containing the subduction zones. 

Parameters at the top of `create_passive_margins.py` control its behaviour. The time range is specified with `start_time`, `end_time` and `time_interval`. For example, you may want to change `time_interval` to 1 (instead of 5) to get a snapshot every 1Myr. The maximum distance of an active margin from a subduction zone (in kms) is specified with `max_distance_of_subduction_zone_from_active_margin_kms`. The grid spacing (in degrees) between points in the grid used for contouring/aggregrating blocks of continental polygons is specified with `continent_contouring_point_spacing_degrees`. The function `continent_contouring_area_threshold_steradians(time)` determines if contour polygons smaller than this will be excluded when contouring/aggregrating blocks of continental polygons. The function `continent_contouring_gap_threshold_radians(time)` determines if gaps between continent polygons smaller than this will be excluded when contouring/aggregrating blocks of continental polygons. With the functions you can write your own to return thresholds based on the `time` argument passed into the function. And finally, the `use_all_cpus` parameter controls CPU usage (whether to use multiple CPUs and how many).

All the passive margins found for all times end up in a single file. Each passive margin feature only exists at its reconstruction time (+/- 0.5 * `time_interval`). And you don’t need a rotation file to view these through time (they don’t have plate IDs) since they are meant to be reconstructed snapshots (from 0Ma to 1000Ma).

> _Note:_ When running the script you’ll get lots of temporary files, but they’ll get cleaned up when the script finishes (provided it’s not interrupted).
