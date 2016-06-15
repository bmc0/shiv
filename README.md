### About:

Shiv is a multithreaded G-code generator for 3D printers.

### Building:

#### Dependencies:

* C++11 compiler

#### Build:

For Unix-like systems:

	$ ./build.sh

Note: If you're compiling under Cygwin, remove `-fopenmp` from `$CXXFLAGS` in
`build.sh`. Libgomp doesn't seem to work very well under Windows so shiv often
ends up running slower with OpenMP enabled.

For Windows (make sure `cl` is in your `%PATH%`):

	> build.bat

OpenMP is enabled by default for both platforms.

### Usage:

#### Synopsis:

	shiv [-hp] [-o output_path] [-c config_path] [-O option=value]
	     [-l layer_height] [-w extrusion_width] [-t tolerance]
	     [-s scale_factor] [-d infill_density] [-n shells]
	     [-r roof_thickness] [-f floor_thickness] [-b brim_width]
	     [-C coarseness] [-x x_translate] [-y y_translate]
	     binary_stl_file

#### Flags:

Flag                 | Description
---------------------|-----------------------------------------
`-h`                 | Show help text
`-p`                 | Preview slices (pipe stdout to gnuplot)
`-o output_path`     | Output gcode path
`-c config_path`     | Configuration file path
`-O option=value`    | Set option to value
`-l layer_height`    | Layer height
`-w extrusion_width` | Constrained extrusion width
`-t tolerance`       | Segment connection tolerance
`-s scale_factor`    | Object scale factor
`-d infill_density`  | Sparse infill density
`-n shells`          | Number of shells
`-r roof_thickness`  | Solid roof thickness
`-f floor_thickness` | Solid floor thickness
`-b brim_width`      | Brim width
`-C coarseness`      | Output coarseness
`-x x_translate`     | Translate object in the x-axis
`-y y_translate`     | Translate object in the y-axis

#### Full options list:

Option                  | Default value | Description
------------------------|--------------:|------------
`layer_height`          |         `0.2` | Layer height.
`tolerance`             |       `0.001` | Segment connection tolerance. Large values can be used to close holes in broken models.
`coarseness`            |         `5.0` | Approximate output coarseness in `1/SCALE_CONSTANT` units (microns with the default `SCALE_CONSTANT` and millimeter input/output units).
`extrusion_width`       |         `0.4` | Constrained extrusion width. Should *generally* be set to a value similar to your nozzle diameter.
`xy_scale_factor`       |       `1.003` | The object is scaled by this ratio in the x and y axes to compensate for shrinkage. Around 1.003 works for PLA. ABS should be somewhere between 1.004 and 1.009.
`z_scale_factor`        |         `1.0` | The object is scaled by this ratio in the z axis to compensate for shrinkage. Should probably be left at 1 unless a high temperature heated chamber is used.
`x_center`              |         `0.0` | X coordinate to center the object on.
`y_center`              |         `0.0` | Y coordinate to center the object on.
`packing_density`       |        `0.99` | Solid packing density. Should be slightly less than 1. 0.99 seems to work for PLA.
`infill_density`        |         `0.2` | Sparse infill density.
`shells`                |           `2` | Number of loops/perimeters/shells (whatever you want to call them).
`roof_thickness`        |         `0.8` | Solid surface thickness when looking upwards.
`floor_thickness`       |         `0.8` | Solid surface thickness when looking downwards.
`material_diameter`     |        `1.75` | Diameter of material.
`flow_multiplier`       |         `1.0` | Flow rate adjustment to compensate for incorrect E-steps or E-step variation between material types.
`perimeter_feed_rate`   |        `25.0` | Outer shell feed rate.
`loop_feed_rate`        |        `40.0` | Inner shell(s) feed rate.
`infill_feed_rate`      |        `50.0` | Infill feed rate.
`travel_feed_rate`      |       `120.0` | Travel feed rate.
`first_layer_mult`      |         `0.5` | First layer feed rates (except travel) are multiplied by this value.
`retract_len`           |         `1.0` | Retraction length.
`retract_speed`         |        `20.0` | Retraction speed.
`restart_speed`         |        `-1.0` | Restart speed. -1 means same as `retract_speed`.
`retract_min_travel`    |         `1.6` | Minimum travel for retraction when not crossing a boundary. Has no effect if `retract_within_island` is false.
`retract_threshold`     |         `8.0` | Unconditional retraction threshold.
`retract_within_island` |       `false` | If false, retraction will not occur unless a boundary is crossed or the travel distance is greater than `retract_threshold`.
`wipe_len`              |         `0.0` | Extra travel distance at the end of a shell.
`cool_layer`            |           `1` | Turn on part cooling at this layer (numbered from zero). Set to a negative number to disable cooling.
`start_gcode`           |        `NULL` | Prepend this G-code to beginning of the output file.
`end_gcode`             |        `NULL` | Append this G-code to the end of the output file.
`cool_on_gcode`         |   `M106 S255` | G-code to turn on cooling.
`cool_off_gcode`        |        `M107` | G-code to turn off cooling.
`temp`                  |       `220.0` | Hotend temperature.
`bed_temp`              |        `65.0` | Bed temperature.
`infill_first`          |       `false` | Do infill before shells.
`align_seams`           |        `true` | Align seams to the lower left corner. The nearest point is picked instead if this is false.
`clean_insets`          |        `true` | Do `ClipperLib::CleanPolygon()` on all insets (only the initial outline is cleaned if this is false).
`fill_inset_gaps`       |        `true` | Fill gaps between shells.
`no_solid`              |       `false` | If true, only generate solid fill on the very top and bottom of the model.
`anchor`                |        `true` | Clip and anchor inset paths.
`outside_first`         |       `false` | Prefer exterior shells.
`solid_infill_first`    |       `false` | Print solid infill before sparse infill.
`separate_z_travel`     |       `false` | Generate a separate z travel move instead of moving all axes together.
`comb`                  |       `false` | Avoid crossing perimeters (disabled by default until more testing is done).
`combine_all`           |       `false` | Orients all outlines counter-clockwise. This can be used to fix certain broken models, but it also fills holes.
`generate_support`      |       `false` | Generate support structure.
`support_everywhere`    |       `false` | False means only touching build plate.
`solid_support_base`    |       `false` | Make supports solid at layer 0.
`poly_fill_type`        |    `non_zero` | Poly fill type for union. Sometimes `even_odd` is useful for broken models with self-intersections and/or incorrect normals.
`fill_threshold`        |         `0.2` | Infill and inset gap fill is removed when it would be narrower than `extrusion_width * fill_threshold`.
`support_angle`         |        `60.0` | Angle threshold for support.
`support_margin`        |         `0.8` | Horizontal spacing between support and model, in units of `edge_width`.
`support_xy_expansion`  |         `2.0` | Expand support map by this amount. Larger values will generate more support material, but the supports will be stronger.
`support_density`       |         `0.3` | Support structure density.
`support_flow_mult`     |        `0.75` | Flow rate is multiplied by this value for the support structure. Smaller values will generate a weaker support structure, but it will be easier to remove.
`min_layer_time`        |         `8.0` | Minimum layer time.
`layer_time_samples`    |           `5` | Number of samples in the layer time moving average.
`min_feed_rate`         |        `10.0` | Minimum feed rate.
`brim_width`            |         `0.0` | Brim width.
`material_density`      |     `0.00125` | Material density in `arbitrary_mass_unit / input_output_unit^3`. The default is approximately correct for PLA and millimeter input/output units.
`material_cost`         |     `0.01499` | Material cost in `arbitrary_currency / arbitrary_mass_unit`. The arbitrary mass unit must be the same as used in `material_density`.

#### Variables:

These variables are expanded to the value of the specified option in any G-code
string option (such as `start_gcode`). Prefix with a `%`.

Variable | Option
---------|-------
`t`      | `temp`
`b`      | `bed_temp`
`R`      | `retract_len`

Use `%%` for a literal `%`.

#### Configuration files:

Configuration files are a simple key-value format. Comment lines start with a
`#`. Lines starting with whitespace are interpreted as a continuation of the
previous line (except for comment lines).

Example:

	### Generic printer configuration
	extrusion_width=0.4
	perimeter_feed_rate=25
	loop_feed_rate=40
	infill_feed_rate=40
	travel_feed_rate=120
	first_layer_mult=0.5
	cool_on_gcode=M106 S255
	cool_off_gcode=M106 S0
	start_gcode=
		G21            ; Metric
		G90            ; Absolute positioning
		M82            ; Absolute extruder position
		G28 X0 Y0      ; Home X and Y
		G28 Z0         ; Home Z
		G1 Z20.0 F3300
		M190 S%b       ; Wait for bed to reach temp
		M109 S%t       ; Wait for hotend to reach temp
		G92 E0
		G1 E5.0 F200   ; Prime extruder
		G92 E0
		G1 E-%R F1800  ; Retract
		G92 E0
	end_gcode=
		M104 S0        ; Zero temps
		M140 S0
		G28 X0 Y0
		M84            ; Disable steppers
	xy_scale_factor=1.003
	z_scale_factor=1.0
	packing_density=0.99
	material_diameter=1.735
	retract_len=0.9
	retract_speed=20
	temp=220
	bed_temp=65
	cool_layer=2
	min_layer_time=10
	min_feed_rate=5

#### Examples:

Slice `infile.stl` and output gcode to `outfile.gcode`:

	shiv -o outfile.gcode infile.stl

Same as above, but with some options set:

	shiv -o outfile.gcode -d 0.5 -n 3 -x 30 -y 30 -O min_layer_time=15 -O temp=200 infile.stl

Preview slices of `infile.stl` in gnuplot:

	shiv -p infile.stl | gnuplot

### Feature wishlist:

* ASCII STL, AMF, OBJ, etc. support (only reads binary STL currently)
* Support structure generation (in progress)
* Bridge detection
* Multi-extrusion support
