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

	shiv [-hp] [-o output_path] [-c config_path] [-S setting=value]
	     [-l layer_height] [-w extrusion_width] [-t tolerance]
	     [-s scale_factor] [-d infill_density] [-n shells]
	     [-r roof_thickness] [-f floor_thickness] [-b brim_width]
	     [-C coarseness] [-x x_translate] [-y y_translate]
	     [-z z_chop] binary_stl_file

#### Flags:

Flag                 | Description
---------------------|-----------------------------------------
`-h`                 | Show help text
`-p`                 | Preview slices (pipe stdout to gnuplot)
`-o output_path`     | Output gcode path
`-c config_path`     | Configuration file path
`-S setting=value`   | Set setting to value
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
`-z z_chop`          | Sink object into build plate

#### Full settings list:

Setting                    | Default value | Description
---------------------------|--------------:|------------
`layer_height`             |         `0.2` | Layer height.
`tolerance`                |       `0.001` | Segment connection tolerance. Large values can be used to close holes in broken models.
`scale_constant`           |   `1000000.0` | Clipper uses integers, so we need to scale floating point values. Precision is `1/scale_constant` units. Coordinates in the range `±4.6e+18/scale_constant` are accepted.
`coarseness`               |        `0.01` | Approximate output coarseness. Useful for simplifying high polygon count meshes.
`extrusion_width`          |        `0.45` | Constrained extrusion width. Should *generally* be set to a value equal to or slightly larger than the nozzle diameter.
`xy_scale_factor`          |       `1.003` | The object is scaled by this ratio in the x and y axes to compensate for shrinkage. Around 1.003 works for PLA. ABS should be somewhere between 1.004 and 1.009.
`z_scale_factor`           |         `1.0` | The object is scaled by this ratio in the z axis to compensate for shrinkage. Should probably be left at 1 unless a high temperature heated chamber is used.
`x_center`                 |         `0.0` | X coordinate to center the object on.
`y_center`                 |         `0.0` | Y coordinate to center the object on.
`packing_density`          |        `0.98` | Solid packing density. Should be slightly less than 1. 0.98 seems to work for PLA.
`edge_packing_density`     |        `0.95` | Packing density of the constrained half of the outer perimeter (relative to `packing_density`).
`shell_clip`               |        `0.15` | Length to clip off the ends of shells in units of `extrusion_width`.
`extra_offset`             |         `0.0` | Offset the object by this distance in the xy plane.
`infill_density`           |         `0.2` | Sparse infill density.
`infill_pattern`           |        `grid` | Sparse infill pattern. Legal values are `grid`, `triangle`, `triangle2`, and `rectilinear`.
`solid_infill_angle`       |        `45.0` | Solid infill angle in degrees.
`sparse_infill_angle`      |        `45.0` | Sparse infill angle in degrees.
`shells`                   |           `2` | Number of loops/perimeters/shells (whatever you want to call them).
`roof_thickness`           |         `0.8` | Solid surface thickness when looking upwards.
`floor_thickness`          |         `0.8` | Solid surface thickness when looking downwards.
`min_shell_contact`        |         `1.0` | Minimum contact patch through roof and floor layers in units of `extrusion_width`. Small values (~0.1 to 1.0) will reduce print time compared to larger values, but may produce weaker parts.
`solid_fill_expansion`     |         `1.0` | Distance to expand solid infill, in units of `extrusion_width`.
`material_diameter`        |        `1.75` | Diameter of material.
`flow_multiplier`          |         `1.0` | Flow rate adjustment to compensate for incorrect E-steps or E-step variation between material types.
`feed_rate`                |        `50.0` | Base feed rate. Feed rates below are actual speeds (in units/s) if set to a positive value, or a multiple of `feed_rate` if set to a negative value. In other words, `40` is 40 units/s, but `-0.5` is `feed_rate * 0.5` units/s.
`perimeter_feed_rate`      |        `-0.5` | Outer shell feed rate.
`loop_feed_rate`           |        `-0.7` | Inner shell(s) feed rate.
`infill_feed_rate`         |        `None` | Sets both `solid_infill_feed_rate` and `sparse_infill_feed_rate`.
`solid_infill_feed_rate`   |        `-1.0` | Solid infill feed rate.
`sparse_infill_feed_rate`  |        `-1.0` | Sparse infill feed rate.
`support_feed_rate`        |        `-1.0` | Support structure feed rate.
`travel_feed_rate`         |       `120.0` | Travel feed rate.
`first_layer_mult`         |         `0.5` | First layer feed rates (except travel) are multiplied by this value.
`coast_len`                |         `0.0` | Length to coast (move with the extruder turned off) at the end of a shell. This can reduce start/end blobs if set correctly, but will cause gaps if set too high.
`wipe_len`                 |         `0.0` | Length to wipe the nozzle after printing a shell. A non-zero value will cause an unconditional retract at the end of each shell.
`retract_len`              |         `1.0` | Retraction length.
`retract_speed`            |        `20.0` | Retraction speed.
`moving_retract_speed`     |        `-0.5` | Retrect speed when doing non-stationary retracts. A negative value means a multiple of `retract_speed`. 
`restart_speed`            |        `-1.0` | Restart speed. A negative value means a multiple of `retract_speed`.
`retract_min_travel`       |         `5.0` | Minimum travel for retraction when not crossing a boundary or when printing shells. Has no effect when printing infill if `retract_within_island` is false.
`retract_threshold`        |        `30.0` | Unconditional retraction threshold.
`retract_within_island`    |       `false` | If false, retraction will not occur unless a boundary is crossed or the travel distance is greater than `retract_threshold`.
`retract_after_shells`     |        `true` | Retract unconditionally after printing the last shell.
`moving_retract`           |       `false` | Do a non-stationary retraction at the end of each shell.
`extra_restart_len`        |         `0.0` | Extra material length on restart.
`cool_layer`               |           `2` | Turn on part cooling at this layer (numbered from zero). Set to `-1` to disable cooling.
`start_gcode`              |        `None` | Prepend this G-code to beginning of the output file.
`end_gcode`                |        `None` | Append this G-code to the end of the output file.
`cool_on_gcode`            |   `M106 S255` | G-code to turn on cooling.
`cool_off_gcode`           |        `M107` | G-code to turn off cooling.
`edge_overlap`             |         `0.5` | Allowable edge path overlap in units of `extrusion_width`.
`comb`                     |        `true` | Avoid crossing boundaries.
`strict_shell_order`       |       `false` | Always do insets in order within an island.
`align_seams`              |        `true` | Align seams to the lower left corner. The nearest point is picked instead if this is false.
`align_interior_seams`     |        `true` | Align interior seams to the lower left corner if `align_seams` is also true. If false, only exterior seams are aligned.
`simplify_insets`          |        `true` | Do `rdp_simplify_path()` operation on all insets (only the initial outline is simplified if this is false)
`fill_inset_gaps`          |        `true` | Fill gaps between shells.
`no_solid`                 |       `false` | If true, only generate solid fill on the very top and bottom of the model.
`anchor`                   |       `false` | Clip and anchor inset paths.
`outside_first`            |       `false` | Prefer exterior shells.
`connect_solid_infill`     |       `false` | Connect the ends of solid infill lines together, forming a zig-zag instead of individual lines.
`solid_infill_first`       |        `true` | Print solid infill before sparse infill. Both infill types will be planned together if this is false. Will be set to `true` automatically if `solid_infill_feed_rate` and `sparse_infill_feed_rate` are not equal or if `connect_solid_infill` is true.
`separate_z_travel`        |       `false` | Generate a separate z travel move instead of moving all axes together.
`combine_all`              |       `false` | Orients all outlines counter-clockwise. This can be used to fix certain broken models, but it also fills holes.
`poly_fill_type`           |    `non_zero` | Poly fill type for union. Sometimes `even_odd` is useful for broken models with self-intersections and/or incorrect normals.
`inset_join_type`          |       `miter` | Join type for negative offsets. Legal values are `miter`, `square`, and `round`. `square` tends to retain tiny details better, but `miter` produces simpler (smaller) gcode.
`outset_join_type`         |       `miter` | Join type for positive offsets. Legal values are `miter`, `square`, and `round`.
`offset_miter_limit`       |         `2.0` | Sets `ClipperOffset.MiterLimit`. See the [ClipperLib documentation](http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Classes/ClipperOffset/Properties/MiterLimit.htm) for details.
`offset_arc_tolerance`     |         `5.0` | Sets `ClipperOffset.ArcTolerance`. See the [ClipperLib documentation](http://www.angusj.com/delphi/clipper/documentation/Docs/Units/ClipperLib/Classes/ClipperOffset/Properties/ArcTolerance.htm) for details.
`fill_threshold`           |        `0.25` | Infill and inset gap fill is removed when it would be narrower than `extrusion_width * fill_threshold`.
`min_sparse_infill_len`    |         `1.0` | Minimum length for sparse infill lines.
`connected_infill_overlap` |        `0.15` | Extra overlap between connected solid infill and shells in units of `extrusion_width`. Extruded volume does not change.
`generate_support`         |       `false` | Generate support structure.
`support_everywhere`       |        `true` | False means only touching build plate.
`solid_support_base`       |        `true` | Make supports solid at layer 0.
`connect_support_lines`    |       `false` | Connect support lines together. Makes the support structure more robust, but harder to remove.
`expand_support_interface` |       `false` | Expand support interface by the distance between support lines.
`support_angle`            |        `70.0` | Angle threshold for support.
`support_margin`           |         `0.6` | Horizontal spacing between support and model, in units of `edge_width`.
`support_vert_margin`      |           `1` | Vertical spacing between support and model, in layers.
`interface_layers`         |           `3` | Number of support interface layers.
`support_xy_expansion`     |         `2.0` | Expand support map by this amount. Larger values will generate more support material, but the supports will be stronger.
`support_density`          |         `0.2` | Support structure density.
`interface_density`        |         `0.7` | Support interface density.
`support_flow_mult`        |        `0.75` | Flow rate is multiplied by this value for the support structure. Smaller values will generate a weaker support structure, but it will be easier to remove. The default works well for PLA, but should be increased for materials that have trouble bridging (like PETG).
`support_wipe_len`         |         `5.0` | Wipe the nozzle over the previously printed line if a boundary will be crossed.
`min_layer_time`           |         `8.0` | Minimum layer time.
`layer_time_samples`       |           `5` | Number of samples in the layer time moving average.
`min_feed_rate`            |        `10.0` | Minimum feed rate.
`brim_width`               |         `0.0` | Brim width.
`brim_adhesion_factor`     |         `0.5` | How stuck to the object the brim is. 0 is just touching and 1 is packed as tightly as normal shells.
`raft_xy_expansion`        |         `5.0` | Expand raft beyond the model by this amount.
`raft_base_layer_height`   |         `0.3` | Layer height for the first layer of the raft.
`raft_base_layer_width`    |         `0.6` | Extrusion width for the first layer of the raft.
`raft_base_layer_density`  |         `0.5` | Fill density for the first layer of the raft.
`raft_vert_margin`         |         `1.0` | Vertical gap between the model and raft, in units of `layer_height`.
`raft_interface_flow_mult` |        `0.75` | Flow rate is multiplied by this value for the raft interface layers.
`raft_interface_layers`    |           `1` | Number of solid interface layers.
`material_density`         |     `0.00125` | Material density in `arbitrary_mass_unit / input_output_unit^3`. The default is approximately correct for PLA and millimeter input/output units.
`material_cost`            |     `0.01499` | Material cost in `arbitrary_currency / arbitrary_mass_unit`. The arbitrary mass unit must be the same as used in `material_density`.
`gcode_variable`           |        `None` | Set a variable that can be expanded within a G-code string option (see "G-code variables" below).
`at_layer`                 |        `None` | Print a string to the output file at the beginning of a given layer (numbered from zero).

#### G-code variables:

Variables and settings can be expanded within the `start_gcode`, `end_gcode`,
`cool_on_gcode`, `cool_off_gcode`, and `at_layer` settings by enclosing the
variable/setting name in curly braces: `{name}`. Alternates can be specified by
separating multiple names with colons: `{var1:var2:var3}`. If a
variable/setting does not exist, the next alternate will be tried. If none of
the alternates exist, it will expand to an empty string and warning will be
printed.

#### Configuration files:

Configuration files are a simple key-value format. Comment lines start with a
`#`. Lines starting with whitespace are interpreted as a continuation of the
previous line (except for comment lines).

Example:

	### Generic printer configuration
	x_center=100
	y_center=100
	extrusion_width=0.4
	feed_rate=60
	travel_feed_rate=120
	first_layer_mult=0.5
	cool_on_gcode=M106 S{fan_speed}
	cool_off_gcode=M106 S0
	start_gcode=
		G21              ; Metric
		G90              ; Absolute positioning
		M82              ; Absolute extruder position
		G28 X0 Y0        ; Home X and Y
		G28 Z0           ; Home Z
		G1 Z20.0 F3300
		M190 S{first_layer_bed_temp:bed_temp} ; Wait for bed to reach temp
		M109 S{first_layer_temp:temp}         ; Wait for hotend to reach temp
		G92 E0
		G1 E20.0 F100    ; Prime extruder
		G92 E0
		G1 E-{retract_len} F{retract_speed} ; Retract
		G92 E0
	end_gcode=
		M104 S0          ; Zero temps
		M140 S0
		G28 X0 Y0
		M84              ; Disable steppers
	at_layer=1=
		M104 S{temp}     ; Set second layer temp
		M140 S{bed_temp} ; Set second layer bed temp
	xy_scale_factor=1.003
	z_scale_factor=1.0
	packing_density=0.98
	edge_packing_density=0.94
	material_diameter=1.735
	retract_len=0.5
	retract_speed=20
	gcode_variable=first_layer_temp=220
	gcode_variable=first_layer_bed_temp=65
	gcode_variable=temp=210
	gcode_variable=bed_temp=60
	gcode_variable=fan_speed=255
	cool_layer=2
	min_layer_time=10
	min_feed_rate=5

#### Examples:

Slice `infile.stl` and output gcode to `outfile.gcode`:

	shiv -o outfile.gcode infile.stl

Same as above, but with some options set:

	shiv -o outfile.gcode -d 0.5 -n 3 -x 30 -y 30 -O min_layer_time=15 -O gcode_variable=temp=200 infile.stl

Preview slices of `infile.stl` in gnuplot:

	shiv -p infile.stl | gnuplot

### Feature wishlist:

* ASCII STL, AMF, OBJ, etc. support (only reads binary STL currently)
* Bridge detection
* Multi-extrusion support
