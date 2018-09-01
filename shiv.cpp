/*
 * Copyright (C) 2016-2018 Michael Barbour <barbour.michael.0@gmail.com>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * 
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <climits>
#include <limits>
#include <chrono>
#include <algorithm>
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "clipper.hpp"
#include "misc_defs.h"
#include "list.h"

#define SIMPLIFY_EPSILON (config.coarseness * config.scale_constant)
#define SHIV_DEBUG 1

#ifdef SHIV_FL_T_IS_FLOAT
#pragma message("fl_t defined as float")
typedef float fl_t;
#else
typedef double fl_t;
#endif

#ifdef SHIV_SINGLE_THREADED_PATH_PLANNING
#pragma message("using single-threaded path planner")
#endif

#define LENGTH(x) (sizeof(x) / sizeof(x[0]))
#define MINIMUM(a, b) (((a) < (b)) ? (a) : (b))
#define MAXIMUM(a, b) (((a) > (b)) ? (a) : (b))
#define MINIMUM_4(a, b, c, d) MINIMUM(MINIMUM((a), (b)), MINIMUM((c), (d)))
#define MAXIMUM_4(a, b, c, d) MAXIMUM(MAXIMUM((a), (b)), MAXIMUM((c), (d)))
#define FL_T_TO_CINT(x) ((ClipperLib::cInt) (llround((x) * config.scale_constant)))
#define CINT_TO_FL_T(x) (((fl_t) (x)) / config.scale_constant)
#define FL_T_TO_INTPOINT(x, y)  ClipperLib::IntPoint(FL_T_TO_CINT(x), FL_T_TO_CINT(y))
#define INTPOINT_TO_FL_T(point) CINT_TO_FL_T((point).X), CINT_TO_FL_T((point).Y)
#define FREE_VECTOR(v) decltype(v)().swap(v)  /* Force memory used by vector to be freed */

#if SHIV_DEBUG
#define DEBUG(...) fprintf(stderr, "DEBUG: " __VA_ARGS__)
#else
#define DEBUG(...) /* no-op */
#endif

static const char e_nomem[] = "fatal: No memory\n";
static const char usage_string[] =
	"usage: shiv [-hp] [-o output_path] [-c config_path] [-S setting=value]\n"
	"            [-l layer_height] [-w extrusion_width] [-t tolerance]\n"
	"            [-s scale_factor] [-d infill_density] [-n shells]\n"
	"            [-r roof_thickness] [-f floor_thickness] [-b brim_width]\n"
	"            [-C coarseness] [-x x_translate] [-y y_translate]\n"
	"            [-z z_chop] binary_stl_file\n"
	"\n"
	"flags:\n"
	"  -h                    show this help\n"
	"  -p                    preview slices (pipe stdout to gnuplot)\n"
	"  -o output_path        output gcode path\n"
	"  -c config_path        configuration file path\n"
	"  -S setting=value      set setting to value\n"
	"  -l layer_height       layer height\n"
	"  -w extrusion_width    constrained extrusion width\n"
	"  -t tolerance          segment connection tolerance\n"
	"  -s scale_factor       object scale factor\n"
	"  -d infill_density     sparse infill density\n"
	"  -n shells             number of shells\n"
	"  -r roof_thickness     solid roof thickness\n"
	"  -f floor_thickness    solid floor thickness\n"
	"  -b brim_width         brim width\n"
	"  -C coarseness         output coarseness\n"
	"  -x x_translate        translate object in the x-axis\n"
	"  -y y_translate        translate object in the y-axis\n"
	"  -z z_chop             sink object into build plate\n";

enum fill_pattern {
	FILL_PATTERN_GRID,
	FILL_PATTERN_TRIANGLE,
	FILL_PATTERN_TRIANGLE2,
	FILL_PATTERN_RECTILINEAR,
};

struct user_var {
	char *key, *value;
};

struct at_layer_gcode {
	int layer;
	char *value;
};

/* default config */
#define DEFAULT_COOL_ON_STR  "M106 S255"
#define DEFAULT_COOL_OFF_STR "M107"
static struct {
	fl_t layer_height             = 0.2;
	fl_t tolerance                = 0.001;      /* Segment connection tolerance */
	fl_t scale_constant           = 1000000.0;  /* Clipper uses integers, so we need to scale floating point values. Precision is 1/scale_constant units. Coordinates in the range `±4.6e+18/scale_constant` are accepted. */
	fl_t coarseness               = 0.01;       /* Approximate output coarseness. Useful for simplifying high polygon count meshes. */
	fl_t extrusion_width          = 0.45;       /* Constrained (solid infill) extrusion width */
	fl_t edge_width;                            /* Unconstrained (edge) extrusion width (calculated from extrusion_width) */
	fl_t extrusion_area;                        /* Cross-sectional area of an extrusion */
	fl_t xy_scale_factor          = 1.003;      /* Scale object in x and y axis by this ratio to compensate for shrinkage */
	fl_t z_scale_factor           = 1.0;        /* Scale object in z axis by this ratio to compensate for shrinkage */
	fl_t x_center                 = 0.0;
	fl_t y_center                 = 0.0;
	fl_t packing_density          = 0.75;       /* Controls how tightly packed the extrusions are. 0 means just touching and 1 means fully packed. */
	fl_t edge_packing_density     = 0.5;        /* Controls how much extra negative offset is applied to the outer perimeter to compensate for reduced packing density of its constrained edge. Automatically set to 1 when outside_first is true or shells is less than 2. */
	fl_t shell_clip               = 0.15;       /* Length to clip off the ends of shells in units of extrusion_width */
	fl_t extra_offset             = 0.0;        /* Offset the object by this distance in the xy plane */
	fl_t edge_offset;                           /* Offset of the outer perimeter (calculated) */
	fl_t infill_density           = 0.2;        /* Sparse infill density */
	fill_pattern infill_pattern   = FILL_PATTERN_GRID;  /* Sparse infill pattern */
	fl_t solid_infill_angle       = 45.0;       /* Solid infill angle (in degrees) */
	fl_t sparse_infill_angle      = 45.0;       /* Solid infill angle (in degrees) */
	int shells                    = 2;          /* Number of loops/perimeters/shells (whatever you want to call them) */
	fl_t roof_thickness           = 0.8;        /* Solid surface thickness when looking upwards */
	int roof_layers;                            /* Calculated from roof_thickness */
	fl_t floor_thickness          = 0.8;        /* Solid surface thickness when looking downwards */
	int floor_layers;                           /* Calculated from floor_thickness */
	fl_t min_shell_contact        = 1.0;        /* Minimum contact patch through roof and floor layers in units of extrusion_width. Small values (~0.1 to 1.0) will reduce print time compared to larger values, but may produce weaker parts. */
	fl_t solid_infill_clip_offset;
	fl_t solid_fill_expansion     = 1.0;        /* Distance to expand solid infill, in units of extrusion_width */
	fl_t material_diameter        = 1.75;       /* Diameter of material */
	fl_t material_area;                         /* Cross-sectional area of filament (calculated from material_diameter) */
	fl_t flow_multiplier          = 1.0;        /* Flow rate adjustment */
	fl_t feed_rate                = 50.0;       /* Base feed rate (all feed rates are in units/s) */
	/* Feed rates below are actual speeds if set to a positive value, or a multiple of 'feed_rate' (or other if specified) if set to a negative value.
	   In other words, '40' is 40 units/s, but '-0.5' is feed_rate * 0.5 units/s. */
	fl_t perimeter_feed_rate      = -0.5;       /* Outer shell feed rate */
	fl_t loop_feed_rate           = -0.7;       /* Inner shell feed rate */
	fl_t solid_infill_feed_rate   = -1.0;
	fl_t sparse_infill_feed_rate  = -1.0;
	fl_t support_feed_rate        = -1.0;
	fl_t iron_feed_rate           = -1.0;       /* Top surface ironing feed rate. A negative value means a multiple of 'solid_infill_feed_rate'. */
	fl_t travel_feed_rate         = 120.0;
	fl_t first_layer_mult         = 0.5;        /* First layer feed rates (except travel) are multiplied by this value */
	fl_t coast_len                = 0.0;        /* Length to coast (move with the extruder turned off) at the end of a shell */
	fl_t wipe_len                 = 0.0;        /* Length to wipe the nozzle after printing a shell. A non-zero value will cause an unconditional retract at the end of each shell. */
	fl_t retract_len              = 1.0;
	fl_t retract_speed            = 20.0;
	fl_t moving_retract_speed     = -0.5;       /* Retrect speed when doing non-stationary retracts. A negative value means a multiple of 'retract_speed'. */
	fl_t restart_speed            = -1.0;       /* A negative value means a multiple of 'retract_speed' */
	fl_t retract_min_travel       = 5.0;        /* Minimum travel for retraction when not crossing a boundary or when printing shells. Has no effect when printing infill if retract_within_island is false. */
	fl_t retract_threshold        = 30.0;       /* Unconditional retraction threshold */
	bool retract_within_island    = false;
	bool retract_after_shells     = false;      /* Retract unconditionally after printing the last shell */
	bool moving_retract           = false;      /* Do a non-stationary retraction at the end of each shell */
	fl_t extra_restart_len        = 0.0;        /* Extra material length on restart */
	int cool_layer                = 2;          /* Turn on part cooling at this layer */
	char *start_gcode             = NULL;
	char *end_gcode               = NULL;
	char *cool_on_gcode           = NULL;       /* Set in main() */
	char *cool_off_gcode          = NULL;       /* Set in main() */
	fl_t edge_overlap             = 0.5;        /* Allowable edge path overlap in units of extrusion_width */
	bool comb                     = true;       /* Avoid crossing boundaries */
	bool strict_shell_order       = false;      /* Always do insets in order within an island */
	bool align_seams              = true;       /* Align seams to the lower left corner */
	bool align_interior_seams     = true;       /* Align interior seams to the lower left corner if 'align_seams' is also true. If false, only exterior seams are aligned. */
	bool simplify_insets          = true;       /* Do simplify_path() operation on all insets (only the initial outline is simplified if this is false) */
	bool fill_inset_gaps          = true;       /* Fill gaps between shells */
	bool no_solid                 = false;      /* If true, only generate solid fill on the very top and bottom of the model */
	bool anchor                   = false;      /* Clip and anchor inset paths */
	bool outside_first            = false;      /* Prefer exterior shells */
	bool iron_top_surface         = false;      /* Run the nozzle over exposed top surfaces a second time */
	bool separate_z_travel        = false;      /* Generate a separate z travel move instead of moving all axes together */
	bool preserve_layer_offset    = false;      /* Preserve layer offset when placing the object on the build plate. Useful for certain multi-part prints. */
	bool combine_all              = false;      /* Orients all outlines counter-clockwise. This can be used to fix certain broken models, but it also fills holes. */
	ClipperLib::PolyFillType poly_fill_type = ClipperLib::pftNonZero;  /* Set poly fill type for union. Sometimes ClipperLib::pftEvenOdd is useful for broken models with self-intersections and/or incorrect normals. */
	ClipperLib::JoinType inset_join_type    = ClipperLib::jtMiter;     /* Join type for negative offsets */
	ClipperLib::JoinType outset_join_type   = ClipperLib::jtMiter;     /* Join type for positive offsets */
	fl_t offset_miter_limit       = 2.0;
	fl_t offset_arc_tolerance     = 5.0;
	fl_t fill_threshold           = 0.25;       /* Infill and inset gap fill is removed when it would be narrower than 'extrusion_width' * 'fill_threshold' */
	fl_t infill_smooth_threshold  = 2.0;        /* Solid infill lines are converted to a smooth curve when the region being filled is narrower than 'extrusion_width' * 'infill_smooth_threshold' */
	fl_t min_sparse_infill_len    = 1.0;        /* Minimum length for sparse infill lines */
	fl_t infill_overlap           = 0.05;       /* Overlap between infill and shells in units of 'extrusion_width' */
	fl_t iron_flow_multiplier     = 0.1;        /* Flow adjustment (relative to normal flow) for top surface ironing */
	fl_t iron_density             = 2.0;        /* Density of passes for top surface ironing */
	bool generate_support         = false;      /* Generate support structure */
	bool support_everywhere       = true;       /* False means only touching build plate */
	bool solid_support_base       = true;       /* Make supports solid at layer 0 */
	bool connect_support_lines    = false;      /* Connect support lines together. Makes the support structure more robust, but harder to remove. */
	bool expand_support_interface = true;       /* Expand support interface by the distance between support lines */
	fl_t support_angle            = 70.0;       /* Angle threshold for support */
	fl_t support_margin           = 0.6;        /* Horizontal spacing between support and model, in units of edge_width */
	int support_vert_margin       = 1;          /* Vertical spacing between support and model, in layers */
	int interface_roof_layers     = 3;          /* Number of support interface layers when looking upwards */
	int interface_floor_layers    = 0;          /* Number of support interface layers when looking downwards */
	fl_t support_xy_expansion     = 2.0;        /* Expand support map by this amount. Larger values will generate more support material, but the supports will be stronger. */
	fl_t support_density          = 0.2;        /* Support structure density */
	fl_t interface_density        = 0.7;        /* Support interface density */
	fl_t interface_clip_offset;
	fl_t support_flow_mult        = 0.75;       /* Flow rate is multiplied by this value for the support structure. Smaller values will generate a weaker support structure, but it will be easier to remove. The default works well for PLA, but should be increased for materials that have trouble bridging (like PETG). */
	fl_t support_wipe_len         = 5.0;        /* Wipe the nozzle over the previously printed line if a boundary will be crossed */
	fl_t min_layer_time           = 8.0;        /* Slow down if the estimated layer time is less than this value */
	int layer_time_samples        = 5;          /* Number of samples in the layer time moving average */
	fl_t min_feed_rate            = 10.0;
	fl_t brim_width               = 0.0;
	int brim_lines;
	fl_t brim_adhesion_factor     = 0.5;        /* How stuck to the object the brim is. 0 is just touching and 1 is packed as tightly as normal shells. */
	bool generate_raft            = false;
	fl_t raft_xy_expansion        = 5.0;
	fl_t raft_base_layer_height   = 0.3;
	fl_t raft_base_layer_width    = 0.6;
	fl_t raft_base_layer_density  = 0.5;
	fl_t raft_vert_margin         = 1.0;        /* Vertical gap between the model and raft, in units of layer_height. */
	fl_t raft_interface_flow_mult = 0.75;       /* Flow rate is multiplied by this value for the raft interface layers. */
	int raft_interface_layers     = 1;          /* Number of solid interface layers. */
	fl_t material_density         = 0.00125;    /* Material density in <arbitrary mass unit> / <input/output unit>^3. The default is correct for PLA and millimeter input/output units */
	fl_t material_cost            = 0.01499;    /* Material cost in <arbitrary currency> / <arbitrary mass unit>. The arbitrary mass unit must be the same as used in material_density */

	std::vector<struct user_var> user_vars;     /* User-set variables */
	std::vector<struct at_layer_gcode> at_layer;

	/* internal stuff */
	fl_t xy_extra                 = 0.0;        /* Extra xy size (brim, raft, extra_offset, support_xy_expansion, etc...). */
	fl_t object_z_extra           = 0.0;        /* Extra z offset to apply to everything but the raft on gcode export. */
} config;

enum setting_type {
	SETTING_TYPE_FL_T,
	SETTING_TYPE_INT,
	SETTING_TYPE_BOOL,
	SETTING_TYPE_FILL_PATTERN,
	SETTING_TYPE_JOIN_TYPE,
	SETTING_TYPE_POLY_FILL_TYPE,
	SETTING_TYPE_STR,
};

struct setting {
	const char *name;
	setting_type type;
	bool read_only, is_feed_rate;
	union {
		struct { int l, h; } i;
		struct { fl_t l, h; } f;
	} range;
	bool range_low_eq, range_high_eq;
	void *data;
};

#define FL_T_INF std::numeric_limits<fl_t>::infinity()

#define SETTING(name, type, read_only, is_feed_rate, range_low, range_high, range_low_eq, range_high_eq) \
	{ #name, type, read_only, is_feed_rate, range_low, range_high, range_low_eq, range_high_eq, (void *) &config.name }

static const struct setting settings[] = {
	/*      name                       type                      read_only is_feed_rate    range_low  range_high    low_eq high_eq */
	SETTING(layer_height,              SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(tolerance,                 SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(scale_constant,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(coarseness,                SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(extrusion_width,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(edge_width,                SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(extrusion_area,            SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(xy_scale_factor,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(z_scale_factor,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(x_center,                  SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(y_center,                  SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(packing_density,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(edge_packing_density,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(shell_clip,                SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(extra_offset,              SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(edge_offset,               SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(infill_density,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(infill_pattern,            SETTING_TYPE_FILL_PATTERN,   false, false, { .i = { 0,         0        } }, false, false),
	SETTING(solid_infill_angle,        SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(sparse_infill_angle,       SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(shells,                    SETTING_TYPE_INT,            false, false, { .i = { 0,         INT_MAX  } }, true,  true),
	SETTING(roof_thickness,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(roof_layers,               SETTING_TYPE_INT,            true,  false, { .i = { 0,         0        } }, false, false),
	SETTING(floor_thickness,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(floor_layers,              SETTING_TYPE_INT,            true,  false, { .i = { 0,         0        } }, false, false),
	SETTING(min_shell_contact,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(solid_infill_clip_offset,  SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(solid_fill_expansion,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(material_diameter,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(material_area,             SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(flow_multiplier,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(feed_rate,                 SETTING_TYPE_FL_T,           false, true,  { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(perimeter_feed_rate,       SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(loop_feed_rate,            SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	/* 'infill_feed_rate' is a special case */
	SETTING(solid_infill_feed_rate,    SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(sparse_infill_feed_rate,   SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(support_feed_rate,         SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(iron_feed_rate,            SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(travel_feed_rate,          SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(first_layer_mult,          SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(coast_len,                 SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(wipe_len,                  SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(retract_len,               SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(retract_speed,             SETTING_TYPE_FL_T,           false, true,  { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(moving_retract_speed,      SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(restart_speed,             SETTING_TYPE_FL_T,           false, true,  { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(retract_min_travel,        SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(retract_threshold,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(retract_within_island,     SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(retract_after_shells,      SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(moving_retract,            SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(extra_restart_len,         SETTING_TYPE_FL_T,           false, false, { .f = { -FL_T_INF, FL_T_INF } }, false, false),
	SETTING(cool_layer,                SETTING_TYPE_INT,            false, false, { .i = { -1,        INT_MAX  } }, true,  true),
	SETTING(start_gcode,               SETTING_TYPE_STR,            false, false, { .i = { 0,         0        } }, false, false),
	SETTING(end_gcode,                 SETTING_TYPE_STR,            false, false, { .i = { 0,         0        } }, false, false),
	SETTING(cool_on_gcode,             SETTING_TYPE_STR,            false, false, { .i = { 0,         0        } }, false, false),
	SETTING(cool_off_gcode,            SETTING_TYPE_STR,            false, false, { .i = { 0,         0        } }, false, false),
	SETTING(edge_overlap,              SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(comb,                      SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(strict_shell_order,        SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(align_seams,               SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(align_interior_seams,      SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(simplify_insets,           SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(fill_inset_gaps,           SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(no_solid,                  SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(anchor,                    SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(outside_first,             SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(iron_top_surface,          SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(separate_z_travel,         SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(preserve_layer_offset,     SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(combine_all,               SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(poly_fill_type,            SETTING_TYPE_POLY_FILL_TYPE, false, false, { .i = { 0,         0        } }, false, false),
	SETTING(inset_join_type,           SETTING_TYPE_JOIN_TYPE,      false, false, { .i = { 0,         0        } }, false, false),
	SETTING(outset_join_type,          SETTING_TYPE_JOIN_TYPE,      false, false, { .i = { 0,         0        } }, false, false),
	SETTING(offset_miter_limit,        SETTING_TYPE_FL_T,           false, false, { .f = { 2.0,       FL_T_INF } }, true,  false),
	SETTING(offset_arc_tolerance,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.25,      FL_T_INF } }, true,  false),
	SETTING(fill_threshold,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(infill_smooth_threshold,   SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       4.0      } }, true,  true),
	SETTING(min_sparse_infill_len,     SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(infill_overlap,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       0.5      } }, true,  true),
	SETTING(iron_flow_multiplier,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(iron_density,              SETTING_TYPE_FL_T,           false, false, { .f = { 1.0,       FL_T_INF } }, true,  false),
	SETTING(generate_support,          SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(support_everywhere,        SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(solid_support_base,        SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(connect_support_lines,     SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(expand_support_interface,  SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(support_angle,             SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       90.0     } }, false, false),
	SETTING(support_margin,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),  /* FIXME: will cause problems with the combing code if set to 0 */
	SETTING(support_vert_margin,       SETTING_TYPE_INT,            false, false, { .i = { 0,         INT_MAX  } }, true,  true),
	SETTING(interface_roof_layers,     SETTING_TYPE_INT,            false, false, { .i = { 0,         INT_MAX  } }, true,  true),
	SETTING(interface_floor_layers,    SETTING_TYPE_INT,            false, false, { .i = { 0,         INT_MAX  } }, true,  true),
	SETTING(support_xy_expansion,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(support_density,           SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, false, true),
	SETTING(interface_density,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, false, true),
	SETTING(interface_clip_offset,     SETTING_TYPE_FL_T,           true,  false, { .f = { 0.0,       0.0      } }, false, false),
	SETTING(support_flow_mult,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, false, true),
	SETTING(support_wipe_len,          SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(min_layer_time,            SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(layer_time_samples,        SETTING_TYPE_INT,            false, false, { .i = { 1,         INT_MAX  } }, true,  true),
	SETTING(min_feed_rate,             SETTING_TYPE_FL_T,           false, true,  { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(brim_width,                SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(brim_lines,                SETTING_TYPE_INT,            true,  false, { .i = { 0,         0        } }, false, false),
	SETTING(brim_adhesion_factor,      SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, true,  true),
	SETTING(generate_raft,             SETTING_TYPE_BOOL,           false, false, { .i = { 0,         0        } }, false, false),
	SETTING(raft_xy_expansion,         SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(raft_base_layer_height,    SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(raft_base_layer_width,     SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(raft_base_layer_density,   SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       1.0      } }, false, true),
	SETTING(raft_vert_margin,          SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, true,  false),
	SETTING(raft_interface_flow_mult,  SETTING_TYPE_FL_T,           false, false, { .f = { 0.0,       FL_T_INF } }, false, false),
	SETTING(raft_interface_layers,     SETTING_TYPE_INT,            false, false, { .i = { 0,         INT_MAX  } }, true,  true),
	SETTING(material_density,          SETTING_TYPE_FL_T,           false, false, { .f = { 0,         FL_T_INF } }, true,  false),
	SETTING(material_cost,             SETTING_TYPE_FL_T,           false, false, { .f = { 0,         FL_T_INF } }, true,  false),
};

struct vertex {
	fl_t x, y, z;
};

struct triangle {
	struct vertex v[3];
};

struct object {
	ssize_t n, n_slices;
	struct vertex c;
	fl_t h, w, d;
	struct triangle *t;
	struct slice *slices;

	ClipperLib::Paths solid_infill_patterns[2];
	std::vector<ClipperLib::Paths> brim;
	ClipperLib::Paths raft[2];
	ClipperLib::Paths support_pattern;
	ClipperLib::Paths support_interface_pattern;
	ClipperLib::Paths raft_base_layer_pattern;
};

struct segment_list {
	struct segment *head, *tail;
};

struct segment {
	struct segment *next, *prev;
	fl_t x[2], y[2];
};

struct cint_rect {
	ClipperLib::cInt x0, y0, x1, y1;
};

struct island {
	ClipperLib::Paths *insets;
	ClipperLib::Paths *inset_gaps;
	ClipperLib::Paths infill_insets;
	ClipperLib::Paths solid_infill;
	ClipperLib::Paths sparse_infill;
	ClipperLib::Paths boundaries;        /* Boundaries when moving inside (slightly outset from the first inset) */
	ClipperLib::Paths comb_paths;        /* Equal to the first inset (needed because insets get erased during planning) */
	ClipperLib::Paths outer_boundaries;  /* Boundaries when moving outside */
	ClipperLib::Paths outer_comb_paths;  /* Slightly outset from outer_boundaries */
	ClipperLib::Paths solid_infill_clip;
	ClipperLib::Paths solid_infill_boundaries;
	ClipperLib::Paths exposed_surface;
	ClipperLib::Paths constraining_edge; /* Slightly inset from infill_insets. Used to determine whether infill lines should be connected. */
	ClipperLib::Paths iron_paths;        /* Paths to follow for top surface ironing */
	struct cint_rect box;  /* bounding box */
};

struct g_move {
	ClipperLib::cInt x, y, z;
	fl_t e, feed_rate;
	bool scalable, is_travel, is_restart;
};

struct slice {
#ifdef _OPENMP
	omp_lock_t lock;
	volatile ssize_t n_seg, s_len;
#else
	ssize_t n_seg, s_len;
#endif
	struct segment *s;
	std::vector<struct island> islands;
	std::vector<struct g_move> moves;
	ClipperLib::PolyTree layer_support_map;
	ClipperLib::Paths support_map;
	ClipperLib::Paths support_boundaries;
	ClipperLib::Paths support_interface_clip;
	ClipperLib::Paths support_lines;
	ClipperLib::Paths support_interface_lines;
	ClipperLib::Paths last_boundaries, last_comb_paths;
	ClipperLib::Paths printed_outer_boundaries, printed_outer_comb_paths;
	fl_t layer_time;
};

struct machine {
	ClipperLib::cInt x, y, z;
	fl_t e, feed_rate;
	bool is_retracted, force_retract;
};

static void die(const char *s, int r)
{
	fputs(s, stderr);
	exit(r);
}

static char * isolate(char *s, char c)
{
	while (*s && *s != c) ++s;
	if (*s != '\0') *s++ = '\0';
	return s;
}

/* We could get the file size with fseek/ftell and read the whole file in one
   go, but that doesn't work with pipes and some other special files... */
static char * get_file_contents(const char *path)
{
	const size_t g = 512;
	size_t s = g, p = 0;
	FILE *f;
	char *c;

	f = fopen(path, "r");
	if (!f)
		return NULL;
	c = (char *) calloc(s, sizeof(char));
	if (!c)
		die(e_nomem, 2);
	while (!ferror(f)) {
		p += fread(&c[p], 1, s - p, f);
		if (p >= s) {
			s += g;
			c = (char *) realloc(c, s * sizeof(char));
			if (!c)
				die(e_nomem, 2);
		}
		if (feof(f)) {
			c[p] = '\0';
			fclose(f);
			return c;
		}
	}
	free(c);
	fclose(f);
	return NULL;
}

static const struct setting * find_config_setting(const char *key)
{
	for (size_t i = 0; i < LENGTH(settings); ++i)
		if (strcmp(key, settings[i].name) == 0)
			return &settings[i];
	return NULL;
}

static int set_config_setting(const char *key, const char *value, int n, const char *path)
{
	char *endptr;
	const struct setting *s = find_config_setting(key);
	if (s) {
		if (s->read_only) {
			if (path) fprintf(stderr, "error: line %d in %s: setting %s is read-only\n", n, path, s->name);
			else      fprintf(stderr, "error: setting %s is read-only\n", s->name);
			return 1;
		}
		errno = 0;
		switch(s->type) {
		case SETTING_TYPE_FL_T:
		case SETTING_TYPE_INT:
			if (s->type == SETTING_TYPE_FL_T)
				*((fl_t *) s->data) = (fl_t) strtod(value, &endptr);
			else
				*((int *) s->data) = strtol(value, &endptr, 10);
			if (errno != 0) {
				if (path) fprintf(stderr, "error: line %d in %s: %s: %s\n", n, path, strerror(errno), value);
				else      fprintf(stderr, "error: %s: %s\n", strerror(errno), value);
				return 1;
			}
			else if (endptr == value) {
				if (path) fprintf(stderr, "error: line %d in %s: invalid input: %s\n", n, path, value);
				else      fprintf(stderr, "error: invalid input: %s\n", value);
				return 1;
			}
			else if (*endptr != '\0') {
				if (path) fprintf(stderr, "warning: line %d in %s: trailing characters: %s\n", n, path, endptr);
				else      fprintf(stderr, "warning: trailing characters: %s\n", endptr);
			}
			if (s->type == SETTING_TYPE_FL_T) {
				if (!(((s->range_low_eq) ? *((fl_t *) s->data) >= s->range.f.l : *((fl_t *) s->data) > s->range.f.l)
						&& ((s->range_high_eq) ? *((fl_t *) s->data) <= s->range.f.h : *((fl_t *) s->data) < s->range.f.h))) {
					if (path) fprintf(stderr, "error: line %d in %s: %s must be ", n, path, s->name);
					else      fprintf(stderr, "error: %s must be ", s->name);
					if (s->range.f.l > -FL_T_INF && s->range.f.h < FL_T_INF)
						fprintf(stderr, "within %c%g,%g%c\n", (s->range_low_eq) ? '[' : '(', s->range.f.l, s->range.f.h, (s->range_high_eq) ? ']' : ')');
					else if (s->range.f.l == -FL_T_INF)
						fprintf(stderr, "%s %g\n", (s->range_high_eq) ? "<=" : "<", s->range.f.h);
					else  /* s->range.f.h == FL_T_INF */
						fprintf(stderr, "%s %g\n", (s->range_low_eq) ? ">=" : ">", s->range.f.l);
					return 1;
				}
			}
			else {
				if (!(((s->range_low_eq) ? *((int *) s->data) >= s->range.i.l : *((int *) s->data) > s->range.i.l)
						&& ((s->range_high_eq) ? *((int *) s->data) <= s->range.i.h : *((int *) s->data) < s->range.i.h))) {
					if (path) fprintf(stderr, "error: line %d in %s: %s must be ", n, path, s->name);
					else      fprintf(stderr, "error: %s must be ", s->name);
					if (s->range.i.l > INT_MIN && s->range.i.h < INT_MAX)
						fprintf(stderr, "within %c%d,%d%c\n", (s->range_low_eq) ? '[' : '(', s->range.i.l, s->range.i.h, (s->range_high_eq) ? ']' : ')');
					else if (s->range.i.l == INT_MIN)
						fprintf(stderr, "%s %d\n", (s->range_high_eq) ? "<=" : "<", s->range.i.h);
					else  /* s->range.i.h == INT_MAX */
						fprintf(stderr, "%s %d\n", (s->range_low_eq) ? ">=" : ">", s->range.i.l);
					return 1;
				}
			}
			break;
		case SETTING_TYPE_BOOL:
			*((bool *) s->data) = (value[0] == 't' || value[0] == 'T' || value[0] == 'y' || value[0] == 'Y' || atoi(value));
			break;
		case SETTING_TYPE_FILL_PATTERN:
			if (strcmp(value, "grid") == 0)
				*((fill_pattern *) s->data) = FILL_PATTERN_GRID;
			else if (strcmp(value, "triangle") == 0)
				*((fill_pattern *) s->data) = FILL_PATTERN_TRIANGLE;
			else if (strcmp(value, "triangle2") == 0)
				*((fill_pattern *) s->data) = FILL_PATTERN_TRIANGLE2;
			else if (strcmp(value, "rectilinear") == 0)
				*((fill_pattern *) s->data) = FILL_PATTERN_RECTILINEAR;
			else {
				if (path) fprintf(stderr, "error: line %d in %s: illegal value for %s: %s\n", n, path, s->name, value);
				else      fprintf(stderr, "error: illegal value for %s: %s\n", s->name, value);
				return 1;
			}
			break;
		case SETTING_TYPE_JOIN_TYPE:
			if (strcmp(value, "miter") == 0)
				*((ClipperLib::JoinType *) s->data) = ClipperLib::jtMiter;
			else if (strcmp(value, "square") == 0)
				*((ClipperLib::JoinType *) s->data) = ClipperLib::jtSquare;
			else if (strcmp(value, "round") == 0)
				*((ClipperLib::JoinType *) s->data) = ClipperLib::jtRound;
			else {
				if (path) fprintf(stderr, "error: line %d in %s: illegal value for %s: %s\n", n, path, s->name, value);
				else      fprintf(stderr, "error: illegal value for %s: %s\n", s->name, value);
				return 1;
			}
			break;
		case SETTING_TYPE_POLY_FILL_TYPE:
			if (strcmp(value, "even_odd") == 0)
				*((ClipperLib::PolyFillType *) s->data) = ClipperLib::pftEvenOdd;
			else if (strcmp(value, "non_zero") == 0)
				*((ClipperLib::PolyFillType *) s->data) = ClipperLib::pftNonZero;
			else if (strcmp(value, "positive") == 0)
				*((ClipperLib::PolyFillType *) s->data) = ClipperLib::pftPositive;
			else if (strcmp(value, "negative") == 0)
				*((ClipperLib::PolyFillType *) s->data) = ClipperLib::pftNegative;
			else {
				if (path) fprintf(stderr, "error: line %d in %s: illegal value for %s: %s\n", n, path, s->name, value);
				else      fprintf(stderr, "error: illegal value for %s: %s\n", s->name, value);
				return 1;
			}
			break;
		case SETTING_TYPE_STR:
			free(*((char **) s->data));
			*((char **) s->data) = strdup(value);
			break;
		default:
			fprintf(stderr, "BUG: unhandled setting type for '%s'\n", s->name);
			return 1;
		}
	}
	else if (strcmp(key, "infill_feed_rate") == 0) {
		if (set_config_setting("solid_infill_feed_rate", value, n, path)
				|| set_config_setting("sparse_infill_feed_rate", value, n, path))
			return 1;
	}
	else if (strcmp(key, "gcode_variable") == 0 || strcmp(key, "v") == 0) {
		char *s = strdup(value);
		struct user_var uv;
		uv.value = isolate(s, '=');
		uv.key = s;
		if (find_config_setting(uv.key)) {
			if (path) fprintf(stderr, "error: line %d in %s: cannot set variable %s: is a setting\n", n, path, uv.key);
			else      fprintf(stderr, "error: cannot set variable %s: is a setting\n", uv.key);
			free(s);
			return 1;
		}
		for (auto it = config.user_vars.begin(); it != config.user_vars.end(); ++it) {
			if (strcmp(it->key, uv.key) == 0) {
				free(it->key);
				config.user_vars.erase(it);
				break;
			}
		}
		config.user_vars.push_back(uv);
	}
	else if (strcmp(key, "at_layer") == 0) {
		char *s = strdup(value);
		struct at_layer_gcode g;
		g.value = isolate(s, '=');
		g.layer = atoi(s);
		config.at_layer.push_back(g);
	}
	else {
		if (path) fprintf(stderr, "error: line %d in %s: invalid setting: %s\n", n, path, key);
		else      fprintf(stderr, "error: invalid setting: %s\n", key);
		return 1;
	}
	return 0;
}

static void print_config_setting(FILE *f, const struct setting *s, bool convert_feed_rates)
{
	switch (s->type) {
	case SETTING_TYPE_FL_T:
		if (convert_feed_rates && s->is_feed_rate) {
			long int feed_rate = lround((double) *((fl_t *) s->data) * 60.0);
			fprintf(f, "%ld", (feed_rate < 1) ? 1 : feed_rate);
		}
		else
			fprintf(f, "%f", (double) *((fl_t *) s->data));
		break;
	case SETTING_TYPE_INT:
		fprintf(f, "%d", *((int *) s->data));
		break;
	case SETTING_TYPE_BOOL:
		fprintf(f, "%s", (*((bool *) s->data)) ? "true" : "false");
		break;
	case SETTING_TYPE_FILL_PATTERN:
		switch (*((fill_pattern *) s->data)) {
		case FILL_PATTERN_GRID:        fputs("grid", f);        break;
		case FILL_PATTERN_TRIANGLE:    fputs("triangle", f);    break;
		case FILL_PATTERN_TRIANGLE2:   fputs("triangle2", f);    break;
		case FILL_PATTERN_RECTILINEAR: fputs("rectilinear", f); break;
		}
		break;
	case SETTING_TYPE_JOIN_TYPE:
		switch (*((ClipperLib::JoinType *) s->data)) {
		case ClipperLib::jtMiter:  fputs("miter", f);  break;
		case ClipperLib::jtSquare: fputs("square", f); break;
		case ClipperLib::jtRound:  fputs("round", f);  break;
		}
		break;
	case SETTING_TYPE_POLY_FILL_TYPE:
		switch (*((ClipperLib::PolyFillType *) s->data)) {
		case ClipperLib::pftEvenOdd:  fputs("even_odd", f); break;
		case ClipperLib::pftNonZero:  fputs("non_zero", f); break;
		case ClipperLib::pftPositive: fputs("positive", f); break;
		case ClipperLib::pftNegative: fputs("negative", f); break;
		}
		break;
	case SETTING_TYPE_STR:
		fputs(*((char **) s->data), f);
		break;
	default:
		fprintf(stderr, "BUG: unhandled setting type for '%s'\n", s->name);
	}
}

static void write_gcode_string(const char *s, FILE *f, bool is_user_var)
{
	bool line_start = true;
	if (!s || strlen(s) == 0)
		return;
	while (*s) {
		if (!is_user_var && *s == '{') {
			const char *end_brace = strchr(s, '}');
			if (end_brace && *(++s) != '\0') {
				char *key_str = strndup(s, end_brace - s), *key = key_str;
				while (*key != '\0') {
					char *next_key = isolate(key, ':');
					const struct setting *setting = find_config_setting(key);
					if (setting) {
						print_config_setting(f, setting, true);
						goto found_var;
					}
					else {
						for (auto it = config.user_vars.begin(); it != config.user_vars.end(); ++it) {
							if (strcmp(key, it->key) == 0) {
								write_gcode_string(it->value, f, true);
								goto found_var;
							}
						}
						if (*next_key == '\0')
							fprintf(stderr, "warning: variable not found: %s\n", key);
					}
					key = next_key;
				}
				found_var:
				free(key_str);
				s = end_brace;
			}
			else {
				fprintf(stderr, "error: syntax: expected '}'\n");
				fputc('{', f);
			}
			line_start = false;
		}
		else if (*s == '\n') {
			if (!line_start)
				fputc('\n', f);
			line_start = true;
		}
		else if (!line_start || (line_start && *s != '\t' && *s != ' ')) {
			line_start = false;
			fputc(*s, f);
		}
		++s;
	}
	if (!is_user_var)
		fputc('\n', f);
}

static int read_config(const char *path)
{
	char *c, *key, *value, *next;
	c = get_file_contents(path);
	if (!c)
		return 1;
	key = c;
	for (ssize_t i = 1; key && *key != '\0'; ++i) {
		next = strchr(key, '\n');
		if (next) {
			while (next && (next[1] == ' ' || next[1] == '\t'))
				next = strchr(next + 1, '\n');
			if (next)
				next = isolate(next, '\n');
		}
		if (*key != '\0' && *key != '#') {
			value = isolate(key, '=');
			if (set_config_setting(key, value, i, path)) {
				free(c);
				return 2;
			}
		}
		key = next;
	}
	free(c);
	return 0;
}

static int read_binary_stl(struct object *o, const char *path)
{
	ssize_t i;
	int k, p;
	bool first = true;
	fl_t top = 0, bottom = 0, left = 0, right = 0, front = 0, back = 0;
	char header[80];
	float stl_poly[13];
	FILE *f;

	if (strcmp(path, "-") == 0)
		f = stdin;
	else
		f = fopen(path, "rb");
	if (!f)
		return 1;
	if (fread(header, 1, 80, f) != 80)
		goto read_fail;
	o->n = 0;
	if (fread(&o->n, 4, 1, f) != 1)
		goto read_fail;
	o->t = (struct triangle *) calloc(o->n, sizeof(struct triangle));
	if (!o->t)
		die(e_nomem, 2);
	for (i = 0; i < o->n; ++i) {
		if (fread(stl_poly, 50, 1, f) != 1)  /* Each polygon entry is 50 bytes */
			goto read_fail;
		p = 3;
		for (k = 0; k < 3; ++k) {
			o->t[i].v[k].x = (fl_t) stl_poly[p++];
			o->t[i].v[k].y = (fl_t) stl_poly[p++];
			o->t[i].v[k].z = (fl_t) stl_poly[p++];
			if (first || o->t[i].v[k].z > top)     top = o->t[i].v[k].z;
			if (first || o->t[i].v[k].z < bottom)  bottom = o->t[i].v[k].z;
			if (first || o->t[i].v[k].x > right)   right = o->t[i].v[k].x;
			if (first || o->t[i].v[k].x < left)    left = o->t[i].v[k].x;
			if (first || o->t[i].v[k].y > back)    back = o->t[i].v[k].y;
			if (first || o->t[i].v[k].y < front)   front = o->t[i].v[k].y;
			first = false;
		}
	}
	fclose(f);
	o->h = top - bottom;
	o->w = right - left;
	o->d = back - front;
	o->c.x = right - o->w / 2.0;
	o->c.y = back - o->d / 2.0;
	o->c.z = top - o->h / 2.0;
	return 0;

	read_fail:
	fclose(f);
	free(o->t);
	return 2;
}

static void translate_object(struct object *o, fl_t x, fl_t y, fl_t z)
{
	if (x == 0.0 && y == 0.0 && z == 0.0)
		return;
	for (ssize_t i = 0; i < o->n; ++i) {
		for (int k = 0; k < 3; ++k) {
			o->t[i].v[k].x += x;
			o->t[i].v[k].y += y;
			o->t[i].v[k].z += z;
		}
	}
	o->c.x += x;
	o->c.y += y;
	o->c.z += z;
}

static void scale_object(struct object *o, fl_t x_ratio, fl_t y_ratio, fl_t z_ratio)
{
	if (x_ratio == 1.0 && y_ratio == 1.0 && z_ratio == 1.0)
		return;
	for (ssize_t i = 0; i < o->n; ++i) {
		for (int k = 0; k < 3; ++k) {
			o->t[i].v[k].x *= x_ratio;
			o->t[i].v[k].y *= y_ratio;
			o->t[i].v[k].z *= z_ratio;
		}
	}
	o->h *= z_ratio;
	o->w *= x_ratio;
	o->d *= y_ratio;
	o->c.x *= x_ratio;
	o->c.y *= y_ratio;
	o->c.z *= z_ratio;
}

/* Project2D code ported from CuraEngine slicer.h (Copyright (C) 2013 David Braam) */
static void project2d(struct segment *seg, const struct vertex *v0, const struct vertex *v1, const struct vertex *v2, fl_t z)
{
	seg->x[0] = v0->x + (v1->x - v0->x) * (z - v0->z) / (v1->z - v0->z);
	seg->y[0] = v0->y + (v1->y - v0->y) * (z - v0->z) / (v1->z - v0->z);
	seg->x[1] = v0->x + (v2->x - v0->x) * (z - v0->z) / (v2->z - v0->z);
	seg->y[1] = v0->y + (v2->y - v0->y) * (z - v0->z) / (v2->z - v0->z);
}

static void find_segments(struct slice *slices, const struct triangle *t)
{
	fl_t min_z, max_z, z;
	ssize_t start, end, i;
	struct segment *s;

	max_z = MAXIMUM(t->v[0].z, t->v[1].z);
	max_z = MAXIMUM(max_z, t->v[2].z);
	max_z = MAXIMUM(0, max_z);  /* Chop off negative z values */

	min_z = MINIMUM(t->v[0].z, t->v[1].z);
	min_z = MINIMUM(min_z, t->v[2].z);
	min_z = MAXIMUM(0, min_z);  /* Chop off negative z values */

	start = (ssize_t) floor(min_z / config.layer_height + (0.4999));
	end = (ssize_t) floor(max_z / config.layer_height + (0.5001));

	for (i = start; i < end; ++i) {
		z = ((fl_t) i) * config.layer_height + config.layer_height / 2;
	#ifdef _OPENMP
		omp_set_lock(&slices[i].lock);
	#endif
		if (slices[i].n_seg >= slices[i].s_len) {
			if (slices[i].s_len == 0) {
				slices[i].s_len = 1024;
				slices[i].s = (struct segment *) malloc(sizeof(struct segment) * slices[i].s_len);
				if (!slices[i].s)
					die(e_nomem, 2);
			}
			else {
				slices[i].s_len = slices[i].s_len << 1;
				slices[i].s = (struct segment *) realloc(slices[i].s, sizeof(struct segment) * slices[i].s_len);
				if (!slices[i].s)
					die(e_nomem, 2);
			}
		}
		s = &slices[i].s[slices[i].n_seg];

		/* Intersection code ported from CuraEngine slicer.cpp (Copyright (C) 2013 David Braam) */
		if (t->v[0].z < z && t->v[1].z >= z && t->v[2].z >= z)
			project2d(s, &t->v[0], &t->v[2], &t->v[1], z);
		else if (t->v[0].z > z && t->v[1].z < z && t->v[2].z < z)
			project2d(s, &t->v[0], &t->v[1], &t->v[2], z);

		else if (t->v[1].z < z && t->v[0].z >= z && t->v[2].z >= z)
			project2d(s, &t->v[1], &t->v[0], &t->v[2], z);
		else if (t->v[1].z > z && t->v[0].z < z && t->v[2].z < z)
			project2d(s, &t->v[1], &t->v[2], &t->v[0], z);

		else if (t->v[2].z < z && t->v[1].z >= z && t->v[0].z >= z)
			project2d(s, &t->v[2], &t->v[1], &t->v[0], z);
		else if (t->v[2].z > z && t->v[1].z < z && t->v[0].z < z)
			project2d(s, &t->v[2], &t->v[0], &t->v[1], z);
		else
			goto invalid_segment;
		if (s->x[0] != s->x[1] || s->y[0] != s->y[1]) /* Ignore zero-length segments */
			++slices[i].n_seg;
		invalid_segment:
	#ifdef _OPENMP
		omp_unset_lock(&slices[i].lock);
	#else
		;
	#endif
	}
}

static void generate_islands(struct slice *slice, const ClipperLib::PolyNode *n)
{
	for (const ClipperLib::PolyNode *c : n->Childs) {
		struct island island = {};
		island.insets = new ClipperLib::Paths[(config.shells > 1) ? config.shells : 1]();
		if (!island.insets)
			die(e_nomem, 2);
		island.insets[0].push_back(c->Contour);
		for (const ClipperLib::PolyNode *cc : c->Childs) {
			island.insets[0].push_back(cc->Contour);
			generate_islands(slice, cc);
		}
		slice->islands.push_back(island);
	}
}

static void find_bounding_box(struct island *island)
{
	bool first = true;
	for (const ClipperLib::Path &path : island->insets[0]) {
		for (const ClipperLib::IntPoint &p : path) {
			if (first) {
				island->box.x0 = p.X;
				island->box.y0 = p.Y;
				island->box.x1 = p.X;
				island->box.y1 = p.Y;
				first = false;
			}
			else {
				if (p.X < island->box.x0)
					island->box.x0 = p.X;
				if (p.X > island->box.x1)
					island->box.x1 = p.X;
				if (p.Y > island->box.y0)
					island->box.y0 = p.Y;
				if (p.Y < island->box.y1)
					island->box.y1 = p.Y;
			}
		}
	}
}

static fl_t distance_to_point(const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1)
{
	const fl_t dx = p1.X - p0.X, dy = p1.Y - p0.Y;
	return sqrt(dx * dx + dy * dy);
}

static fl_t distance_to_line(const ClipperLib::IntPoint &p, const ClipperLib::IntPoint &l0, const ClipperLib::IntPoint &l1)
{
	const fl_t dx = l1.X - l0.X, dy = l1.Y - l0.Y;
	const fl_t len = dx * dx + dy * dy;
	if (len == 0.0)
		return distance_to_point(p, l0);
	/* Project p onto the line parameterized as l0 + t(l1 - l0) */
	const fl_t t = ((p.X - l0.X) * dx + (p.Y - l0.Y) * dy) / len;
	if (t < 0.0)  /* Falls beyond the first point */
		return distance_to_point(p, l0);
	else if (t > 1.0) /* Falls beyond the second point */
		return distance_to_point(p, l1);
	const ClipperLib::IntPoint proj(l0.X + t * dx, l0.Y + t * dy);
	return distance_to_point(p, proj);
}

static fl_t perpendicular_distance_to_line(const ClipperLib::IntPoint &p, const ClipperLib::IntPoint &l0, const ClipperLib::IntPoint &l1)
{
	const fl_t dx = l1.X - l0.X, dy = l1.Y - l0.Y;
	const fl_t len = dx * dx + dy * dy;
	if (len == 0.0)
		return distance_to_point(p, l0);
	const fl_t n = dx * (fl_t) (l0.Y - p.Y) - (fl_t) (l0.X - p.X) * dy;
	return fabs(n) / sqrt(len);
}

/* Note: epsilon is in scaled units */
static ClipperLib::Path rdp_simplify_path(const ClipperLib::Path &p, const fl_t epsilon)
{
	ClipperLib::Path res;
	fl_t max_dist = 0.0;
	size_t index = 0;
	for (size_t i = 1; i < p.size(); ++i) {
		fl_t dist = distance_to_line(p[i], p.front(), p.back());
		if (dist > max_dist) {
			index = i;
			max_dist = dist;
		}
	}
	if (max_dist > epsilon) {
		ClipperLib::Path set, r;
		set.insert(set.end(), p.begin(), p.begin() + index + 1);
		r = rdp_simplify_path(set, epsilon);
		res.insert(res.end(), r.begin(), r.end() - 1);

		set.clear();
		set.insert(set.end(), p.begin() + index, p.end());
		r = rdp_simplify_path(set, epsilon);
		res.insert(res.end(), r.begin(), r.end());
	}
	else {
		res.push_back(p.front());
		res.push_back(p.back());
	}
	return res;
}

static void simplify_path(ClipperLib::Path &p, const fl_t epsilon)
{
	if (p.size() < 3) return;
	p.push_back(p[0]);      /* the first and last point must be the same for rdp_simplify_path() to work correctly */
	p = rdp_simplify_path(p, epsilon);
	p.erase(p.end() - 1);   /* remove duplicate end point again */
}

static void simplify_paths(ClipperLib::Paths &paths, const fl_t epsilon)
{
	for (ClipperLib::Path &p : paths)
		simplify_path(p, epsilon);
}

static void generate_outlines(struct slice *slice, ssize_t slice_index)
{
	struct segment_list iseg = { NULL, NULL }, oseg = { NULL, NULL };
	const fl_t tolerance_sq = config.tolerance * config.tolerance;
	ClipperLib::Paths outlines;

	for (ssize_t i = 0; i < slice->n_seg; ++i)
		LIST_ADD_TAIL(&iseg, &slice->s[i]);

	while (iseg.head) {
		ssize_t segment_count = 0, flip_count = 0;
		/* ssize_t inexact_count = 0; */
		/* Add first segment to the polygon */
		struct segment *s = iseg.head;
		LIST_REMOVE_HEAD(&iseg);
		LIST_ADD_HEAD(&oseg, s);

		next_segment:
		++segment_count;
		bool flip_points = false;
		fl_t best_dist = FL_T_INF;
		struct segment *best = NULL, *begin = oseg.head, *end = oseg.tail;

		/* Check whether the polygon is closed */
		if (begin != end) {
			if (begin->x[0] == end->x[1] && begin->y[0] == end->y[1])
				goto add_poly;
		}

		/* Link up connected segments */
		LIST_FOREACH(&iseg, s) {
			if (s->x[0] == end->x[1] && s->y[0] == end->y[1]) {
				LIST_REMOVE(&iseg, s);
				LIST_ADD_TAIL(&oseg, s);
				goto next_segment;
			}
			if (s->x[1] == end->x[1] && s->y[1] == end->y[1]) {
				DEBUG("flipped segment %zd at layer %zd\n", segment_count, slice_index + 1);
				++flip_count;
				/* Flip point order */
				fl_t t = s->x[0];
				s->x[0] = s->x[1];
				s->x[1] = t;
				t = s->y[0];
				s->y[0] = s->y[1];
				s->y[1] = t;

				LIST_REMOVE(&iseg, s);
				LIST_ADD_TAIL(&oseg, s);
				goto next_segment;
			}
		}

		/* Find closest segment if none are connected exactly */
		LIST_FOREACH(&iseg, s) {
			fl_t dist0 = (s->x[0] - end->x[1]) * (s->x[0] - end->x[1]) + (s->y[0] - end->y[1]) * (s->y[0] - end->y[1]);
			fl_t dist1 = (s->x[1] - end->x[1]) * (s->x[1] - end->x[1]) + (s->y[1] - end->y[1]) * (s->y[1] - end->y[1]);
			fl_t dist = MINIMUM(dist0, dist1);
			if (dist < best_dist) {
				flip_points = (dist1 < dist0);
				best_dist = dist;
				best = s;
			}
		}

		/* Check whether the polygon is closed (within tolerance_sq and less than best_dist) */
		if (begin != end) {
			fl_t close_dist = (begin->x[0] - end->x[1]) * (begin->x[0] - end->x[1]) + (begin->y[0] - end->y[1]) * (begin->y[0] - end->y[1]);
			if (close_dist <= tolerance_sq && close_dist < best_dist)
				goto add_poly;
		}

		/* Connect nearest segment (within tolerance_sq) */
		if (best && best_dist <= tolerance_sq) {
			s = best;
			if (flip_points) {
				DEBUG("flipped segment %zd at layer %zd\n", segment_count, slice_index + 1);
				++flip_count;
				fl_t t = s->x[0];
				s->x[0] = s->x[1];
				s->x[1] = t;
				t = s->y[0];
				s->y[0] = s->y[1];
				s->y[1] = t;
			}
			LIST_REMOVE(&iseg, s);
			LIST_ADD_TAIL(&oseg, s);
			/* ++inexact_count; */
			goto next_segment;
		}

		/* If there are any segments left and more than one output segment, there is probably a hole in the mesh */
		if (iseg.head && oseg.head != oseg.tail)
			fprintf(stderr, "warning: there is (probably) a hole in the mesh at layer %zd (best_dist = %f)\n", slice_index + 1, sqrt(best_dist));
		goto next_poly;

		add_poly:
		if (oseg.head) {
			ClipperLib::Path poly;
			poly.reserve(segment_count);
			LIST_FOREACH(&oseg, s)
				poly.push_back(FL_T_TO_INTPOINT(s->x[0], s->y[0]));
			if (SIMPLIFY_EPSILON > 0.0)
				simplify_path(poly, SIMPLIFY_EPSILON);
			if (config.combine_all) {
				ClipperLib::Paths polys;
				/* Remove self-intersections */
				ClipperLib::SimplifyPolygon(poly, polys, ClipperLib::pftEvenOdd);
				for (ClipperLib::Path &p : polys) {
					/* Reverse any clockwise outlines */
					if (!ClipperLib::Orientation(p))
						ClipperLib::ReversePath(p);
					outlines.push_back(p);
				}
			}
			else {
				/* If we flipped more than half of the segments, the first segment in the outline was probably oriented incorrectly */
				if (flip_count > segment_count / 2) {
					DEBUG("reversed outline order at layer %zd\n", slice_index + 1);
					ClipperLib::ReversePath(poly);
				}
				outlines.push_back(poly);
			}
			/* DEBUG("inexact_count = %zd / %zd\n", inexact_count, segment_count); */
		}
		next_poly:
		oseg.head = oseg.tail = NULL;
	}
	free(slice->s);
	ClipperLib::SimplifyPolygons(outlines, config.poly_fill_type);
	ClipperLib::PolyTree tree;
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	co.AddPaths(outlines, config.inset_join_type, ClipperLib::etClosedPolygon);
	if (1.0 - config.edge_overlap > 0.0) {
		ClipperLib::Paths tmp;
		fl_t extra = -config.extrusion_width * (1.0 - config.edge_overlap) / 2.0;
		co.Execute(tmp, FL_T_TO_CINT(config.edge_offset + config.extra_offset + extra));
		co.Clear();
		co.AddPaths(tmp, config.outset_join_type, ClipperLib::etClosedPolygon);
		co.Execute(tree, FL_T_TO_CINT(-extra));
	}
	else
		co.Execute(tree, FL_T_TO_CINT(config.edge_offset + config.extra_offset));
	generate_islands(slice, &tree);
	if (config.simplify_insets && SIMPLIFY_EPSILON > 0.0)
		for (struct island &island : slice->islands)
			simplify_paths(island.insets[0], SIMPLIFY_EPSILON);
	for (struct island &island : slice->islands)
		find_bounding_box(&island);
}

static void remove_overlap(ClipperLib::Paths &src, ClipperLib::Paths &dest, fl_t ratio)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	co.AddPaths(src, config.inset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(dest, FL_T_TO_CINT(config.extrusion_width * ratio / -2.0));
	co.Clear();
	co.AddPaths(dest, config.outset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(dest, FL_T_TO_CINT(config.extrusion_width * ratio / 2.0));
}

static void do_offset(ClipperLib::Paths &src, ClipperLib::Paths &dest, fl_t dist, fl_t overlap_removal_ratio)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	co.AddPaths(src, (dist > 0.0) ? config.outset_join_type : config.inset_join_type, ClipperLib::etClosedPolygon);
	if (overlap_removal_ratio > 0.0) {
		fl_t extra = (dist > 0.0) ? (config.extrusion_width * overlap_removal_ratio / 2.0) : (config.extrusion_width * overlap_removal_ratio / -2.0);
		co.Execute(dest, FL_T_TO_CINT(dist + extra));
		co.Clear();
		co.AddPaths(dest, (dist > 0.0) ? config.inset_join_type : config.outset_join_type, ClipperLib::etClosedPolygon);
		co.Execute(dest, FL_T_TO_CINT(-extra));
	}
	else
		co.Execute(dest, FL_T_TO_CINT(dist));
}

static void do_offset_square(ClipperLib::Paths &src, ClipperLib::Paths &dest, fl_t dist, fl_t overlap_removal_ratio)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	co.AddPaths(src, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
	if (overlap_removal_ratio > 0.0) {
		fl_t extra = (dist > 0.0) ? (config.extrusion_width * overlap_removal_ratio / 2.0) : (config.extrusion_width * overlap_removal_ratio / -2.0);
		co.Execute(dest, FL_T_TO_CINT(dist + extra));
		co.Clear();
		co.AddPaths(dest, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
		co.Execute(dest, FL_T_TO_CINT(-extra));
	}
	else
		co.Execute(dest, FL_T_TO_CINT(dist));
}

#define BOUND_OFFSET (config.extrusion_width / 8.0)
#define BOUND_SIMPLIFY_EPSILON (BOUND_OFFSET / 2.0 * config.scale_constant)
static void generate_insets(struct slice *slice)
{
	for (struct island &island : slice->islands) {
		if (config.shells > 0) {
			for (int i = 1; i < config.shells; ++i) {
				do_offset(island.insets[i - 1], island.insets[i], -config.extrusion_width, 1.0);
				if (config.simplify_insets && SIMPLIFY_EPSILON > 0.0)
					simplify_paths(island.insets[i], SIMPLIFY_EPSILON);
				if (island.insets[i].size() == 0)  /* break if nothing is being generated */
					goto done;
			}
			do_offset(island.insets[config.shells - 1], island.infill_insets, (0.5 - config.infill_overlap) * -config.extrusion_width, 0.0);
			if (SIMPLIFY_EPSILON > 0.0)
				simplify_paths(island.infill_insets, SIMPLIFY_EPSILON);
		}
		else {
			/* The offset distance here is not *technically* correct, but I'm not sure one can expect high dimensional accuracy when only printing infill anyway... */
			island.infill_insets = island.insets[0];
		}

		done:
		do_offset(island.insets[0], island.boundaries, BOUND_OFFSET, 0.0);
		simplify_paths(island.boundaries, BOUND_SIMPLIFY_EPSILON);
		if (config.solid_infill_clip_offset > 0.0)
			do_offset(island.infill_insets, island.solid_infill_clip, config.solid_infill_clip_offset, 0.0);
		else
			island.solid_infill_clip = island.infill_insets;
		if (config.comb || config.generate_support) {
			do_offset(island.insets[0], island.outer_boundaries, 0.5 * config.edge_width - config.edge_offset, 0.0);
			simplify_paths(island.outer_boundaries, BOUND_SIMPLIFY_EPSILON);
		}
		if (config.comb) {
			island.comb_paths = island.insets[0];
			do_offset(island.outer_boundaries, island.outer_comb_paths, BOUND_OFFSET, 0.0);
			simplify_paths(island.outer_comb_paths, BOUND_SIMPLIFY_EPSILON);
		}
		if (config.shells > 1 && config.fill_inset_gaps) {
			ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
			ClipperLib::Paths hole;
			island.inset_gaps = new ClipperLib::Paths[config.shells - 1]();
			if (!island.inset_gaps)
				die(e_nomem, 2);
			for (int i = 0; i < config.shells - 1 && island.insets[i].size() > 0; ++i) {
				co.AddPaths(island.insets[i], config.inset_join_type, ClipperLib::etClosedPolygon);
				hole = island.insets[i + 1];
				ClipperLib::ReversePaths(hole);
				co.AddPaths(hole, config.inset_join_type, ClipperLib::etClosedPolygon);
				if (config.fill_threshold > 0.0) {
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT((0.5 + config.fill_threshold / 2.0) * -config.extrusion_width));
					co.Clear();
					co.AddPaths(island.inset_gaps[i], config.outset_join_type, ClipperLib::etClosedPolygon);
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT((config.infill_overlap + config.fill_threshold / 2.0) * config.extrusion_width));
				}
				else {
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT((0.5 - config.infill_overlap) * -config.extrusion_width));
				}
				co.Clear();
			}
		}
		do_offset(island.infill_insets, island.constraining_edge, -BOUND_OFFSET, 0.0);
		if (config.align_seams) {
			for (int i = 0; i < ((config.align_interior_seams) ? config.shells : 1); ++i) {
				for (ClipperLib::Path &p : island.insets[i]) {
					if (p.size() >= 3) {
						fl_t lowest = FL_T_INF;
						auto best = p.begin();
						for (auto it = best; it != p.end(); ++it) {
							fl_t v = (*it).X + (*it).Y;
							if (v < lowest) {
								best = it;
								lowest = v;
							}
						}
						std::rotate(p.begin(), best, p.end());
					}
				}
			}
		}
	}
}

/* FIXME (maybe?): This function uses (0,0) as the origin. Perhaps the object center would be a better choice. */
static void generate_line_fill_at_angle(ClipperLib::Paths &p, fl_t x0, fl_t y0, fl_t x1, fl_t y1, fl_t density, fl_t angle)
{
	const fl_t sin_angle = sin(angle), cos_angle = cos(angle);
	const fl_t sin_neg_angle = sin(-angle), cos_neg_angle = cos(-angle);
	const fl_t move = config.extrusion_width / density;
	/* compute rotated bounding box (all 4 corners) */
	/* upper left */
	const fl_t c0_x = x0 * cos_neg_angle - y0 * sin_neg_angle;
	const fl_t c0_y = x0 * sin_neg_angle + y0 * cos_neg_angle;
	/* lower left */
	const fl_t c1_x = x0 * cos_neg_angle - y1 * sin_neg_angle;
	const fl_t c1_y = x0 * sin_neg_angle + y1 * cos_neg_angle;
	/* lower right */
	const fl_t c2_x = x1 * cos_neg_angle - y1 * sin_neg_angle;
	const fl_t c2_y = x1 * sin_neg_angle + y1 * cos_neg_angle;
	/* upper right */
	const fl_t c3_x = x1 * cos_neg_angle - y0 * sin_neg_angle;
	const fl_t c3_y = x1 * sin_neg_angle + y0 * cos_neg_angle;
	/* find start and end indices */
	const ssize_t start = lround(floor(MINIMUM_4(c0_y, c1_y, c2_y, c3_y) / move));
	const ssize_t end = lround(ceil(MAXIMUM_4(c0_y, c1_y, c2_y, c3_y) / move));
	const fl_t min_x = MINIMUM_4(c0_x, c1_x, c2_x, c3_x);
	const fl_t max_x = MAXIMUM_4(c0_x, c1_x, c2_x, c3_x);
	for (ssize_t i = start; i <= end; ++i) {
		ClipperLib::Path line(2);
		const fl_t y = move * i;
		line[0].X = FL_T_TO_CINT(cos_angle * min_x - sin_angle * y);
		line[0].Y = FL_T_TO_CINT(sin_angle * min_x + cos_angle * y);
		line[1].X = FL_T_TO_CINT(cos_angle * max_x - sin_angle * y);
		line[1].Y = FL_T_TO_CINT(sin_angle * max_x + cos_angle * y);
		p.push_back(line);
	}
}

/* TODO: Generate support patterns per-region instead of in this function. */
static void generate_infill_patterns(struct object *o)
{
	const fl_t x_len_2 = (o->w + config.xy_extra) / 2.0, y_len_2 = (o->d + config.xy_extra) / 2.0;
	const fl_t x0 = o->c.x - x_len_2, y0 = o->c.y - y_len_2, x1 = o->c.x + x_len_2, y1 = o->c.y + y_len_2;
	const fl_t solid_infill_angle_rad = config.solid_infill_angle / 180.0 * M_PI;

	if (config.generate_raft || (config.generate_support && config.solid_support_base)) {
		/* generate_line_fill_at_angle(o->solid_infill_patterns[0], x0, y0, x1, y1, 1.0, solid_infill_angle_rad); */  /* The raft and support code only uses solid_infill_patterns[1] */
		generate_line_fill_at_angle(o->solid_infill_patterns[1], x0, y0, x1, y1, 1.0, solid_infill_angle_rad + M_PI_2);
	}
	if (config.generate_support) {  /* +- 45deg from solid infill angle */
		generate_line_fill_at_angle(o->support_pattern, x0, y0, x1, y1, config.support_density, solid_infill_angle_rad - M_PI_4);
		generate_line_fill_at_angle(o->support_interface_pattern, x0, y0, x1, y1, config.interface_density, solid_infill_angle_rad + M_PI_4);
	}
	if (config.generate_raft)
		generate_line_fill_at_angle(o->raft_base_layer_pattern, x0, y0, x1, y1, (config.extrusion_width / config.raft_base_layer_width) * config.raft_base_layer_density, solid_infill_angle_rad);
}

static void generate_infill_for_box(ClipperLib::Paths &p, const struct cint_rect &box, fl_t density, fl_t angle, fill_pattern pattern, ssize_t slice_index)
{
	if (density > 0.0) {
		const fl_t angle_rad = angle / 180.0 * M_PI;
		const fl_t x0 = CINT_TO_FL_T(box.x0), y0 = CINT_TO_FL_T(box.y0), x1 = CINT_TO_FL_T(box.x1), y1 = CINT_TO_FL_T(box.y1);
		switch (pattern) {
		case FILL_PATTERN_GRID:
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density / 2.0, angle_rad);
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density / 2.0, angle_rad + M_PI_2);
			break;
		case FILL_PATTERN_TRIANGLE:
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density / 3.0, angle_rad);
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density / 3.0, angle_rad + M_PI / 3.0);
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density / 3.0, angle_rad + 2.0 * M_PI / 3.0);
			break;
		case FILL_PATTERN_TRIANGLE2:
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density, angle_rad + (fl_t) slice_index * M_PI / 3.0);
			break;
		case FILL_PATTERN_RECTILINEAR:
		default:
			generate_line_fill_at_angle(p, x0, y0, x1, y1, density, angle_rad + (fl_t) slice_index * M_PI / 2.0);
		}
	}
}

#define BOUNDING_BOX_INTERSECTS(a, b) (!((b).x0 > (a).x1 || (b).x1 < (a).x0 || (b).y0 < (a).y1 || (b).y1 > (a).y0))

static void generate_infill(struct object *o, ssize_t slice_index)
{
	for (struct island &island : o->slices[slice_index].islands) {
		ClipperLib::Clipper c;
		ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);  /* for building island.solid_infill_boundaries */
		ClipperLib::PolyTree s;
		ClipperLib::Paths s_tmp;
		ClipperLib::Paths solid_infill_pattern, sparse_infill_pattern;
		if (config.roof_layers > 0) {
			if (slice_index + 1 == o->n_slices)
				island.exposed_surface = island.infill_insets;
			else {
				c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
				for (const struct island &clip_island : o->slices[slice_index + 1].islands)
					if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
						c.AddPaths(clip_island.insets[0], ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctDifference, island.exposed_surface, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
			}
			if (island.exposed_surface.size() > 0)
				do_offset(island.exposed_surface, island.exposed_surface, -config.extrusion_width, 0.0);
			if (config.iron_top_surface) {
				ClipperLib::Paths iron_areas;
				if (slice_index + 1 == o->n_slices)
					do_offset(island.insets[0], iron_areas, -config.extrusion_width / 2.0, 0.0);
				else {
					do_offset(island.insets[0], iron_areas, -config.extrusion_width / 2.0, 0.0);
					c.AddPaths(iron_areas, ClipperLib::ptSubject, true);
					for (const struct island &clip_island : o->slices[slice_index + 1].islands)
						if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
							c.AddPaths(clip_island.insets[0], ClipperLib::ptClip, true);
					c.Execute(ClipperLib::ctDifference, iron_areas, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
					c.Clear();
				}
				if (iron_areas.size() > 0) {
					remove_overlap(iron_areas, iron_areas, 1.0);
					ClipperLib::Paths iron_pattern;
					generate_infill_for_box(iron_pattern, island.box, config.iron_density, config.solid_infill_angle, FILL_PATTERN_RECTILINEAR, slice_index + 1);
					c.AddPaths(iron_pattern, ClipperLib::ptSubject, false);
					c.AddPaths(iron_areas, ClipperLib::ptClip, true);
					c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
					c.Clear();
					ClipperLib::OpenPathsFromPolyTree(s, island.iron_paths);
					s.Clear();
				}
			}
		}
		if (config.infill_density == 1.0 || slice_index < config.floor_layers || slice_index + config.roof_layers >= o->n_slices) {
			if (config.fill_threshold > 0.0) {
				remove_overlap(island.infill_insets, s_tmp, config.fill_threshold);
				c.AddPaths(s_tmp, ClipperLib::ptClip, true);
				co.AddPaths(s_tmp, config.outset_join_type, ClipperLib::etClosedPolygon);
			}
			else {
				c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
				co.AddPaths(island.infill_insets, config.outset_join_type, ClipperLib::etClosedPolygon);
			}
			generate_infill_for_box(solid_infill_pattern, island.box, 1.0, config.solid_infill_angle, FILL_PATTERN_RECTILINEAR, slice_index);
			c.AddPaths(solid_infill_pattern, ClipperLib::ptSubject, false);
			if (config.fill_inset_gaps) {
				for (int i = 0; i < config.shells - 1; ++i) {
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
					co.AddPaths(island.inset_gaps[i], config.outset_join_type, ClipperLib::etClosedPolygon);
				}
			}
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
			co.Execute(island.solid_infill_boundaries, FL_T_TO_CINT(BOUND_OFFSET));
			simplify_paths(island.solid_infill_boundaries, BOUND_SIMPLIFY_EPSILON);
		}
		else if (!config.no_solid && (config.floor_layers > 0 || config.roof_layers > 0)) {
			c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			for (int i = -config.floor_layers; i <= config.roof_layers; ++i) {
				if (i != 0) {
					for (const struct island &clip_island : o->slices[slice_index + i].islands)
						if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
							c.AddPaths(clip_island.solid_infill_clip, ClipperLib::ptClip, true);
					c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
					c.Clear();
					if (i != config.roof_layers)
						c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
				}
			}
			c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			c.Clear();
			if (config.fill_threshold > 0.0)
				remove_overlap(s_tmp, s_tmp, config.fill_threshold);
			if (config.solid_fill_expansion > 0.0 || config.solid_infill_clip_offset > 0.0) {
				do_offset_square(s_tmp, s_tmp, config.solid_infill_clip_offset + config.solid_fill_expansion * config.extrusion_width, 0.0);
				c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
				c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
			}
			generate_infill_for_box(solid_infill_pattern, island.box, 1.0, config.solid_infill_angle, FILL_PATTERN_RECTILINEAR, slice_index);
			c.AddPaths(solid_infill_pattern, ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			co.AddPaths(s_tmp, config.outset_join_type, ClipperLib::etClosedPolygon);
			if (config.fill_inset_gaps) {
				for (int i = 0; i < config.shells - 1; ++i) {
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
					co.AddPaths(island.inset_gaps[i], config.outset_join_type, ClipperLib::etClosedPolygon);
				}
			}
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			c.Clear();
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
			co.Execute(island.solid_infill_boundaries, FL_T_TO_CINT(BOUND_OFFSET));
			simplify_paths(island.solid_infill_boundaries, BOUND_SIMPLIFY_EPSILON);

			if (config.infill_density > 0.0) {
				c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
				c.AddPaths(s_tmp, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
				if (config.fill_threshold > 0.0)
					remove_overlap(s_tmp, s_tmp, config.fill_threshold);
				generate_infill_for_box(sparse_infill_pattern, island.box, config.infill_density, config.sparse_infill_angle, config.infill_pattern, slice_index);
				c.AddPaths(sparse_infill_pattern, ClipperLib::ptSubject, false);
				c.AddPaths(s_tmp, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				ClipperLib::OpenPathsFromPolyTree(s, island.sparse_infill);
			}
		}
		else {
			if (config.infill_density > 0.0) {
				if (config.fill_threshold > 0.0) {
					remove_overlap(island.infill_insets, s_tmp, config.fill_threshold);
					c.AddPaths(s_tmp, ClipperLib::ptClip, true);
				}
				else
					c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
				generate_infill_for_box(sparse_infill_pattern, island.box, config.infill_density, config.sparse_infill_angle, config.infill_pattern, slice_index);
				c.AddPaths(sparse_infill_pattern, ClipperLib::ptSubject, false);
				c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				ClipperLib::OpenPathsFromPolyTree(s, island.sparse_infill);
			}
			if (config.fill_inset_gaps) {
				c.Clear();
				generate_infill_for_box(solid_infill_pattern, island.box, 1.0, config.solid_infill_angle, FILL_PATTERN_RECTILINEAR, slice_index);
				c.AddPaths(solid_infill_pattern, ClipperLib::ptSubject, false);
				for (int i = 0; i < config.shells - 1; ++i) {
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
					co.AddPaths(island.inset_gaps[i], config.outset_join_type, ClipperLib::etClosedPolygon);
				}
				c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
				co.Execute(island.solid_infill_boundaries, FL_T_TO_CINT(BOUND_OFFSET));
				simplify_paths(island.solid_infill_boundaries, BOUND_SIMPLIFY_EPSILON);
			}
		}
		if (config.min_sparse_infill_len > 0.0) {
			const fl_t min_len = config.min_sparse_infill_len * config.scale_constant;
			for (size_t i = 0; i < island.sparse_infill.size();) {
				if (distance_to_point(island.sparse_infill[i][0], island.sparse_infill[i][1]) < min_len)
					island.sparse_infill.erase(island.sparse_infill.begin() + i);
				else
					++i;
			}
		}
	}
}

static void generate_layer_support_map(struct object *o, ssize_t slice_index)
{
	if (slice_index < config.support_vert_margin + 1)
		return;
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	ClipperLib::Clipper c;
	ClipperLib::Paths clip_paths;
	for (const struct island &island : o->slices[slice_index - 1].islands)
		co.AddPaths(island.insets[0], config.outset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(clip_paths, FL_T_TO_CINT(tan(config.support_angle / 180.0 * M_PI) * config.layer_height));
	co.Clear();
	for (const struct island &island : o->slices[slice_index].islands)
		c.AddPaths(island.insets[0], ClipperLib::ptSubject, true);
	c.AddPaths(clip_paths, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctDifference, clip_paths, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	c.Clear();
	co.AddPaths(clip_paths, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
	co.Execute(o->slices[slice_index].layer_support_map, FL_T_TO_CINT(config.support_xy_expansion + (0.5 + config.support_margin) * config.edge_width - config.edge_offset));
}

static void generate_support_boundaries(struct slice *slice)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	for (const struct island &island : slice->islands)
		co.AddPaths(island.insets[0], config.outset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(slice->support_boundaries, FL_T_TO_CINT((0.5 + config.support_margin) * config.edge_width - config.edge_offset));
	simplify_paths(slice->support_boundaries, BOUND_SIMPLIFY_EPSILON);
}

static void extend_support_downward(struct object *o, const ClipperLib::PolyNode *n, ssize_t slice_index)
{
	ssize_t k;
	ClipperLib::Paths p;
	p.push_back(n->Contour);
	for (const ClipperLib::PolyNode *c : n->Childs)
		p.push_back(c->Contour);
	ClipperLib::Paths *clipped_paths = new ClipperLib::Paths[slice_index + 1];
	for (k = slice_index; k >= 0; --k) {
		ClipperLib::Clipper c;
		c.AddPaths(p, ClipperLib::ptSubject, true);
		for (ssize_t i = (k >= config.support_vert_margin) ? -config.support_vert_margin : -k; k + i <= o->n_slices && i <= config.support_vert_margin; ++i)
			c.AddPaths(o->slices[k + i].support_boundaries, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, clipped_paths[k], ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		if (clipped_paths[k].size() == 0)
			break;
	}
	if (config.support_everywhere || k == -1) {
		for (++k; k <= slice_index; ++k) {
		#ifdef _OPENMP
			omp_set_lock(&o->slices[k].lock);
		#endif
			o->slices[k].support_map.insert(o->slices[k].support_map.end(), clipped_paths[k].begin(), clipped_paths[k].end());
		#ifdef _OPENMP
			omp_unset_lock(&o->slices[k].lock);
		#endif
		}
	}
	delete[] clipped_paths;
}

static void generate_support_maps(struct object *o, const ClipperLib::PolyNode *n, ssize_t slice_index)
{
	for (const ClipperLib::PolyNode *c : n->Childs) {
		extend_support_downward(o, c, slice_index);
		for (const ClipperLib::PolyNode *cc : c->Childs)
			generate_support_maps(o, cc, slice_index);
	}
}

static void union_support_maps(struct slice *slice)
{
	ClipperLib::SimplifyPolygons(slice->support_map, ClipperLib::pftNonZero);
}

static void remove_supports_not_touching_build_plate(struct object *o)
{
	if (o->n_slices < 1)
		return;
	ClipperLib::Clipper c;
	ClipperLib::Paths clip_paths;
	clip_paths.insert(clip_paths.end(), o->slices[0].support_boundaries.begin(), o->slices[0].support_boundaries.end());
	for (ssize_t i = 1; i < o->n_slices; ++i) {
		clip_paths.insert(clip_paths.end(), o->slices[i].support_boundaries.begin(), o->slices[i].support_boundaries.end());
		ClipperLib::SimplifyPolygons(clip_paths, ClipperLib::pftNonZero);
		c.AddPaths(o->slices[i].support_map, ClipperLib::ptSubject, true);
		c.AddPaths(clip_paths, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, o->slices[i].support_map, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
	}
}

static void generate_support_interface_clip_regions(struct slice *slice)
{
	do_offset_square(slice->support_map, slice->support_interface_clip, config.interface_clip_offset, 0.0);
}

static void generate_support_lines(struct object *o, struct slice *slice, ssize_t slice_index)
{
	ClipperLib::Clipper c;
	ClipperLib::PolyTree s;
	if (config.solid_support_base && slice_index == 0) {
		c.AddPaths(o->solid_infill_patterns[1], ClipperLib::ptSubject, false);
		c.AddPaths(slice->support_map, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		ClipperLib::OpenPathsFromPolyTree(s, slice->support_interface_lines);
	}
	else if (config.interface_roof_layers > 0 || config.interface_floor_layers > 0) {
		ClipperLib::Paths s_tmp;
		c.AddPaths(slice->support_map, ClipperLib::ptSubject, true);
		for (int i = (slice_index > config.interface_floor_layers) ? -config.interface_floor_layers : -slice_index; slice_index + i < o->n_slices && i <= config.interface_roof_layers; ++i) {
			if (i != 0) {
				c.AddPaths(o->slices[slice_index + i].support_interface_clip, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
				if (i < config.interface_roof_layers)
					c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
			}
		}
		c.Clear();
		c.AddPaths(slice->support_map, ClipperLib::ptSubject, true);
		c.AddPaths(s_tmp, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
		if (config.expand_support_interface) {
			do_offset_square(s_tmp, s_tmp, config.extrusion_width / config.support_density, 0.0);
			c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
			c.AddPaths(slice->support_map, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			c.Clear();
		}
		c.AddPaths(o->support_interface_pattern, ClipperLib::ptSubject, false);
		c.AddPaths(s_tmp, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
		ClipperLib::OpenPathsFromPolyTree(s, slice->support_interface_lines);

		c.AddPaths(slice->support_map, ClipperLib::ptSubject, true);
		c.AddPaths(s_tmp, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
		c.AddPaths(o->support_pattern, ClipperLib::ptSubject, false);
		c.AddPaths(s_tmp, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		ClipperLib::OpenPathsFromPolyTree(s, slice->support_lines);
	}
	else {
		c.AddPaths(o->support_pattern, ClipperLib::ptSubject, false);
		c.AddPaths(slice->support_map, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		ClipperLib::OpenPathsFromPolyTree(s, slice->support_lines);
	}
}

static void generate_brim(struct object *o)
{
	if (o->n_slices < 1)
		return;
	o->brim.reserve(config.brim_lines);
	for (int i = 1; i <= config.brim_lines; ++i) {
		ClipperLib::Paths tmp;
		for (const struct island &island : o->slices[0].islands)
			tmp.insert(tmp.end(), island.insets[0].begin(), island.insets[0].end());
		if (config.generate_support) {
			tmp.insert(tmp.end(), o->slices[0].support_map.begin(), o->slices[0].support_map.end());
			ClipperLib::SimplifyPolygons(tmp, ClipperLib::pftNonZero);
		}
		do_offset_square(tmp, tmp, config.extrusion_width * i + (config.edge_offset * -2.0 - config.extrusion_width) * (1.0 - config.brim_adhesion_factor) * 2.0, 1.0);
		if (SIMPLIFY_EPSILON > 0.0)
			simplify_paths(tmp, SIMPLIFY_EPSILON);
		o->brim.push_back(tmp);
	}
}

static void generate_raft(struct object *o)
{
	if (o->n_slices < 1)
		return;
	ClipperLib::Clipper c;
	ClipperLib::PolyTree s;
	ClipperLib::Paths tmp;
	if (config.brim_lines > 0) {
		for (const ClipperLib::Paths &p : o->brim)
			tmp.insert(tmp.end(), p.begin(), p.end());
		ClipperLib::SimplifyPolygons(tmp, ClipperLib::pftNonZero);
	}
	else {
		for (const struct island &island : o->slices[0].islands)
			tmp.insert(tmp.end(), island.insets[0].begin(), island.insets[0].end());
		if (config.generate_support) {
			tmp.insert(tmp.end(), o->slices[0].support_map.begin(), o->slices[0].support_map.end());
			ClipperLib::SimplifyPolygons(tmp, ClipperLib::pftNonZero);
		}
	}
	do_offset_square(tmp, tmp, config.raft_xy_expansion, 0.0);
	c.AddPaths(o->raft_base_layer_pattern, ClipperLib::ptSubject, false);
	c.AddPaths(tmp, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	ClipperLib::OpenPathsFromPolyTree(s, o->raft[0]);
	c.Clear();
	c.AddPaths(o->solid_infill_patterns[1], ClipperLib::ptSubject, false);
	c.AddPaths(tmp, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	ClipperLib::OpenPathsFromPolyTree(s, o->raft[1]);
}

static void slice_object(struct object *o)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	ssize_t i;
	o->n_slices = (ssize_t) ceil((o->c.z + o->h / 2.0) / config.layer_height);
	o->slices = new struct slice[o->n_slices]();
	if (!o->slices)
		die(e_nomem, 2);

	start = std::chrono::high_resolution_clock::now();
	fputs("  find segments...", stderr);
#ifdef _OPENMP
	for (i = 0; i < o->n_slices; ++i)
		omp_init_lock(&o->slices[i].lock);
	#pragma omp parallel for
#endif
	for (i = 0; i < o->n; ++i)
		find_segments(o->slices, &o->t[i]);
	fputs(" done\n", stderr);
	free(o->t);
	fputs("  generate outlines...", stderr);
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < o->n_slices; ++i)
		generate_outlines(&o->slices[i], i);
	fputs(" done\n", stderr);
	fputs("  generate insets...", stderr);
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < o->n_slices; ++i)
		generate_insets(&o->slices[i]);
	fputs(" done\n", stderr);
	fputs("  generate infill...", stderr);
	generate_infill_patterns(o);
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (i = 0; i < o->n_slices; ++i)
		generate_infill(o, i);
	fputs(" done\n", stderr);

	if (config.generate_support) {
		fputs("  generate support...", stderr);
	#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
	#endif
		for (i = 0; i < o->n_slices; ++i) {
			generate_layer_support_map(o, i);
			generate_support_boundaries(&o->slices[i]);
		}
	#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
	#endif
		for (i = 0; i < o->n_slices; ++i)
			generate_support_maps(o, &o->slices[i].layer_support_map, i);
	#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
	#endif
		for (i = 0; i < o->n_slices; ++i)
			union_support_maps(&o->slices[i]);
		if (!config.support_everywhere)
			remove_supports_not_touching_build_plate(o);
		if (config.interface_roof_layers > 0 || config.interface_floor_layers > 0) {
		#ifdef _OPENMP
			#pragma omp parallel for schedule(dynamic)
		#endif
			for (i = 0; i < o->n_slices; ++i)
				generate_support_interface_clip_regions(&o->slices[i]);
		}
	#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
	#endif
		for (i = 0; i < o->n_slices; ++i)
			generate_support_lines(o, &o->slices[i], i);
		fputs(" done\n", stderr);
	}
	if (config.brim_lines > 0) {
		fputs("  generate brim...", stderr);
		generate_brim(o);
		fputs(" done\n", stderr);
	}
	if (config.generate_raft) {
		fputs("  generate raft...", stderr);
		generate_raft(o);
		fputs(" done\n", stderr);
	}
#ifdef _OPENMP
	for (i = 0; i < o->n_slices; ++i)
		omp_destroy_lock(&o->slices[i].lock);
#endif

	fprintf(stderr, "sliced in %fs\n",
		(double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0);
}

static void preview_slices(const struct object *o)
{
	ssize_t i;
	fprintf(stdout, "set size ratio -1\n");
	fprintf(stdout, "set xrange [%e:%e]\n", o->c.x - (o->w + config.xy_extra) / 2.0, o->c.x + (o->w + config.xy_extra) / 2.0);
	fprintf(stdout, "set yrange [%e:%e]\n", o->c.y - (o->d + config.xy_extra) / 2.0, o->c.y + (o->d + config.xy_extra) / 2.0);
	for (i = 0; i < o->n_slices; ++i) {
		fprintf(stderr, "layer %zd/%zd: intersections = %zd; islands = %zd\n", i + 1, o->n_slices, o->slices[i].n_seg, o->slices[i].islands.size());
		fprintf(stdout, "set title 'Layer %zd/%zd'\n", i + 1, o->n_slices);
		/* Draw outlines */
		fprintf(stdout, "plot \"-\" lc \"red\" dt 3 with lines title \"boundaries\", \"-\" lc \"blue\" with lines title \"insets\", \"-\" lc \"green\" with lines title \"infill\"\n");
		for (const struct island &island : o->slices[i].islands) {
			for (const ClipperLib::Path &path : island.boundaries) {
				if (path.size() >= 3) {
					for (const ClipperLib::IntPoint &p : path)
						fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
					fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / config.scale_constant, ((double) path[0].Y) / config.scale_constant);
					putc('\n', stdout);
				}
			}
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x0) / config.scale_constant, ((double) island.box.y0) / config.scale_constant);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x1) / config.scale_constant, ((double) island.box.y0) / config.scale_constant);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x1) / config.scale_constant, ((double) island.box.y1) / config.scale_constant);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x0) / config.scale_constant, ((double) island.box.y1) / config.scale_constant);
			fprintf(stdout, "%.4e %.4e\n\n", ((double) island.box.x0) / config.scale_constant, ((double) island.box.y0) / config.scale_constant);
		}
		/* Draw support map */
		for (const ClipperLib::Path &path : o->slices[i].support_map) {
			if (path.size() >= 3) {
				for (const ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
				fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / config.scale_constant, ((double) path[0].Y) / config.scale_constant);
				putc('\n', stdout);
			}
		}
		fprintf(stdout, "e\n");
		/* Draw insets */
		for (const struct island &island : o->slices[i].islands) {
			for (int k = 0; k < config.shells; ++k) {
				for (const ClipperLib::Path &path : island.insets[k]) {
					if (path.size() >= 3) {
						for (const ClipperLib::IntPoint &p : path)
							fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
						fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / config.scale_constant, ((double) path[0].Y) / config.scale_constant);
						putc('\n', stdout);
					}
				}
			}
		}
		/* Draw brim */
		if (i == 0) {
			for (const ClipperLib::Paths &paths : o->brim) {
				for (const ClipperLib::Path &path : paths) {
					if (path.size() >= 3) {
						for (const ClipperLib::IntPoint &p : path)
							fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
						fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / config.scale_constant, ((double) path[0].Y) / config.scale_constant);
						putc('\n', stdout);
					}
				}
			}
		}
		fprintf(stdout, "e\n");
		/* Draw infill */
		for (const struct island &island : o->slices[i].islands) {
			for (const ClipperLib::Path &path : island.solid_infill) {
				for (const ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
				putc('\n', stdout);
			}
			for (const ClipperLib::Path &path : island.sparse_infill) {
				for (const ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
				putc('\n', stdout);
			}
		}
		/* Draw support lines */
		for (const ClipperLib::Path &path : o->slices[i].support_lines) {
			for (const ClipperLib::IntPoint &p : path)
				fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
			putc('\n', stdout);
		}
		for (const ClipperLib::Path &path : o->slices[i].support_interface_lines) {
			for (const ClipperLib::IntPoint &p : path)
				fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / config.scale_constant, ((double) p.Y) / config.scale_constant);
			putc('\n', stdout);
		}
		fprintf(stdout, "e\n");
		fflush(stdout);
		fprintf(stderr, "press enter for next layer...");
		getchar();
	}
}

/* 0 is colinear, 1 is counter-clockwise and -1 is clockwise */
static int triplet_orientation(const ClipperLib::IntPoint &a, const ClipperLib::IntPoint &b, const ClipperLib::IntPoint &c)
{
	const ClipperLib::cInt v = (b.X - a.X) * (c.Y - a.Y) - (c.X - a.X) * (b.Y - a.Y);  /* Calculate signed area * 2 */
	return (v == 0) ? 0 : (v > 0) ? 1 : -1;
}

static int is_on_segment(const ClipperLib::IntPoint &a, const ClipperLib::IntPoint &b, const ClipperLib::IntPoint &c)
{
	if (b.X <= MAXIMUM(a.X, c.X) && b.X >= MINIMUM(a.X, c.X) && b.Y <= MAXIMUM(a.Y, c.Y) && b.Y >= MINIMUM(a.Y, c.Y))
		return true;
	return false;
}

static bool intersects(const ClipperLib::IntPoint &a, const ClipperLib::IntPoint &b, const ClipperLib::IntPoint &c, const ClipperLib::IntPoint &d)
{
	const int o1 = triplet_orientation(a, b, c);
	const int o2 = triplet_orientation(a, b, d);
	const int o3 = triplet_orientation(c, d, a);
	const int o4 = triplet_orientation(c, d, b);
	if (o1 != o2 && o3 != o4)
		return true;
	if (o1 == 0 && is_on_segment(a, c, b))
		return true;
	if (o2 == 0 && is_on_segment(a, d, b))
		return true;
	if (o3 == 0 && is_on_segment(c, a, d))
		return true;
	if (o4 == 0 && is_on_segment(c, b, d))
		return true;
	return false;
}

static ssize_t get_boundary_crossing(const ClipperLib::Path &p, const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1)
{
	for (size_t i = 1; i < p.size(); ++i) {
		if (intersects(p[i - 1], p[i], p0, p1))
			return (ssize_t) i - 1;
	}
	if (intersects(p[p.size() - 1], p[0], p0, p1))
		return p.size() - 1;
	return -1;
}

static ssize_t crosses_boundary(const struct machine *m, const ClipperLib::Paths &bounds, ClipperLib::cInt x, ClipperLib::cInt y)
{
	const ClipperLib::IntPoint p0(m->x, m->y);
	const ClipperLib::IntPoint p1(x, y);
	for (size_t i = 0; i < bounds.size(); ++i) {
		if (get_boundary_crossing(bounds[i], p0, p1) >= 0)
			return (ssize_t) i;
	}
	return -1;
}

static bool crosses_exposed_surface(const struct island *island, ClipperLib::cInt x0, ClipperLib::cInt y0, ClipperLib::cInt x1, ClipperLib::cInt y1)
{
	const ClipperLib::IntPoint p0(x0, y0);
	const ClipperLib::IntPoint p1(x1, y1);
	bool in_outer = false;
	for (const ClipperLib::Path &p : island->exposed_surface) {
		if (get_boundary_crossing(p, p0, p1) >= 0)
			return true;
		else if (ClipperLib::PointInPolygon(p0, p) || ClipperLib::PointInPolygon(p1, p)) {
			if (!ClipperLib::Orientation(p))
				return false;  /* p0 and p1 are inside a hole, so no exposed surface was crossed */
			else
				in_outer = true;
		}
	}
	return in_outer;
}

/* Looks for the nearest line segment and then returns the nearest endpoint on that segment. */
static size_t find_nearest_segment_endpoint_on_closed_path(const ClipperLib::Path &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist)
{
	size_t best = 0;
	fl_t best_dist = FL_T_INF;
	const ClipperLib::IntPoint p0(x, y);
	for (size_t i = 0, i2 = 1; i < p.size(); ++i) {
		i2 = (i2 == p.size()) ? 0 : i2;
		const fl_t dist = distance_to_line(p0, p[i], p[i2]);
		if (dist < best_dist) {
			best_dist = dist;
			best = (distance_to_point(p0, p[i]) < distance_to_point(p0, p[i2])) ? i : i2;
		}
		++i2;
	}
	if (r_dist)
		*r_dist = distance_to_point(p0, p[best]) / config.scale_constant;
	return best;
}

static size_t find_nearest_path(const ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist, size_t *r_start)
{
	size_t best = 0, start = 0;
	fl_t best_dist = FL_T_INF;
	const fl_t x0 = CINT_TO_FL_T(x), y0 = CINT_TO_FL_T(y);
	for (size_t i = 0; i < p.size(); ++i) {
		for (size_t k = 0; k < p[i].size(); ++k) {
			fl_t x1 = CINT_TO_FL_T(p[i][k].X), y1 = CINT_TO_FL_T(p[i][k].Y);
			fl_t dist = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
			if (dist < best_dist) {
				best_dist = dist;
				best = i;
				start = k;
			}
		}
	}
	if (r_dist)
		*r_dist = sqrt(best_dist);
	if (r_start)
		*r_start = start;
	return best;
}

static size_t find_nearest_aligned_path(const ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist)
{
	size_t best = 0;
	fl_t best_dist = FL_T_INF;
	const fl_t x0 = CINT_TO_FL_T(x), y0 = CINT_TO_FL_T(y);
	for (size_t i = 0; i < p.size(); ++i) {
		fl_t x1 = CINT_TO_FL_T(p[i][0].X), y1 = CINT_TO_FL_T(p[i][0].Y);
		fl_t dist = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
		if (dist < best_dist) {
			best_dist = dist;
			best = i;
		}
	}
	if (r_dist)
		*r_dist = sqrt(best_dist);
	return best;
}

static size_t find_nearest_segment(const ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist, bool *r_flip)
{
	size_t best = 0;
	fl_t best_dist = FL_T_INF;
	const ClipperLib::IntPoint p0(x, y);
	for (size_t i = 0; i < p.size(); ++i) {
		if (p[i].size() > 2)
			fprintf(stderr, "error: bug in clipper: line segment has more than two points!\n");
		const fl_t dist = distance_to_line(p0, p[i][0], p[i][1]);
		if (dist < best_dist) {
			best_dist = dist;
			best = i;
		}
	}
	const fl_t dist0 = distance_to_point(p0, p[best][0]);
	const fl_t dist1 = distance_to_point(p0, p[best][1]);
	const bool flip = dist0 > dist1;
	if (r_dist)
		*r_dist = ((flip) ? dist1 : dist0) / config.scale_constant;
	if (r_flip)
		*r_flip = flip;
	return best;
}

static void append_g_move(struct slice *slice, const struct g_move &move, fl_t len)
{
	/* FIXME: should probably take acceleration into account... */
	slice->layer_time += len / move.feed_rate;
	slice->moves.push_back(move);
}

static fl_t get_partial_path_len(const ClipperLib::Path &p, size_t start, size_t end, bool reverse)
{
	fl_t l = 0.0;
	fl_t x0 = CINT_TO_FL_T(p[start].X), y0 = CINT_TO_FL_T(p[start].Y);
	size_t i = start;
	do {
		if (reverse) i = (i > 0) ? i - 1 : p.size() - 1;
		else i = (i < p.size() - 1) ? i + 1 : 0;
		fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		x0 = x1;
		y0 = y1;
	} while (i != end);
	return l;
}

static bool crosses_boundary_2pt(const ClipperLib::Path &p, const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1, fl_t *r_dist)
{
	fl_t best_dist = FL_T_INF;
	size_t intersections = 0;
	for (size_t k = 0; k < p.size(); ++k) {
		const size_t k2 = (k == 0) ? p.size() - 1 : k - 1;
		if (intersects(p[k2], p[k], p0, p1)) {
			const fl_t dist = distance_to_line(p0, p[k2], p[k]);
			if (dist < best_dist)
				best_dist = dist;
			++intersections;
			if (p[k] == p0 || p[k] == p1)
				++k;  /* So the same point isn't registered as two intersections */
		}
	}
	if (r_dist)
		*r_dist = best_dist / config.scale_constant;
	return (intersections > 1);
}

static ssize_t nearest_boundary_crossing_2pt(const ClipperLib::Paths &b, const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1)
{
	fl_t best_dist = FL_T_INF;
	ssize_t b_idx = -1;
	for (size_t i = 0; i < b.size(); ++i) {
		/* Find all intersections, selecting the nearest one. Ignore if there are fewer than two intersections. */
		fl_t tmp_dist = FL_T_INF;
		if (crosses_boundary_2pt(b[i], p0, p1, &tmp_dist) && tmp_dist < best_dist) {
			b_idx = (ssize_t) i;
			best_dist = tmp_dist;
		}
	}
	return b_idx;
}

static size_t find_best_travel_point(const ClipperLib::Paths &b, size_t b_idx, const ClipperLib::IntPoint &p0, size_t start_idx, size_t end_idx, bool reverse)
{
	const ClipperLib::Path &p = b[b_idx];
	size_t i = end_idx, r;
	do {
		r = i;
		if (nearest_boundary_crossing_2pt(b, p0, p[i]) < 0) return i;
		if (reverse) i = (i < p.size() - 1) ? i + 1 : 0;
		else i = (i > 0) ? i - 1 : p.size() - 1;
	} while (i != start_idx);
	return r;
}

static void do_retract(struct slice *slice, struct machine *m)
{
	if (!m->is_retracted && config.retract_len > 0.0) {
		struct g_move retract_move = { m->x, m->y, m->z, -config.retract_len, config.retract_speed, false, false, false };
		append_g_move(slice, retract_move, config.retract_len);
		m->is_retracted = true;
	}
}

static void append_linear_travel(struct slice *slice, struct machine *m, ClipperLib::cInt x, ClipperLib::cInt y, ClipperLib::cInt z, fl_t feed_rate)
{
	if (x != m->x || y != m->y || z != m->z) {
		const fl_t f_x = CINT_TO_FL_T(x), f_y = CINT_TO_FL_T(y), f_z = CINT_TO_FL_T(z);
		const fl_t f_mx = CINT_TO_FL_T(m->x), f_my = CINT_TO_FL_T(m->y), f_mz = CINT_TO_FL_T(m->z);
		const fl_t len = sqrt((f_mx - f_x) * (f_mx - f_x) + (f_my - f_y) * (f_my - f_y) + (f_mz - f_z) * (f_mz - f_z));
		const struct g_move move = { x, y, z, 0.0, feed_rate, false, true, false };
		append_g_move(slice, move, len);
		m->x = x;
		m->y = y;
		m->z = z;
	}
}

static fl_t append_comb_move(const struct machine *m, const struct island *island, ClipperLib::Path &comb_moves, const ClipperLib::IntPoint &p0, const ClipperLib::IntPoint &p1, bool *force_retract)
{
	if (!(*force_retract) && !m->is_retracted && island && crosses_exposed_surface(island, p0.X, p0.Y, p1.X, p1.Y))
		*force_retract = true;
	comb_moves.push_back(p1);
	return distance_to_point(p0, p1) / config.scale_constant;
}

static void combed_travel(struct slice *slice, const struct island *island, struct machine *m, const ClipperLib::Paths &bounds, const ClipperLib::Paths &paths, ClipperLib::cInt x, ClipperLib::cInt y, fl_t feed_rate, fl_t retract_threshold)
{
	if (x == m->x || y == m->y || paths.size() == 0)
		return;
	ClipperLib::Paths b = bounds;
	ssize_t last_bound_idx = -1;
	fl_t closest_dist = FL_T_INF, comb_dist = 0.0;
	ClipperLib::IntPoint p0(m->x, m->y), p1(x, y);
	ClipperLib::Path comb_moves;
	bool force_retract = false;

	while (b.size() > 0) {
		ssize_t bound_idx = nearest_boundary_crossing_2pt(b, p0, p1);
		if (bound_idx < 0)  /* No boundary crossings, so we can move directly there */
			break;
		if (bound_idx == last_bound_idx) {
			b.erase(b.begin() + bound_idx);
			last_bound_idx = -1;
			force_retract = true;
			DEBUG("combed_travel(): warning: removed a boundary at z = %f\n", CINT_TO_FL_T(m->z));
			continue;
		}

		ClipperLib::Path &p = b[bound_idx];
		/* Find the boundary point closest to the start point */
		size_t start_idx = find_nearest_segment_endpoint_on_closed_path(p, p0.X, p0.Y, NULL);
		/* Find the boundary point closest to the end point */
		size_t end_idx = find_nearest_segment_endpoint_on_closed_path(p, x, y, NULL);
		/* fprintf(stderr, "bound_idx=%zd last_bound_idx=%zd start_idx=%zd end_idx=%zd\n", bound_idx, last_bound_idx, start_idx, end_idx); */
		if (distance_to_point(p[end_idx], p1) >= closest_dist) {
			b.erase(b.begin() + bound_idx);
			last_bound_idx = -1;
			force_retract = true;
			DEBUG("combed_travel(): warning: useless indirection at z = %f\n", CINT_TO_FL_T(m->z));
			continue;
		}
		if (start_idx == end_idx) {
			size_t path_pt_idx, path_idx = find_nearest_path(paths, p[end_idx].X, p[end_idx].Y, NULL, &path_pt_idx);
			comb_dist += append_comb_move(m, island, comb_moves, p0, paths[path_idx][path_pt_idx], &force_retract);
			p0 = paths[path_idx][path_pt_idx];
		}
		else {
			/* Find shortest direction */
			const fl_t f_len = get_partial_path_len(p, start_idx, end_idx, false);
			const fl_t r_len = get_partial_path_len(p, start_idx, end_idx, true);
			bool reverse = (r_len < f_len) ? true : false;

			size_t i;
			/* Start one point back from start_idx so we check all relevant segments */
			if (reverse) i = (start_idx < p.size() - 1) ? start_idx + 1 : 0;
			else i = (start_idx > 0) ? start_idx - 1 : p.size() - 1;
			do {
				i = find_best_travel_point(b, bound_idx, p0, i, end_idx, reverse);
				size_t path_pt_idx, path_idx = find_nearest_path(paths, p[i].X, p[i].Y, NULL, &path_pt_idx);
				comb_dist += append_comb_move(m, island, comb_moves, p0, paths[path_idx][path_pt_idx], &force_retract);
				p0 = paths[path_idx][path_pt_idx];
				if (!crosses_boundary_2pt(p, p0, p1, NULL) && distance_to_point(p0, p1) < closest_dist)
					break;  /* Progress was made and there are no more boundary crossings for current boundary */
			} while (i != end_idx);
		}
		const fl_t dist = distance_to_point(p0, p1);
		if (dist >= closest_dist) {
			b.erase(b.begin() + bound_idx);
			last_bound_idx = -1;
			force_retract = true;
			DEBUG("combed_travel(): warning: ended up farther away at z = %f\n", CINT_TO_FL_T(m->z));
		}
		else {
			closest_dist = dist;
			last_bound_idx = bound_idx;
		}
	}
	comb_dist += distance_to_point(p0, p1) / config.scale_constant;
	if (force_retract || comb_dist >= retract_threshold)
		do_retract(slice, m);
	for (ClipperLib::IntPoint &pt : comb_moves)
		append_linear_travel(slice, m, pt.X, pt.Y, m->z, feed_rate);
}

/* Move to the point nearest to the target location */
static void move_to_island_exit(struct slice *slice, struct machine *m, ClipperLib::cInt x, ClipperLib::cInt y, fl_t feed_rate)
{
	size_t path_pt_idx, path_idx = find_nearest_path(slice->last_comb_paths, x, y, NULL, &path_pt_idx);
	const ClipperLib::IntPoint &point = slice->last_comb_paths[path_idx][path_pt_idx];
	combed_travel(slice, NULL, m, slice->last_boundaries, slice->last_comb_paths, point.X, point.Y, feed_rate, 0.0);
	append_linear_travel(slice, m, point.X, point.Y, m->z, feed_rate);
}

static void linear_move(struct slice *slice, const struct island *island, struct machine *m, ClipperLib::cInt x, ClipperLib::cInt y, ClipperLib::cInt z, fl_t extra_e_len, fl_t feed_rate, fl_t flow_adjust, bool scalable, bool is_travel, bool doing_infill)
{
	const fl_t f_x = CINT_TO_FL_T(x), f_y = CINT_TO_FL_T(y), f_z = CINT_TO_FL_T(z);
	const fl_t f_mx = CINT_TO_FL_T(m->x), f_my = CINT_TO_FL_T(m->y), f_mz = CINT_TO_FL_T(m->z);
	struct g_move move = { x, y, z, 0.0, feed_rate, scalable, is_travel, false };
	const fl_t len = sqrt((f_mx - f_x) * (f_mx - f_x) + (f_my - f_y) * (f_my - f_y) + (f_mz - f_z) * (f_mz - f_z));
	if (is_travel) {
		fl_t retract_threshold = ((doing_infill && !config.retract_within_island) ? config.retract_threshold : config.retract_min_travel);
		if (m->force_retract)
			do_retract(slice, m);
		if (z == m->z && config.comb) {
			if (slice->last_boundaries.size() > 0) {
				/* Inside an island and moving to a point outside of it */
				do_retract(slice, m);
				if (slice->last_comb_paths.size() > 0)
					move_to_island_exit(slice, m, x, y, feed_rate);
				slice->last_boundaries.clear();
				slice->last_comb_paths.clear();
				combed_travel(slice, island, m, slice->printed_outer_boundaries, slice->printed_outer_comb_paths, x, y, feed_rate, retract_threshold);
			}
			else if (island) {
				/* Moving within an island */
				combed_travel(slice, island, m, island->boundaries, island->comb_paths, x, y, feed_rate, retract_threshold);
			}
			else {
				/* Moving between two points that are not in an island */
				combed_travel(slice, island, m, slice->printed_outer_boundaries, slice->printed_outer_comb_paths, x, y, feed_rate, retract_threshold);
			}
		}
		else if (!m->is_retracted && config.retract_len > 0.0
			&& (slice->last_boundaries.size() > 0
				|| len > retract_threshold
				|| (island && crosses_boundary(m, island->boundaries, x, y) >= 0)
				|| (island && len > config.extrusion_width * 2.0 && crosses_exposed_surface(island, m->x, m->y, x, y)))) {
			do_retract(slice, m);
		}
	}
	else {
		if (m->is_retracted && config.retract_len > 0.0) {
			struct g_move restart_move = { m->x, m->y, m->z, config.retract_len, config.restart_speed, false, false, true };
			if (config.extra_restart_len < 0.0)
				restart_move.e += config.extra_restart_len;
			else
				extra_e_len += config.extra_restart_len;
			append_g_move(slice, restart_move, restart_move.e);
			m->is_retracted = false;
		}
		move.e = len * config.extrusion_area * config.flow_multiplier * flow_adjust / config.material_area;
	}
	if (extra_e_len != 0.0) {
		struct g_move restart_move = { m->x, m->y, m->z, extra_e_len, feed_rate * config.extrusion_area / config.material_area, true, false, true };
		append_g_move(slice, restart_move, fabs(extra_e_len));
	}
	if (x != m->x || y != m->y || z != m->z || move.e != 0.0) {
		append_g_move(slice, move, len);
		m->x = x;
		m->y = y;
		m->z = z;
	}
	m->force_retract = false;
}

static bool path_len_is_greater_than(const ClipperLib::Path &p, fl_t len)
{
	fl_t l = 0.0;
	fl_t x0 = CINT_TO_FL_T(p[0].X), y0 = CINT_TO_FL_T(p[0].Y);
	for (size_t i = 1; i < p.size(); ++i) {
		const fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		if (l > len)
			return true;
		x0 = x1;
		y0 = y1;
	}
	const fl_t x1 = CINT_TO_FL_T(p[0].X), y1 = CINT_TO_FL_T(p[0].Y);
	l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
	if (l > len)
		return true;
	return false;
}

/* Requires the last point to be the actual end point. For a closed
   path, the last point should be the same as the first point. */
static void clip_path_from_end(ClipperLib::Path &p, ClipperLib::Path *clipped_points, fl_t clip)
{
	fl_t x0 = CINT_TO_FL_T((p.end() - 1)->X), y0 = CINT_TO_FL_T((p.end() - 1)->Y), l = 0.0;
	for (;;) {
		if (clipped_points)
			clipped_points->push_back(*(p.end() - 1));
		p.erase(p.end() - 1);
		const fl_t x1 = CINT_TO_FL_T((p.end() - 1)->X), y1 = CINT_TO_FL_T((p.end() - 1)->Y);
		const fl_t xv = x1 - x0, yv = y1 - y0;
		const fl_t norm = sqrt(xv * xv + yv * yv);
		l += norm;
		if (l == clip)
			break;
		else if (l > clip) {
			const fl_t new_x = x1 - (l - clip) * (xv / norm), new_y = y1 - (l - clip) * (yv / norm);
			p.push_back(FL_T_TO_INTPOINT(new_x, new_y));
			break;
		}
		x0 = x1;
		y0 = y1;
	}
	if (clipped_points)
		std::reverse(clipped_points->begin(), clipped_points->end());
}

/* Do a non-stationary retract at the end of a shell. Expects a standard closed path where the first point is also the end point (in other words, no repeated points). */
static size_t moving_retract(const ClipperLib::Path &p, struct slice *slice, struct machine *m, ClipperLib::cInt z, size_t start_idx, fl_t feed_rate)
{
	const fl_t len_ratio = config.moving_retract_speed / feed_rate;
	const fl_t move_len = config.retract_len / len_ratio;
	fl_t x0 = CINT_TO_FL_T(m->x), y0 = CINT_TO_FL_T(m->y), l = 0.0, rl = 0.0;
	size_t i = start_idx;
	for (;; ++i) {
		if (i >= p.size())
			i = 0;
		const fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		const fl_t xv = x1 - x0, yv = y1 - y0;
		const fl_t norm = sqrt(xv * xv + yv * yv);
		l += norm;
		if (rl + norm * len_ratio >= config.retract_len) {
			const fl_t new_x = x1 - (l - move_len) * (xv / norm), new_y = y1 - (l - move_len) * (yv / norm);
			const struct g_move move = { FL_T_TO_CINT(new_x), FL_T_TO_CINT(new_y), z, -(config.retract_len - rl), feed_rate, false, false, false };
			append_g_move(slice, move, move_len - (l - norm));
			m->x = move.x;
			m->y = move.y;
			m->z = move.z;
			break;
		}
		else if (norm > 0.0) {
			const struct g_move move = { p[i].X, p[i].Y, z, -norm * len_ratio, feed_rate, false, false, false };
			append_g_move(slice, move, norm);
			/* No need to update m->{x,y,z} */
		}
		rl += norm * len_ratio;
		x0 = x1;
		y0 = y1;
	}
	m->is_retracted = true;
	return i;
}

/* Wipe the nozzle at the end of a shell */
static void shell_wipe(const ClipperLib::Path &p, struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, size_t start_idx, fl_t feed_rate)
{
	fl_t x0 = CINT_TO_FL_T(m->x), y0 = CINT_TO_FL_T(m->y), l = 0.0;
	for (size_t i = start_idx;; ++i) {
		if (i >= p.size())
			i = 0;
		const fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		const fl_t xv = x1 - x0, yv = y1 - y0;
		const fl_t norm = sqrt(xv * xv + yv * yv);
		l += norm;
		if (l >= config.wipe_len) {
			const fl_t new_x = x1 - (l - config.wipe_len) * (xv / norm), new_y = y1 - (l - config.wipe_len) * (yv / norm);
			linear_move(slice, island, m, FL_T_TO_CINT(new_x), FL_T_TO_CINT(new_y), z, 0.0, feed_rate, 1.0, false, true, false);
			break;
		}
		else if (norm > 0.0) {
			linear_move(slice, island, m, p[i].X, p[i].Y, z, 0.0, feed_rate, 1.0, false, true, false);
		}
		x0 = x1;
		y0 = y1;
	}
}

static void generate_closed_path_moves(const ClipperLib::Path &p, size_t start_idx, struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, fl_t feed_rate)
{
	if (p.size() < 3)
		return;
	fl_t total_clip = 0.0;
	bool first_point = true, do_anchor = false;
	if (config.shell_clip > 0.0 && path_len_is_greater_than(p, config.shell_clip * config.extrusion_width * 2.0))
		total_clip += config.shell_clip * config.extrusion_width;
	if (config.anchor && path_len_is_greater_than(p, total_clip + config.extrusion_width)) {
		do_anchor = true;
		total_clip += config.extrusion_width / 2.0 * M_PI_4;
	}
	ClipperLib::Path lp = p;
	if (start_idx != 0)
		std::rotate(lp.begin(), lp.begin() + start_idx, lp.end());
	lp.push_back(lp[0]);
	if (total_clip > 0.0)
		clip_path_from_end(lp, NULL, total_clip);
	ClipperLib::Path coast_path;
	if (config.coast_len > 0.0 && path_len_is_greater_than(p, total_clip + config.coast_len * 2.0))
		clip_path_from_end(lp, &coast_path, config.coast_len);
	for (const ClipperLib::IntPoint &point : lp) {
		if (first_point) {
			linear_move(slice, island, m, point.X, point.Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, false);
			first_point = false;
		}
		else {
			fl_t anchor_e_len = 0.0;
			if (do_anchor) {
				anchor_e_len = config.extrusion_width / 2.0 * M_PI_4 * config.extrusion_area * config.flow_multiplier / config.material_area;
				do_anchor = false;
			}
			linear_move(slice, island, m, point.X, point.Y, z, anchor_e_len, feed_rate, 1.0, true, false, false);
		}
	}
	m->is_retracted = true;  /* Make sure we don't retract */
	for (const ClipperLib::IntPoint &point : coast_path)
		linear_move(slice, island, m, point.X, point.Y, z, 0.0, feed_rate, 1.0, true, true, false);
	m->is_retracted = false;
	if (config.moving_retract && config.retract_len > 0.0)
		start_idx = moving_retract(p, slice, m, z, start_idx, feed_rate);
	if (config.wipe_len > 0.0) {
		m->force_retract = true;
		shell_wipe(p, slice, island, m, z, start_idx, feed_rate);
	}
}

static void plan_brim(struct object *o, struct machine *m, ClipperLib::cInt z)
{
	for (ClipperLib::Paths &p : o->brim) {
		while (!p.empty()) {
			size_t best = 0, start = 0;
			best = find_nearest_path(p, m->x, m->y, NULL, &start);
			generate_closed_path_moves(p[best], start, &o->slices[0], NULL, m, z, config.perimeter_feed_rate);
			p.erase(p.begin() + best);
		}
	}
	m->force_retract = true;
}

static void do_support_wipe(struct slice *slice, ClipperLib::Path &last_line, struct machine *m, ClipperLib::cInt z)
{
	if (config.support_wipe_len > 0.0) {
		m->force_retract = true;
		const fl_t xv = last_line[1].X - last_line[0].X;
		const fl_t yv = last_line[1].Y - last_line[0].Y;
		const fl_t norm = sqrt(xv * xv + yv * yv);
		if (norm > config.support_wipe_len * config.scale_constant)
			linear_move(slice, NULL, m, last_line[1].X - config.support_wipe_len * config.scale_constant * (xv / norm), last_line[1].Y - config.support_wipe_len * config.scale_constant * (yv / norm), z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
		else
			linear_move(slice, NULL, m, last_line[0].X, last_line[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
	}
}

static void plan_support(struct slice *slice, ClipperLib::Paths &lines, struct machine *m, ClipperLib::cInt z, fl_t min_len, fl_t connect_threshold, fl_t flow_adjust, fl_t feed_rate)
{
	ClipperLib::Path last_line(2);
	bool first = true;
	while (!lines.empty()) {
		bool flip_points;
		fl_t best_dist;
		size_t best = find_nearest_segment(lines, m->x, m->y, &best_dist, &flip_points);
		ClipperLib::Path &p = lines[best];
		const fl_t x0 = CINT_TO_FL_T(p[0].X), y0 = CINT_TO_FL_T(p[0].Y);
		const fl_t x1 = CINT_TO_FL_T(p[1].X), y1 = CINT_TO_FL_T(p[1].Y);
		const fl_t len = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		if (len > min_len) {
			bool cross_bound = false;
			if (!first) {
				const ClipperLib::IntPoint p0(m->x, m->y);
				for (const struct island &island : slice->islands) {
					for (const ClipperLib::Path &bound : island.outer_boundaries) {
						if (get_boundary_crossing(bound, last_line[0], last_line[1]) >= 0 || get_boundary_crossing(bound, p0, (flip_points) ? p[1] : p[0]) >= 0) {
							cross_bound = true;
							m->force_retract = true;
							break;
						}
					}
				}
			}
			bool connect = (!first && !cross_bound && best_dist < connect_threshold);
			if (!first && cross_bound)
				do_support_wipe(slice, last_line, m, z);
			if (flip_points)
				std::swap(p[0], p[1]);
			if (connect)
				linear_move(slice, NULL, m, p[0].X, p[0].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
			else
				linear_move(slice, NULL, m, p[0].X, p[0].Y, z, 0.0, config.travel_feed_rate, flow_adjust, false, true, true);
			linear_move(slice, NULL, m, p[1].X, p[1].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
			last_line = p;
			first = false;
		}
		lines.erase(lines.begin() + best);
	}
	if (!first)
		do_support_wipe(slice, last_line, m, z);
}

static void plan_insets_weighted(struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, bool outside_first)
{
	for (;;) {
		bool done = true;
		fl_t best_dist = FL_T_INF;
		size_t best = 0, inset = 0, start = 0;
		for (int i = 0; i < config.shells; ++i) {
			if (!island->insets[i].empty()) {
				fl_t dist;
				size_t r, start_tmp = 0;
				if (config.align_seams && (config.align_interior_seams || i == 0))
					r = find_nearest_aligned_path(island->insets[i], m->x, m->y, &dist);
				else
					r = find_nearest_path(island->insets[i], m->x, m->y, &dist, &start_tmp);
				if (outside_first) {
					if (i != 0)
						dist = dist * (i + 1) + config.retract_min_travel;  /* prefer exterior */
				}
				else {
					if (i != config.shells - 1)
						dist = dist * (config.shells - i) + config.retract_min_travel;  /* prefer interior */
				}
				if (dist < best_dist) {
					best_dist = dist;
					best = r;
					inset = i;
					start = start_tmp;
					done = false;
				}
			}
		}
		if (done)
			break;
		generate_closed_path_moves(island->insets[inset][best], start, slice, island, m, z, (inset == 0) ? config.perimeter_feed_rate : config.loop_feed_rate);
		island->insets[inset].erase(island->insets[inset].begin() + best);
	}
}

static void plan_insets_strict_order(struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, bool outside_first)
{
	int i = (outside_first) ? 0 : config.shells - 1;
	while (i >= 0 && i < config.shells) {
		if (island->insets[i].empty()) {
			i = (outside_first) ? i + 1 : i - 1;
			continue;
		}
		size_t best, start = 0;
		if (config.align_seams && (config.align_interior_seams || i == 0))
			best = find_nearest_aligned_path(island->insets[i], m->x, m->y, NULL);
		else
			best = find_nearest_path(island->insets[i], m->x, m->y, NULL, &start);
		generate_closed_path_moves(island->insets[i][best], start, slice, island, m, z, (i == 0) ? config.perimeter_feed_rate : config.loop_feed_rate);
		island->insets[i].erase(island->insets[i].begin() + best);
	}
}

static void plan_insets(struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, bool outside_first)
{
	if (config.shells > 0) {
		if (config.strict_shell_order)
			plan_insets_strict_order(slice, island, m, z, outside_first);
		else
			plan_insets_weighted(slice, island, m, z, outside_first);
	}
	if (config.retract_after_shells)
		m->force_retract = true;
}

static void plan_infill_simple(ClipperLib::Paths &lines, struct slice *slice, struct island *island, struct machine *m, fl_t feed_rate, fl_t flow_adjust, ClipperLib::cInt z)
{
	while (!lines.empty()) {
		bool flip_points;
		const size_t best = find_nearest_segment(lines, m->x, m->y, NULL, &flip_points);
		ClipperLib::Path &p = lines[best];
		if (flip_points)
			std::swap(p[0], p[1]);
		linear_move(slice, island, m, p[0].X, p[0].Y, z, 0.0, config.travel_feed_rate, flow_adjust, false, true, true);
		linear_move(slice, island, m, p[1].X, p[1].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
		lines.erase(lines.begin() + best);
	}
}

static size_t find_next_solid_infill_segment(const ClipperLib::Paths &p, ClipperLib::Path &line0, fl_t *r_dist, bool *r_flip, bool *r_is_adjacent)
{
	const fl_t adjacent_dist_fudge = config.extrusion_width / 8.0;
	bool best_flip = false, best_is_adjacent = false;
	size_t best = 0;
	fl_t best_dist = FL_T_INF, best_adj_dist = FL_T_INF;

	for (size_t i = 0; i < p.size(); ++i) {
		if (p[i].size() > 2)
			fprintf(stderr, "error: bug in clipper: line segment has more than two points!\n");
		ClipperLib::Path line1 = p[i];
		const fl_t l_dist0 = distance_to_line(line0[0], line1[0], line1[1]);
		const fl_t l_dist1 = distance_to_line(line0[1], line1[0], line1[1]);
		const fl_t l_dist2 = distance_to_line(line1[0], line0[0], line0[1]);
		const fl_t l_dist3 = distance_to_line(line1[1], line0[0], line0[1]);
		const fl_t min_dist = MINIMUM_4(l_dist0, l_dist1, l_dist2, l_dist3);
		const fl_t scaled_min_dist = min_dist / config.scale_constant;
		const fl_t scaled_p_dist = perpendicular_distance_to_line(line0[1], line1[0], line1[1]) / config.scale_constant;
		const fl_t pt_dist0 = distance_to_point(line0[1], line1[0]);
		const fl_t pt_dist1 = distance_to_point(line0[1], line1[1]);
		const bool is_adjacent = scaled_p_dist < config.extrusion_width + adjacent_dist_fudge && scaled_p_dist > config.extrusion_width - adjacent_dist_fudge && scaled_min_dist < config.extrusion_width * 2.0;
		if (pt_dist0 > pt_dist1)
			std::swap(line1[0], line1[1]);
		const bool is_opposite_dir = ((line0[0].X < line0[1].X) != (line1[0].X < line1[1].X)) || ((line0[0].Y < line0[1].Y) != (line1[0].Y < line1[1].Y));
		fl_t adj_dist = l_dist1;
		if (!is_opposite_dir)
			adj_dist *= 2.0;
		if (!is_adjacent)
			adj_dist *= 2.0;
		if (adj_dist < best_adj_dist) {
			best_adj_dist = adj_dist;
			best_flip = (is_adjacent && !is_opposite_dir) || (pt_dist0 > pt_dist1);
			best_dist = (best_flip) ? pt_dist1 : pt_dist0;
			best_is_adjacent = is_adjacent;
			best = i;
		}
	}
	if (r_dist)
		*r_dist = best_dist / config.scale_constant;
	if (r_flip)
		*r_flip = best_flip;
	if (r_is_adjacent)
		*r_is_adjacent = best_is_adjacent;
	return best;
}

static void plan_smoothed_solid_infill(ClipperLib::Paths &lines, struct slice *slice, struct island *island, struct machine *m, fl_t feed_rate, ClipperLib::cInt z)
{
	if (lines.empty())
		return;
	bool flip_points, last_was_smoothed = false, needs_travel = true;
	size_t best = find_nearest_segment(lines, m->x, m->y, NULL, &flip_points);
	ClipperLib::Path line0 = lines[best];
	lines.erase(lines.begin() + best);
	if (flip_points)
		std::swap(line0[0], line0[1]);
	while (!lines.empty()) {
		fl_t best_dist;
		bool is_adjacent;
		best = find_next_solid_infill_segment(lines, line0, &best_dist, &flip_points, &is_adjacent);
		ClipperLib::Path line1 = lines[best];
		lines.erase(lines.begin() + best);
		if (flip_points)
			std::swap(line1[0], line1[1]);
		bool cross_bound = false;
		for (const ClipperLib::Path &bound : island->solid_infill_boundaries) {
			   if (get_boundary_crossing(bound, line0[1], line1[0]) >= 0) {
					   cross_bound = true;
					   break;
			   }
		}
		bool is_constrained = false, in_outer = false, in_hole = false;
		if (island->constraining_edge.empty()) {
			is_constrained = true;  /* All solid infill is inset gap fill if there's no constraining edge */
		}
		else {
			for (const ClipperLib::Path &bound : island->constraining_edge) {
				/* Note: is_constrained will always be true for inset gap fill */
				bool in_bound = ClipperLib::PointInPolygon(line0[1], bound) || ClipperLib::PointInPolygon(line1[0], bound);
				bool bound_is_hole = !ClipperLib::Orientation(bound);
				if (in_bound == bound_is_hole) {
					is_constrained = true;
					in_hole = in_hole || bound_is_hole;
				}
				else if (in_bound && !bound_is_hole) {
					in_outer = true;
				}
			}
			if (is_constrained && in_outer && !in_hole)
				is_constrained = false;
		}
		const ClipperLib::IntPoint line0_midpoint((line0[0].X + line0[1].X) / 2, (line0[0].Y + line0[1].Y) / 2);
		const ClipperLib::IntPoint line1_midpoint((line1[0].X + line1[1].X) / 2, (line1[0].Y + line1[1].Y) / 2);
		const fl_t len_line0 = distance_to_point(line0[0], line0[1]);
		const fl_t len_line1 = distance_to_point(line1[0], line1[1]);
		const fl_t len_mid = distance_to_point(line0_midpoint, line1_midpoint);
		const fl_t xv0 = line0[1].X - line0[0].X, yv0 = line0[1].Y - line0[0].Y;
		const fl_t xv1 = line1[1].X - line1[0].X, yv1 = line1[1].Y - line1[0].Y;
		const fl_t xv_mid = line1_midpoint.X - line0_midpoint.X, yv_mid = line1_midpoint.Y - line0_midpoint.Y;
		/* det(a, b) / (||a|| * ||b||) * ||a|| gives the width of the region covered by a (measured perpendicular to b) */
		const fl_t region_width0 = fabs((xv0 * yv_mid - yv0 * xv_mid) / len_mid);
		const fl_t region_width1 = fabs((xv1 * yv_mid - yv1 * xv_mid) / len_mid);
		const fl_t p_dist = perpendicular_distance_to_line(line0[1], line1[0], line1[1]) / config.scale_constant;
		const fl_t shortening_dist = best_dist / p_dist * config.extrusion_width / 2.0;
		const bool is_opposite_dir = ((line0[0].X < line0[1].X) != (line1[0].X < line1[1].X)) || ((line0[0].Y < line0[1].Y) != (line1[0].Y < line1[1].Y));
		const fl_t connect_min_len = MAXIMUM(shortening_dist, config.extrusion_width / 2.0) * config.scale_constant;
		if (config.infill_smooth_threshold > 0.0
				&& !cross_bound
				&& is_adjacent
				&& len_line0 <= config.extrusion_width * config.infill_smooth_threshold * 2.0 * config.scale_constant
				&& len_line1 <= config.extrusion_width * config.infill_smooth_threshold * 2.0 * config.scale_constant
				&& region_width0 <= config.extrusion_width * config.infill_smooth_threshold * config.scale_constant
				&& region_width1 <= config.extrusion_width * config.infill_smooth_threshold * config.scale_constant) {
			if (!last_was_smoothed) {
				/* move to line0 start */
				if (needs_travel)
					linear_move(slice, island, m, line0[0].X, line0[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
				/* extrude to line0_midpoint */
				linear_move(slice, island, m, line0_midpoint.X, line0_midpoint.Y, z, 0.0, feed_rate, 1.0, true, false, true);
			}
			/* extrude to line1_midpoint */
			const fl_t extrude_ratio = (len_line0 + len_line1) / 2.0 / distance_to_point(line0_midpoint, line1_midpoint);
			const fl_t scaled_feed_rate = (feed_rate / extrude_ratio < config.travel_feed_rate) ? feed_rate / extrude_ratio : config.travel_feed_rate;
			linear_move(slice, island, m, line1_midpoint.X, line1_midpoint.Y, z, 0.0, scaled_feed_rate, extrude_ratio, true, false, true);
			last_was_smoothed = true;
			needs_travel = false;
		}
		else if (!cross_bound
				&& !is_constrained
				&& is_adjacent
				&& is_opposite_dir
				&& best_dist < config.extrusion_width * 3.864
				&& ((last_was_smoothed) ? len_line0 / 2.0 : len_line0) > connect_min_len
				&& len_line1 / 2.0 > connect_min_len) {
			ClipperLib::IntPoint pt0, pt1;
			/* adjusted line0 end point */
			pt0.X = line0[1].X - llround(shortening_dist * config.scale_constant * (xv0 / len_line0));
			pt0.Y = line0[1].Y - llround(shortening_dist * config.scale_constant * (yv0 / len_line0));
			/* adjusted line1 start point */
			pt1.X = line1[0].X - llround(shortening_dist * config.scale_constant * (-xv1 / len_line1));
			pt1.Y = line1[0].Y - llround(shortening_dist * config.scale_constant * (-yv1 / len_line1));
			if (needs_travel)
				linear_move(slice, island, m, line0[0].X, line0[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
			/* extrude line0 */
			linear_move(slice, island, m, pt0.X, pt0.Y, z, 0.0, feed_rate, 1.0, true, false, true);
			/* extrude connection */
			linear_move(slice, island, m, pt1.X, pt1.Y, z, 0.0, feed_rate, 1.0, true, false, true);
			last_was_smoothed = false;
			needs_travel = false;
		}
		else {
			if (needs_travel)
				linear_move(slice, island, m, line0[0].X, line0[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
			linear_move(slice, island, m, line0[1].X, line0[1].Y, z, 0.0, feed_rate, 1.0, true, false, true);
			last_was_smoothed = false;
			needs_travel = true;
		}
		line0 = line1;
	}
	if (needs_travel)
		linear_move(slice, island, m, line0[0].X, line0[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
	linear_move(slice, island, m, line0[1].X, line0[1].Y, z, 0.0, feed_rate, 1.0, true, false, true);
}

static void plan_moves(struct object *o, struct slice *slice, ssize_t layer_num, struct machine *m)
{
	const ClipperLib::cInt z = FL_T_TO_CINT(((fl_t) layer_num) * config.layer_height + config.layer_height + config.object_z_extra);
	if (layer_num == 0 && config.brim_lines > 0)
		plan_brim(o, m, z);
	if (config.generate_support) {
		const fl_t support_flow_adjust = (layer_num > 0) ? config.support_flow_mult : 1.0;
		const fl_t support_feed_rate = (layer_num > 0) ? config.support_feed_rate : config.perimeter_feed_rate;
		plan_support(slice, slice->support_interface_lines, m, z, config.extrusion_width, (layer_num == 0 || config.connect_support_lines) ? (layer_num == 0 && config.solid_support_base) ? config.extrusion_width * 1.9 : config.extrusion_width / config.interface_density * 1.9 : 0.0, support_flow_adjust, support_feed_rate);
		plan_support(slice, slice->support_lines, m, z, config.extrusion_width * 2.0, (layer_num == 0 || config.connect_support_lines) ? config.extrusion_width / config.support_density * 10.0 : 0.0, support_flow_adjust, support_feed_rate);
	}
	while (slice->islands.size() > 0) {
		size_t best = 0;
		fl_t best_dist = FL_T_INF;
		for (size_t i = 0; i < slice->islands.size(); ++i) {
			fl_t dist;
			if (config.align_seams)
				find_nearest_aligned_path(slice->islands[i].insets[0], m->x, m->y, &dist);
			else
				find_nearest_path(slice->islands[i].insets[0], m->x, m->y, &dist, NULL);
			if (dist < best_dist) {
				best = i;
				best_dist = dist;
			}
		}
		struct island &island = slice->islands[best];
		plan_insets(slice, &island, m, z, config.outside_first || layer_num == 0);
		plan_smoothed_solid_infill(island.solid_infill, slice, &island, m, config.solid_infill_feed_rate, z);
		plan_infill_simple(island.iron_paths, slice, &island, m, config.iron_feed_rate, config.iron_flow_multiplier, z);
		plan_infill_simple(island.sparse_infill, slice, &island, m, config.sparse_infill_feed_rate, 1.0, z);
		delete[] island.insets;
		delete[] island.inset_gaps;
		if (config.comb) {
			/* Insert outer boundaries and comb paths for the island we just printed */
			slice->printed_outer_boundaries.insert(slice->printed_outer_boundaries.end(), island.outer_boundaries.begin(), island.outer_boundaries.end());
			slice->printed_outer_comb_paths.insert(slice->printed_outer_comb_paths.end(), island.outer_comb_paths.begin(), island.outer_comb_paths.end());
			/* Set last_boundaries and last_comb_paths */
			slice->last_boundaries = island.boundaries;
			slice->last_comb_paths = island.comb_paths;
		}
		slice->islands.erase(slice->islands.begin() + best);
	}
	m->force_retract = true;  /* Force retract on layer change */
	if (config.comb) {
		FREE_VECTOR(slice->last_boundaries);
		FREE_VECTOR(slice->last_comb_paths);
		FREE_VECTOR(slice->printed_outer_boundaries);
		FREE_VECTOR(slice->printed_outer_comb_paths);
	}
}

static void plan_raft(struct object *o, struct slice *slice, struct machine *m)
{
	ClipperLib::cInt z = FL_T_TO_CINT(config.raft_base_layer_height);

	fl_t flow_adjust = (config.raft_base_layer_height * config.raft_base_layer_width) / (config.layer_height * config.extrusion_width);
	fl_t feed_rate = config.solid_infill_feed_rate * config.first_layer_mult;
	plan_support(slice, o->raft[0], m, z, config.extrusion_width * 2.0, config.raft_base_layer_width / config.raft_base_layer_density * 1.9, flow_adjust, feed_rate);

	flow_adjust = config.raft_interface_flow_mult;
	feed_rate = config.solid_infill_feed_rate;
	for (int i = 1; i <= config.raft_interface_layers; ++i) {
		ClipperLib::Paths lines = o->raft[1];
		z = FL_T_TO_CINT(config.raft_base_layer_height + config.layer_height * i);
		plan_support(slice, lines, m, z, config.extrusion_width * 2.0, config.extrusion_width * 1.9, flow_adjust, feed_rate);
	}
	m->force_retract = true;
}

static void write_gcode_move(FILE *f, const struct g_move *move, struct machine *m, fl_t feed_rate_mult, bool force_xyz)
{
	fl_t feed_rate = move->feed_rate;
	if (move->scalable) {
		feed_rate *= feed_rate_mult;
		const fl_t min_feed_rate = (move->is_restart) ? config.min_feed_rate * config.extrusion_area / config.material_area : config.min_feed_rate;
		if (feed_rate < min_feed_rate)
			feed_rate = min_feed_rate;
	}
	if (move->is_travel && move->z != m->z && config.separate_z_travel) {
		fprintf(f, "G1 Z%.3f", CINT_TO_FL_T(move->z));
		if (feed_rate != m->feed_rate)
			fprintf(f, " F%ld", (feed_rate * 60.0 <= 1.0) ? 1 : lround(feed_rate * 60.0));
		fputc('\n', f);
		m->z = move->z;
	}
	fputs("G1", f);
	if (force_xyz || move->x != m->x)
		fprintf(f, " X%.3f", CINT_TO_FL_T(move->x));
	if (force_xyz || move->y != m->y)
		fprintf(f, " Y%.3f", CINT_TO_FL_T(move->y));
	if (force_xyz || move->z != m->z)
		fprintf(f, " Z%.3f", CINT_TO_FL_T(move->z));
	if (move->e != 0.0)
		fprintf(f, " E%.5f", m->e + move->e);
	if (feed_rate != m->feed_rate)
		fprintf(f, " F%ld", (feed_rate * 60.0 <= 1.0) ? 1 : lround(feed_rate * 60.0));
	fputc('\n', f);
	m->x = move->x;
	m->y = move->y;
	m->z = move->z;
	m->e += move->e;
	m->feed_rate = feed_rate;
}

#define NEW_PLAN_MACHINE(name, obj) struct machine name = { FL_T_TO_CINT(obj->c.x - (obj->w + config.xy_extra) / 2.0), FL_T_TO_CINT(obj->c.y - (obj->d + config.xy_extra) / 2.0), 0, 0.0, 0.0, true, false }

static int write_gcode(const char *path, struct object *o)
{
	FILE *f;
	if (strcmp(path, "-") == 0)
		f = stdout;
	else
		f = fopen(path, "w");
	if (!f)
		return 1;
	bool is_first_move = true;
#if defined(SHIV_SINGLE_THREADED_PATH_PLANNING)
	NEW_PLAN_MACHINE(plan_m, o);
#endif
	struct machine export_m = {};
	fl_t total_e = 0.0;
	fl_t feed_rate_mult = config.first_layer_mult;
	fprintf(stderr, "plan moves and write gcode to %s...", path);
	write_gcode_string(config.start_gcode, f, false);
	if (config.generate_raft) {
		struct slice raft_dummy_slice;
	#if !defined(SHIV_SINGLE_THREADED_PATH_PLANNING)
		NEW_PLAN_MACHINE(plan_m, o);
	#endif
		plan_raft(o, &raft_dummy_slice, &plan_m);
		linear_move(&raft_dummy_slice, NULL, &plan_m, plan_m.x, plan_m.y, plan_m.z, 0.0, config.travel_feed_rate, 1.0, false, true, false);  /* do retract, if needed */
		fputs("; raft\n", f);
		for (const struct g_move &move : raft_dummy_slice.moves) {
			write_gcode_move(f, &move, &export_m, 1.0, is_first_move);
			is_first_move = false;
		}
	}
#if !defined(SHIV_SINGLE_THREADED_PATH_PLANNING) && defined(_OPENMP)
	#pragma omp parallel for ordered schedule(dynamic)
#endif
	for (ssize_t i = 0; i < o->n_slices; ++i) {
		struct slice *slice = &o->slices[i];
	#if !defined(SHIV_SINGLE_THREADED_PATH_PLANNING)
		NEW_PLAN_MACHINE(plan_m, o);
	#endif
		plan_moves(o, slice, i, &plan_m);
		linear_move(slice, NULL, &plan_m, plan_m.x, plan_m.y, plan_m.z, 0.0, config.travel_feed_rate, 1.0, false, true, false);  /* do retract, if needed */
	#if !defined(SHIV_SINGLE_THREADED_PATH_PLANNING) && defined(_OPENMP)
		#pragma omp ordered
		{
	#endif
		fprintf(f, "; layer %zd (z = %f)\n", i, ((fl_t) i) * config.layer_height + config.layer_height + config.object_z_extra);
		for (struct at_layer_gcode &g : config.at_layer)
			if (g.layer == i)
				write_gcode_string(g.value, f, false);
		if (i == config.cool_layer)
			write_gcode_string(config.cool_on_gcode, f, false);
		fl_t average_layer_time = slice->layer_time / feed_rate_mult;
		for (int k = 1; k < config.layer_time_samples; ++k)
			average_layer_time += (k < i) ? o->slices[i - k].layer_time : o->slices[0].layer_time / config.first_layer_mult;
		average_layer_time /= config.layer_time_samples;
		if (average_layer_time < config.min_layer_time)
			feed_rate_mult *= average_layer_time / config.min_layer_time;
		for (const struct g_move &move : slice->moves) {
			write_gcode_move(f, &move, &export_m, feed_rate_mult, is_first_move);
			is_first_move = false;
		}
		feed_rate_mult = 1.0;
		total_e += export_m.e;
		export_m.e = 0.0;
		fputs("G92 E0\n", f);
	#if !defined(SHIV_SINGLE_THREADED_PATH_PLANNING) && defined(_OPENMP)
		}
	#endif
		FREE_VECTOR(slice->moves);
	}
	write_gcode_string(config.cool_off_gcode, f, false);
	write_gcode_string(config.end_gcode, f, false);
	fputs(" done\n", stderr);
	const fl_t mass = config.material_area * total_e * config.material_density / config.flow_multiplier;
	fprintf(f, "; material length = %.4f\n", total_e / config.flow_multiplier);
	fprintf(f, "; material mass   = %.4f\n", mass);
	fprintf(f, "; material cost   = %.4f\n", mass * config.material_cost);
	fprintf(stderr, "material length = %.4f\n", total_e / config.flow_multiplier);
	fprintf(stderr, "material mass   = %.4f\n", mass);
	fprintf(stderr, "material cost   = %.4f\n", mass * config.material_cost);
	const long int bytes = ftell(f);
	if (bytes >= 2048 * 1024)
		fprintf(stderr, "wrote %.2fMiB\n", bytes / 1024.0 / 1024.0);
	else if (bytes >= 2048)
		fprintf(stderr, "wrote %.2fKiB\n", bytes / 1024.0);
	else
		fprintf(stderr, "wrote %ldB\n", bytes);
	fclose(f);
	return 0;
}

#define GET_FEED_RATE(x, m) (((x) >= 0.0) ? (x) : (m) * -(x))

int main(int argc, char *argv[])
{
	int opt, ret;
	char *path, *output_path = NULL;
	struct object o;
	fl_t scale_factor = 1.0, x_translate = 0.0, y_translate = 0.0, z_chop = 0.0;
	bool do_preview = false;

	/* Parse options */
	while ((opt = getopt(argc, argv, ":hpo:c:O:S:l:w:t:s:d:n:r:f:b:C:x:y:z:")) != -1) {
		char *key, *value;
		int ret;
		switch (opt) {
		case 'h':
			fputs(usage_string, stderr);
			return 0;
		case 'p':
			do_preview = true;
			break;
		case 'o':
			output_path = optarg;
			break;
		case 'c':
			ret = read_config(optarg);
			if (ret == 1)
				fprintf(stderr, "error: failed to open config file: %s: %s\n", optarg, strerror(errno));
			if (ret)
				return 1;
			fprintf(stderr, "loaded config file: %s\n", optarg);
			break;
		case 'O':
			fprintf(stderr, "warning: -O is deprecated; please use -S instead\n");
		case 'S':
			key = strdup(optarg);
			value = isolate(key, '=');
			if (set_config_setting(key, value, 0, NULL))
				return 1;
			free(key);
			break;
		case 'l':
			if (set_config_setting("layer_height", optarg, 0, NULL))
				return 1;
			break;
		case 'w':
			if (set_config_setting("extrusion_width", optarg, 0, NULL))
				return 1;
			break;
		case 't':
			if (set_config_setting("tolerance", optarg, 0, NULL))
				return 1;
			break;
		case 's':
			scale_factor = atof(optarg);
			if (scale_factor == 0.0) {
				fputs("error: scale_factor cannot be 0\n", stderr);
				return 1;
			}
			break;
		case 'd':
			if (set_config_setting("infill_density", optarg, 0, NULL))
				return 1;
			break;
		case 'n':
			if (set_config_setting("shells", optarg, 0, NULL))
				return 1;
			break;
		case 'r':
			if (set_config_setting("roof_thickness", optarg, 0, NULL))
				return 1;
			break;
		case 'f':
			if (set_config_setting("floor_thickness", optarg, 0, NULL))
				return 1;
			break;
		case 'b':
			if (set_config_setting("brim_width", optarg, 0, NULL))
				return 1;
			break;
		case 'C':
			if (set_config_setting("coarseness", optarg, 0, NULL))
				return 1;
			break;
		case 'x':
			x_translate = atof(optarg);
			break;
		case 'y':
			y_translate = atof(optarg);
			break;
		case 'z':
			z_chop = atof(optarg);
			break;
		default:
			if (opt == ':')
				fprintf(stderr, "error: expected argument to option '%c'\n", optopt);
			else
				fprintf(stderr, "error: illegal option '%c'\n", optopt);
			return 1;
		}
	}
	if (optind + 1 == argc)
		path = argv[optind];
	else if (optind + 1 < argc) {
		fputs("error: only one input may be given\n", stderr);
		return 1;
	}
	else {
		fputs("error: expected path\n", stderr);
		fputs(usage_string, stderr);
		return 1;
	}

	if (config.layer_height > config.extrusion_width) {
		fputs("error: layer_height must not be greater than extrusion_width\n", stderr);
		return 1;
	}

	config.roof_layers = lround(config.roof_thickness / config.layer_height);
	config.floor_layers = lround(config.floor_thickness / config.layer_height);
	if (config.outside_first || config.shells < 2)
		config.edge_packing_density = 1.0;
	config.extrusion_area = config.extrusion_width * config.layer_height - (config.layer_height * config.layer_height - config.layer_height * config.layer_height * M_PI_4) * (1.0 - config.packing_density);
	config.edge_width = (config.extrusion_area - (config.layer_height * config.layer_height * M_PI_4)) / config.layer_height + config.layer_height;
	config.edge_offset = (config.edge_width + (config.edge_width - config.extrusion_width) * (1.0 - config.edge_packing_density)) / -2.0;
	config.material_area = config.material_diameter * config.material_diameter * M_PI_4;
	if (config.cool_on_gcode == NULL)
		config.cool_on_gcode = strdup(DEFAULT_COOL_ON_STR);
	if (config.cool_off_gcode == NULL)
		config.cool_off_gcode = strdup(DEFAULT_COOL_OFF_STR);
	config.x_center += x_translate;
	config.y_center += y_translate;
	config.brim_lines = lround(config.brim_width / config.extrusion_width);
	config.solid_infill_clip_offset = (0.5 + config.shells - config.fill_threshold - config.min_shell_contact) * config.extrusion_width;
	config.solid_infill_clip_offset = MAXIMUM(config.solid_infill_clip_offset, 0.0);
	config.xy_extra = (config.extra_offset + config.extrusion_width * config.brim_lines) * 2.0;
	if (config.generate_support)
		config.xy_extra += (config.support_xy_expansion + (0.5 + config.support_margin) * config.edge_width - config.edge_offset) * 2.0;
	const fl_t interface_clip_offset_1 = config.extrusion_width * (1.0 - config.edge_overlap) / 2.0 + (0.5 + config.support_margin) * config.edge_width - config.edge_offset - config.extrusion_width / 8.0;
	const fl_t interface_clip_offset_2 = tan(config.support_angle / 180.0 * M_PI) * config.layer_height;
	config.interface_clip_offset = MINIMUM(interface_clip_offset_1, interface_clip_offset_2);
	if (config.generate_raft) {
		config.xy_extra += config.raft_xy_expansion * 2.0;
		config.object_z_extra += config.raft_base_layer_height + config.layer_height * (config.raft_vert_margin + config.raft_interface_layers);
	}
	/* set feed rates */
	config.perimeter_feed_rate = GET_FEED_RATE(config.perimeter_feed_rate, config.feed_rate);
	config.loop_feed_rate = GET_FEED_RATE(config.loop_feed_rate, config.feed_rate);
	config.solid_infill_feed_rate = GET_FEED_RATE(config.solid_infill_feed_rate, config.feed_rate);
	config.sparse_infill_feed_rate = GET_FEED_RATE(config.sparse_infill_feed_rate, config.feed_rate);
	config.support_feed_rate = GET_FEED_RATE(config.support_feed_rate, config.feed_rate);
	config.iron_feed_rate = GET_FEED_RATE(config.iron_feed_rate, config.solid_infill_feed_rate);
	config.travel_feed_rate = GET_FEED_RATE(config.travel_feed_rate, config.feed_rate);
	config.moving_retract_speed = GET_FEED_RATE(config.moving_retract_speed, config.retract_speed);
	config.restart_speed = GET_FEED_RATE(config.restart_speed, config.retract_speed);

	fprintf(stderr, "configuration:\n");
	fprintf(stderr, "  %-24s = %f\n", "scale_factor (-s)", scale_factor);
	for (size_t i = 0; i < LENGTH(settings); ++i) {
		if (settings[i].type != SETTING_TYPE_STR) {
			fprintf(stderr, " %c%-24s = ", (settings[i].read_only) ? '*' : ' ', settings[i].name);
			print_config_setting(stderr, &settings[i], false);
			putc('\n', stderr);
		}
	}
#ifdef _OPENMP
	fprintf(stderr, "OpenMP enabled (%d threads)\n", omp_get_max_threads());
#endif

	fprintf(stderr, "load object...\n");
	memset(&o, 0, sizeof(o));
	ret = read_binary_stl(&o, path);
	if (ret) {
		fprintf(stderr, "error: failed to read stl: %s: %s\n", path, (ret == 2) ? "short read" : strerror(errno));
		return 1;
	}

	fprintf(stderr, "  polygons = %zd\n", o.n);
	fprintf(stderr, "  center   = (%f, %f, %f)\n", o.c.x, o.c.y, o.c.z);
	fprintf(stderr, "  height   = %f\n", o.h);
	fprintf(stderr, "  width    = %f\n", o.w);
	fprintf(stderr, "  depth    = %f\n", o.d);

	fprintf(stderr, "scale and translate object...\n");
	scale_object(&o, config.xy_scale_factor * scale_factor, config.xy_scale_factor * scale_factor, config.z_scale_factor * scale_factor);
	const fl_t z_translate = (config.preserve_layer_offset) ? round((o.h / 2.0 - o.c.z) / config.layer_height) * config.layer_height : o.h / 2.0 - o.c.z;
	translate_object(&o, -o.c.x + config.x_center, -o.c.y + config.y_center, z_translate - z_chop);
	fprintf(stderr, "  center   = (%f, %f, %f)\n", o.c.x, o.c.y, o.c.z);
	fprintf(stderr, "  height   = %f\n", o.h);
	fprintf(stderr, "  width    = %f\n", o.w);
	fprintf(stderr, "  depth    = %f\n", o.d);

	fprintf(stderr, "slice object...\n");
	slice_object(&o);
	if (do_preview)
		preview_slices(&o);
	if (output_path) {
		if (write_gcode(output_path, &o)) {
			fprintf(stderr, "error: failed to write gcode output: %s: %s\n", output_path, strerror(errno));
			return 1;
		}
	}

	return 0;
}
