/*
 * Copyright (C) 2016 Michael Barbour <barbour.michael.0@gmail.com>
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
#include <chrono>
#include <algorithm>
#include <getopt.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "clipper.hpp"
#include "misc_defs.h"
#include "list.h"

#define SCALE_CONSTANT             1000  /* Clipper uses integers, so we need to scale floating point values. Precision is 1 / SCALE_CONSTANT units. */
#define CLEAN_DIST                 (1.41421356 * config.coarseness)
#define USE_BOUNDING_BOX           1
#define SHIV_DEBUG                 1
typedef double fl_t;

#define MINIMUM(a, b) ((a < b) ? a : b)
#define MAXIMUM(a, b) ((a > b) ? a : b)
#define FL_T_TO_CINT(x) ((ClipperLib::cInt) (lround((x) * SCALE_CONSTANT)))
#define CINT_TO_FL_T(x) (((fl_t) (x)) / SCALE_CONSTANT)
#define FL_T_TO_INTPOINT(x, y)  ClipperLib::IntPoint(FL_T_TO_CINT(x), FL_T_TO_CINT(y))
#define INTPOINT_TO_FL_T(point) CINT_TO_FL_T((point).X), CINT_TO_FL_T((point).Y)

#if SHIV_DEBUG
#define DEBUG(...) fprintf(stderr, "DEBUG: " __VA_ARGS__)
#else
#define DEBUG(...) /* no-op */
#endif

static const char e_nomem[] = "fatal: No memory\n";
static const char usage_string[] =
	"usage: shiv [-hp] [-o output_path] [-c config_path] [-O option=value]\n"
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
	"  -O option=value       set option to value\n"
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

/* default config */
#define DEFAULT_COOL_ON_STR  "M106 S255"
#define DEFAULT_COOL_OFF_STR "M107"
static struct {
	fl_t layer_height          = 0.2;
	fl_t tolerance             = 0.001;    /* Segment connection tolerance */
	fl_t coarseness            = 5.0;      /* Approximate output coarseness in 1/SCALE_CONSTANT units (microns with the default SCALE_CONSTANT and mm units) (argument to ClipperLib::CleanPolygon) */
	fl_t extrusion_width       = 0.4;      /* Constrained (solid infill) extrusion width */
	fl_t edge_width;                       /* Unconstrained (edge) extrusion width (calculated from extrusion_width) */
	fl_t extrusion_area;                   /* Cross-sectional area of an extrusion */
	fl_t xy_scale_factor       = 1.003;    /* Scale object in x and y axis by this ratio to compensate for shrinkage */
	fl_t z_scale_factor        = 1.0;      /* Scale object in z axis by this ratio to compensate for shrinkage */
	fl_t x_center              = 0.0;
	fl_t y_center              = 0.0;
	fl_t packing_density       = 0.98;     /* Solid packing density (should be slightly less than 1; 0.98 seems to work well for PLA) */
	fl_t edge_packing_density  = 0.95;     /* Packing density of the contranied half of the outer perimeter */
	fl_t seam_packing_density  = 0.95;     /* Packing density of the ends of each shell (the seam) */
	fl_t extra_offset          = 0.0;      /* Offset the object by this distance in the xy plane */
	fl_t edge_offset;                      /* Offset of the outer perimeter (calculated) */
	fl_t shell_clip;                       /* Shells are clipped by this much (calculated from seam_packing_density) */
	fl_t infill_density        = 0.2;      /* Sparse infill density */
	int shells                 = 2;        /* Number of loops/perimeters/shells (whatever you want to call them) */
	fl_t roof_thickness        = 0.8;      /* Solid surface thickness when looking upwards */
	int roof_layers;                       /* Calculated from roof_thickness */
	fl_t floor_thickness       = 0.8;      /* Solid surface thickness when looking downwards */
	int floor_layers;                      /* Calculated from floor_thickness */
	fl_t min_shell_contact     = 1.0;      /* Minimum contact patch through roof and floor layers in units of extrusion_width. Small values (~0.1 to 1.0) will reduce print time compared to larger values, but may produce weaker parts. */
	fl_t solid_infill_clip_offset;
	fl_t solid_fill_expansion  = 1.0;      /* Distance to expand solid infill, in units of extrusion_width */
	fl_t material_diameter     = 1.75;     /* Diameter of material */
	fl_t material_area;                    /* Cross-sectional area of filament (calculated from material_diameter) */
	fl_t flow_multiplier       = 1.0;      /* Flow rate adjustment */
	fl_t feed_rate             = 50.0;     /* Base feed rate */
	/* Feed rates below are actual speeds if set to a positive value, or a multiple of 'feed_rate' if set to a negative value.
	   In other words, '40' is 40 units/s, but '-0.5' is feed_rate * 0.5 units/s. */
	fl_t perimeter_feed_rate   = -0.5;     /* Outer shell feed rate */
	fl_t loop_feed_rate        = -1.0;     /* Inner shell feed rate */
	fl_t infill_feed_rate      = -1.0;
	fl_t support_feed_rate     = -1.0;
	fl_t travel_feed_rate      = 120.0;
	fl_t first_layer_mult      = 0.5;      /* First layer feed rates (except travel) are multiplied by this value */
	fl_t coast_len             = 0.0;      /* Length to coast (move with the extruder turned off) at the end of a shell */
	fl_t retract_len           = 1.0;
	fl_t retract_speed         = 20.0;
	fl_t restart_speed         = -1.0;     /* A negative value means a multiple of 'retract_speed' */
	fl_t retract_min_travel    = 10.0;     /* Minimum travel for retraction when not crossing a boundary or when printing shells. Has no effect when printing infill if retract_within_island is false. */
	fl_t retract_threshold     = 30.0;     /* Unconditional retraction threshold */
	bool retract_within_island = false;
	fl_t extra_restart_len     = 0.0;      /* Extra material length on restart */
	int cool_layer             = 1;        /* Turn on part cooling at this layer */
	char *start_gcode          = NULL;
	char *end_gcode            = NULL;
	char *cool_on_gcode        = NULL;
	char *cool_off_gcode       = NULL;
	fl_t temp                  = 220.0;    /* Hotend temperature */
	fl_t bed_temp              = 65.0;
	fl_t edge_overlap          = 0.5;      /* Allowable edge path overlap in units of extrusion_width */
	bool strict_shell_order    = false;    /* Always do insets in order within an island */
	bool infill_first          = false;    /* Do infill before shells */
	bool align_seams           = true;     /* Align seams to the lower left corner */
	bool clean_insets          = true;     /* Do ClipperLib::CleanPolygon operation on all insets (only the initial outline is cleaned if this is false) */
	bool fill_inset_gaps       = true;     /* Fill gaps between shells */
	bool no_solid              = false;    /* If true, only generate solid fill on the very top and bottom of the model */
	bool outside_first         = false;    /* Prefer exterior shells */
	bool solid_infill_first    = false;    /* Print solid infill before sparse infill */
	bool separate_z_travel     = false;    /* Generate a separate z travel move instead of moving all axes together */
	bool combine_all           = false;    /* Orients all outlines counter-clockwise. This can be used to fix certain broken models, but it also fills holes. */
	bool generate_support      = false;    /* Generate support structure */
	bool support_everywhere    = false;    /* False means only touching build plate */
	bool solid_support_base    = false;    /* Make supports solid at layer 0 */
	bool connect_support_lines = false;     /* Connect support lines together. Makes the support structure more robust, but harder to remove. */
	ClipperLib::PolyFillType poly_fill_type = ClipperLib::pftNonZero;  /* Set poly fill type for union. Sometimes ClipperLib::pftEvenOdd is useful for broken models with self-intersections and/or incorrect normals. */
	ClipperLib::JoinType inset_join_type    = ClipperLib::jtMiter;     /* Join type for negative offsets */
	ClipperLib::JoinType outset_join_type   = ClipperLib::jtMiter;     /* Join type for positive offsets */
	fl_t offset_miter_limit    = 5.0;
	fl_t offset_arc_tolerance  = 5.0;
	fl_t fill_threshold        = 0.5;      /* Remove infill or inset gap fill when it would be narrower than extrusion_width * fill_threshold */
	fl_t support_angle         = 70.0;     /* Angle threshold for support */
	fl_t support_margin        = 0.6;      /* Horizontal spacing between support and model, in units of edge_width */
	int support_vert_margin    = 1;        /* Vertical spacing between support and model, in layers */
	int interface_layers       = 0;        /* Number of solid support interface layers */
	fl_t support_xy_expansion  = 2.0;      /* Expand support map by this amount. Larger values will generate more support material, but the supports will be stronger. */
	fl_t support_density       = 0.3;      /* Support structure density */
	fl_t support_flow_mult     = 0.75;     /* Flow rate is multiplied by this value for the support structure. Smaller values will generate a weaker support structure, but it will be easier to remove. */
	fl_t min_layer_time        = 8.0;      /* Slow down if the estimated layer time is less than this value */
	int layer_time_samples     = 5;        /* Number of samples in the layer time moving average */
	fl_t min_feed_rate         = 10.0;
	fl_t brim_width            = 0.0;
	int brim_lines;
	fl_t material_density      = 0.00125;  /* Material density in <arbitrary mass unit> / <input/output unit>^3. The default is correct for PLA and millimeter input/output units */
	fl_t material_cost         = 0.01499;  /* Material cost in <arbitrary currency> / <arbitrary mass unit>. The arbitrary mass unit must be the same as used in material_density */

	/* internal stuff */
	fl_t xy_extra              = 0.0;      /* extra xy size (brim, extra_offset, support_xy_expansion, etc...) */
} config;

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
	ClipperLib::Paths sparse_infill_pattern;
	ClipperLib::Paths brim;
	ClipperLib::Paths support_pattern;
	ClipperLib::Paths support_interface_pattern;
};

struct segment {
	struct node n;
	fl_t x[2], y[2];
};

struct cint_rect {
	ClipperLib::cInt x0, y0, x1, y1;
};

struct island {
	ClipperLib::Paths outlines;
	ClipperLib::Paths *insets;
	ClipperLib::Paths *inset_gaps;
	ClipperLib::Paths infill_insets;
	ClipperLib::Paths solid_infill;
	ClipperLib::Paths sparse_infill;
	ClipperLib::Paths boundaries;
	ClipperLib::Paths solid_infill_clip;
	ClipperLib::Paths exposed_surface;
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
	ClipperLib::Paths layer_support_map;
	ClipperLib::Paths support_map;
	ClipperLib::Paths support_boundaries;
	ClipperLib::Paths support_lines;
	ClipperLib::Paths support_interface_lines;
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
	*s = '\0';
	return s + 1;
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

#define CHECK_VALUE(cond, name, strcond) \
	if (!(cond)) { \
 		fputs("error: " name " must be " strcond "\n", stderr); \
		return 1; \
	} \

#define PARSE_BOOL(value) (value[0] == 't' || value[0] == 'T' || value[0] == 'y' || value[0] == 'Y' || atoi(value))

static int parse_join_type(const char *s, ClipperLib::JoinType *jt)
{
	if (strcmp(s, "miter") == 0)
		*jt = ClipperLib::jtMiter;
	else if (strcmp(s, "square") == 0)
		*jt = ClipperLib::jtSquare;
	else if (strcmp(s, "round") == 0)
		*jt = ClipperLib::jtRound;
	else
		return 1;
	return 0;
}

static int set_config_option(const char *key, const char *value, int n, const char *path)
{
	if (strcmp(key, "layer_height") == 0) {
		config.layer_height = atof(value);
		CHECK_VALUE(config.layer_height > 0.0, "layer height", "> 0");
	}
	else if (strcmp(key, "tolerance") == 0) {
		config.tolerance = atof(value);
		CHECK_VALUE(config.tolerance >= 0.0, "tolerance", ">= 0");
	}
	else if (strcmp(key, "coarseness") == 0) {
		config.coarseness = atof(value);
		CHECK_VALUE(config.coarseness >= 0.0, "coarseness", ">= 0");
	}
	else if (strcmp(key, "extrusion_width") == 0) {
		config.extrusion_width = atof(value);
		CHECK_VALUE(config.extrusion_width > 0.0, "extrusion width", "> 0");
	}
	else if (strcmp(key, "xy_scale_factor") == 0) {
		config.xy_scale_factor = atof(value);
		CHECK_VALUE(config.xy_scale_factor > 0.0, "x/y scale factor", "> 0");
	}
	else if (strcmp(key, "z_scale_factor") == 0) {
		config.z_scale_factor = atof(value);
		CHECK_VALUE(config.z_scale_factor > 0.0, "z scale factor", "> 0");
	}
	else if (strcmp(key, "x_center") == 0) {
		config.x_center = atof(value);
	}
	else if (strcmp(key, "y_center") == 0) {
		config.y_center = atof(value);
	}
	else if (strcmp(key, "packing_density") == 0) {
		config.packing_density = atof(value);
		CHECK_VALUE(config.packing_density > 0.0 && config.packing_density <= 1.0, "packing density", "within (0, 1]");
		if (config.packing_density < M_PI / 4)
			fprintf(stderr, "warning: packing density probably shouldn't be < %f\n", M_PI / 4);
	}
	else if (strcmp(key, "edge_packing_density") == 0) {
		config.edge_packing_density = atof(value);
		CHECK_VALUE(config.edge_packing_density > 0.0 && config.edge_packing_density <= 1.0, "edge packing density", "within (0, 1]");
		if (config.edge_packing_density < M_PI / 4)
			fprintf(stderr, "warning: edge packing density probably shouldn't be < %f\n", M_PI / 4);
	}
	else if (strcmp(key, "seam_packing_density") == 0) {
		config.seam_packing_density = atof(value);
		CHECK_VALUE(config.seam_packing_density > 0.0 && config.seam_packing_density <= 1.0, "seam packing density", "within (0, 1]");
		if (config.seam_packing_density < M_PI / 4)
			fprintf(stderr, "warning: seam packing density probably shouldn't be < %f\n", M_PI / 4);
	}
	else if (strcmp(key, "extra_offset") == 0) {
		config.extra_offset = atof(value);
	}
	else if (strcmp(key, "infill_density") == 0) {
		config.infill_density = atof(value);
		CHECK_VALUE(config.infill_density >= 0.0 && config.infill_density <= 1.0, "infill density", "within [0, 1]");
	}
	else if (strcmp(key, "shells") == 0) {
		config.shells = atoi(value);
		CHECK_VALUE(config.shells >= 0, "number of shells", ">= 0");
	}
	else if (strcmp(key, "roof_thickness") == 0) {
		config.roof_thickness = atof(value);
		CHECK_VALUE(config.roof_thickness >= 0.0, "roof thickness", ">= 0");
	}
	else if (strcmp(key, "floor_thickness") == 0) {
		config.floor_thickness = atof(value);
		CHECK_VALUE(config.floor_thickness >= 0.0, "floor thickness", ">= 0");
	}
	else if (strcmp(key, "min_shell_contact") == 0) {
		config.min_shell_contact = atof(value);
		CHECK_VALUE(config.min_shell_contact >= 0.0, "min_shell_contact", ">= 0");
	}
	else if (strcmp(key, "solid_fill_expansion") == 0) {
		config.solid_fill_expansion = atof(value);
		CHECK_VALUE(config.solid_fill_expansion >= 0.0, "solid fill expansion", ">= 0");
	}
	else if (strcmp(key, "material_diameter") == 0) {
		config.material_diameter = atof(value);
		CHECK_VALUE(config.material_diameter > 0.0, "material diameter", "> 0");
	}
	else if (strcmp(key, "flow_multiplier") == 0) {
		config.flow_multiplier = atof(value);
		CHECK_VALUE(config.flow_multiplier >= 0.0, "flow multiplier", ">= 0");
	}
	else if (strcmp(key, "feed_rate") == 0) {
		config.feed_rate = atof(value);
		CHECK_VALUE(config.feed_rate > 0.0, "feed rate", "> 0");
	}
	else if (strcmp(key, "perimeter_feed_rate") == 0) {
		config.perimeter_feed_rate = atof(value);
	}
	else if (strcmp(key, "loop_feed_rate") == 0) {
		config.loop_feed_rate = atof(value);
	}
	else if (strcmp(key, "infill_feed_rate") == 0) {
		config.infill_feed_rate = atof(value);
	}
	else if (strcmp(key, "support_feed_rate") == 0) {
		config.support_feed_rate = atof(value);
	}
	else if (strcmp(key, "travel_feed_rate") == 0) {
		config.travel_feed_rate = atof(value);
	}
	else if (strcmp(key, "first_layer_mult") == 0) {
		config.first_layer_mult = atof(value);
		CHECK_VALUE(config.first_layer_mult > 0.0, "first layer multiplier", "> 0");
	}
	else if (strcmp(key, "coast_len") == 0) {
		config.coast_len = atof(value);
		CHECK_VALUE(config.coast_len >= 0.0, "coast length", ">= 0");
	}
	else if (strcmp(key, "retract_len") == 0) {
		config.retract_len = atof(value);
		CHECK_VALUE(config.retract_len >= 0.0, "retract length", ">= 0");
	}
	else if (strcmp(key, "retract_speed") == 0) {
		config.retract_speed = atof(value);
		CHECK_VALUE(config.retract_speed > 0.0, "retract speed", "> 0");
	}
	else if (strcmp(key, "restart_speed") == 0) {
		config.restart_speed = atof(value);
	}
	else if (strcmp(key, "retract_min_travel") == 0) {
		config.retract_min_travel = atof(value);
		CHECK_VALUE(config.retract_min_travel >= 0.0, "minimum travel for retraction", ">= 0");
	}
	else if (strcmp(key, "retract_threshold") == 0) {
		config.retract_threshold = atof(value);
		CHECK_VALUE(config.retract_threshold >= 0.0, "unconditional retraction threshold", ">= 0");
	}
	else if (strcmp(key, "retract_within_island") == 0) {
		config.retract_within_island = PARSE_BOOL(value);
	}
	else if (strcmp(key, "extra_restart_len") == 0) {
		config.extra_restart_len = atof(value);
	}
	else if (strcmp(key, "cool_layer") == 0) {
		config.cool_layer = atoi(value);
	}
	else if (strcmp(key, "start_gcode") == 0) {
		free(config.start_gcode);
		config.start_gcode = strdup(value);
	}
	else if (strcmp(key, "end_gcode") == 0) {
		free(config.end_gcode);
		config.end_gcode = strdup(value);
	}
	else if (strcmp(key, "cool_on_gcode") == 0) {
		free(config.cool_on_gcode);
		config.cool_on_gcode = strdup(value);
	}
	else if (strcmp(key, "cool_off_gcode") == 0) {
		free(config.cool_off_gcode);
		config.cool_off_gcode = strdup(value);
	}
	else if (strcmp(key, "temp") == 0) {
		config.temp = atof(value);
	}
	else if (strcmp(key, "bed_temp") == 0) {
		config.bed_temp = atof(value);
	}
	else if (strcmp(key, "edge_overlap") == 0) {
		config.edge_overlap = atof(value);
		CHECK_VALUE(config.edge_overlap >= 0.0 && config.edge_overlap <= 1.0, "edge overlap", "within [0, 1]");
	}
	else if (strcmp(key, "strict_shell_order") == 0) {
		config.strict_shell_order = PARSE_BOOL(value);
	}
	else if (strcmp(key, "infill_first") == 0) {
		config.infill_first = PARSE_BOOL(value);
	}
	else if (strcmp(key, "align_seams") == 0) {
		config.align_seams = PARSE_BOOL(value);
	}
	else if (strcmp(key, "clean_insets") == 0) {
		config.clean_insets = PARSE_BOOL(value);
	}
	else if (strcmp(key, "fill_inset_gaps") == 0) {
		config.fill_inset_gaps = PARSE_BOOL(value);
	}
	else if (strcmp(key, "no_solid") == 0) {
		config.no_solid = PARSE_BOOL(value);
	}
	else if (strcmp(key, "outside_first") == 0) {
		config.outside_first = PARSE_BOOL(value);
	}
	else if (strcmp(key, "solid_infill_first") == 0) {
		config.solid_infill_first = PARSE_BOOL(value);
	}
	else if (strcmp(key, "separate_z_travel") == 0) {
		config.separate_z_travel = PARSE_BOOL(value);
	}
	else if (strcmp(key, "combine_all") == 0) {
		config.combine_all = PARSE_BOOL(value);
	}
	else if (strcmp(key, "generate_support") == 0) {
		config.generate_support = PARSE_BOOL(value);
	}
	else if (strcmp(key, "support_everywhere") == 0) {
		config.support_everywhere = PARSE_BOOL(value);
	}
	else if (strcmp(key, "solid_support_base") == 0) {
		config.solid_support_base = PARSE_BOOL(value);
	}
	else if (strcmp(key, "connect_support_lines") == 0) {
		config.connect_support_lines = PARSE_BOOL(value);
	}
	else if (strcmp(key, "poly_fill_type") == 0) {
		if (strcmp(value, "even_odd") == 0)
			config.poly_fill_type = ClipperLib::pftEvenOdd;
		else if (strcmp(value, "non_zero") == 0)
			config.poly_fill_type = ClipperLib::pftNonZero;
		else if (strcmp(value, "positive") == 0)
			config.poly_fill_type = ClipperLib::pftPositive;
		else if (strcmp(value, "negative") == 0)
			config.poly_fill_type = ClipperLib::pftNegative;
	}
	else if (strcmp(key, "inset_join_type") == 0) {
		if (parse_join_type(value, &config.inset_join_type))
			fprintf(stderr, "error: illegal value for inset_join_type: %s", value);
	}
	else if (strcmp(key, "outset_join_type") == 0) {
		if (parse_join_type(value, &config.outset_join_type))
			fprintf(stderr, "error: illegal value for outset_join_type: %s", value);
	}
	else if (strcmp(key, "offset_miter_limit") == 0) {
		config.offset_miter_limit = atof(value);
		CHECK_VALUE(config.offset_miter_limit >= 2.0, "offset miter limit", ">= 2");
	}
	else if (strcmp(key, "offset_arc_tolerance") == 0) {
		config.offset_arc_tolerance = atof(value);
		CHECK_VALUE(config.offset_arc_tolerance >= 0.25, "offset arc tolerance", ">= 0.25");
	}
	else if (strcmp(key, "fill_threshold") == 0) {
		config.fill_threshold = atof(value);
		CHECK_VALUE(config.fill_threshold >= 0.0, "fill threshold", ">= 0");
	}
	else if (strcmp(key, "support_angle") == 0) {
		config.support_angle = atof(value);
		CHECK_VALUE(config.support_angle > 0.0 && config.support_angle < 90.0, "support angle", "within (0, 90)");
	}
	else if (strcmp(key, "support_margin") == 0) {
		config.support_margin = atof(value);
		CHECK_VALUE(config.support_margin >= 0.0, "support margin", ">= 0");
	}
	else if (strcmp(key, "support_vert_margin") == 0) {
		config.support_vert_margin = atoi(value);
		CHECK_VALUE(config.support_vert_margin >= 0, "support vertical margin", ">= 0");
	}
	else if (strcmp(key, "interface_layers") == 0) {
		config.interface_layers = atoi(value);
		CHECK_VALUE(config.interface_layers >= 0, "interface layers", ">= 0");
	}
	else if (strcmp(key, "support_xy_expansion") == 0) {
		config.support_xy_expansion = atof(value);
		CHECK_VALUE(config.support_xy_expansion >= 0.0, "support xy expansion", ">= 0");
	}
	else if (strcmp(key, "support_density") == 0) {
		config.support_density = atof(value);
		CHECK_VALUE(config.support_density > 0.0 && config.support_density <= 1.0, "support density", "within (0, 1]");
	}
	else if (strcmp(key, "support_flow_mult") == 0) {
		config.support_flow_mult = atof(value);
		CHECK_VALUE(config.support_flow_mult > 0.0 && config.support_flow_mult <= 1.0, "support flow multiplier", "within (0, 1]");
	}
	else if (strcmp(key, "min_layer_time") == 0) {
		config.min_layer_time = atof(value);
		CHECK_VALUE(config.min_layer_time >= 0.0, "minimum layer time", ">= 0");
	}
	else if (strcmp(key, "layer_time_samples") == 0) {
		config.layer_time_samples = atoi(value);
		CHECK_VALUE(config.layer_time_samples >= 1, "number of layer time samples", ">= 1");
	}
	else if (strcmp(key, "min_feed_rate") == 0) {
		config.min_feed_rate = atof(value);
		CHECK_VALUE(config.min_feed_rate > 0.0, "minimum feed rate", "> 0");
	}
	else if (strcmp(key, "brim_width") == 0) {
		config.brim_width = atof(value);
		CHECK_VALUE(config.brim_width >= 0.0, "brim width", ">= 0");
	}
	else if (strcmp(key, "material_density") == 0) {
		config.material_density = atof(value);
		CHECK_VALUE(config.material_density >= 0.0, "material density", ">= 0");
	}
	else if (strcmp(key, "material_cost") == 0) {
		config.material_cost = atof(value);
		CHECK_VALUE(config.material_cost >= 0.0, "material cost", ">= 0");
	}
	else if (path) {
		fprintf(stderr, "error: line %d in %s: invalid option: %s\n", n, path, key);
		return 1;
	}
	else {
		fprintf(stderr, "error: invalid option: %s\n", key);
		return 1;
	}
	return 0;
}

static int read_config(const char *path)
{
	char *c, *key, *value, *next;
	c = get_file_contents(path);
	if (!c)
		return 1;
	key = c;
	for (ssize_t i = 1; *key != '\0'; ++i) {
		next = strchr(key, '\n');
		if (next) {
			while (next && (next[1] == ' ' || next[1] == '\t'))
				next = strchr(next + 1, '\n');
			next = isolate(next, '\n');
		}
		if (*key != '\0' && *key != '#') {
			value = isolate(key, '=');
			if (set_config_option(key, value, i, path))
				return 2;
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
	float stl_poly[13];
	FILE *f;

	if (strcmp(path, "-") == 0)
		f = stdin;
	else
		f = fopen(path, "rb");
	if (!f)
		return 1;
	fseek(f, 80, SEEK_SET);
	o->n = 0;
	fread(&o->n, 4, 1, f);
	o->t = (struct triangle *) calloc(o->n, sizeof(struct triangle));
	if (!o->t)
		die(e_nomem, 2);
	for (i = 0; i < o->n; ++i) {
		fread(stl_poly, 50, 1, f);  /* Each polygon entry is 50 bytes */
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
static void project2d(struct segment *seg, struct vertex *v0, struct vertex *v1, struct vertex *v2, fl_t z)
{
	seg->x[0] = v0->x + (v1->x - v0->x) * (z - v0->z) / (v1->z - v0->z);
	seg->y[0] = v0->y + (v1->y - v0->y) * (z - v0->z) / (v1->z - v0->z);
	seg->x[1] = v0->x + (v2->x - v0->x) * (z - v0->z) / (v2->z - v0->z);
	seg->y[1] = v0->y + (v2->y - v0->y) * (z - v0->z) / (v2->z - v0->z);
}

static void find_segments(struct slice *slices, struct triangle *t)
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

static void generate_islands(struct slice *slice, ClipperLib::PolyNode *n)
{
	for (ClipperLib::PolyNode *c : n->Childs) {
		struct island island = {};
		island.outlines.push_back(c->Contour);
		for (ClipperLib::PolyNode *cc : c->Childs) {
			island.outlines.push_back(cc->Contour);
			generate_islands(slice, cc);
		}
		slice->islands.push_back(island);
	}
}

#if USE_BOUNDING_BOX
static void find_bounding_box(struct island *island)
{
	bool first = true;
	for (ClipperLib::IntPoint &p : island->outlines[0]) {
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
#endif

static void generate_outlines(struct slice *slice, ssize_t slice_index)
{
	ssize_t i;
	struct list iseg = NEW_LIST, oseg = NEW_LIST;
	struct segment *s, *best, *begin, *end;
	bool flip_points;
	fl_t best_dist;
	ClipperLib::Paths outlines;

	for (i = 0; i < slice->n_seg; ++i)
		LIST_ADD_TAIL(&iseg, (struct node *) &slice->s[i]);

	while (iseg.h) {
		ssize_t segment_count = 0, flip_count = 0;
		/* ssize_t inexact_count = 0; */
		/* Add first segment to the polygon */
		s = (struct segment *) iseg.h;
		LIST_REMOVE_HEAD(&iseg);
		LIST_ADD_HEAD(&oseg, (struct node *) s);

		next_segment:
		++segment_count;
		best = NULL;
		flip_points = false;
		best_dist = HUGE_VAL;
		begin = (struct segment *) oseg.h;
		end = (struct segment *) oseg.t;

		/* Check whether the polygon is closed */
		if (begin != end) {
			if (begin->x[0] == end->x[1] && begin->y[0] == end->y[1])
				goto add_poly;
		}

		/* Link up connected segments */
		LIST_FOREACH(&iseg, s) {
			if (s->x[0] == end->x[1] && s->y[0] == end->y[1]) {
				LIST_REMOVE(&iseg, (struct node *) s);
				LIST_ADD_TAIL(&oseg, (struct node *) s);
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

				LIST_REMOVE(&iseg, (struct node *) s);
				LIST_ADD_TAIL(&oseg, (struct node *) s);
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

		/* Check whether the polygon is closed (within config.tolerance and less than best_dist) */
		if (begin != end) {
			fl_t close_dist = (begin->x[0] - end->x[1]) * (begin->x[0] - end->x[1]) + (begin->y[0] - end->y[1]) * (begin->y[0] - end->y[1]);
			if (close_dist <= config.tolerance && close_dist < best_dist)
				goto add_poly;
		}

		/* Connect nearest segment (within config.tolerance) */
		if (best && best_dist <= config.tolerance) {
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
			LIST_REMOVE(&iseg, (struct node *) s);
			LIST_ADD_TAIL(&oseg, (struct node *) s);
			/* ++inexact_count; */
			goto next_segment;
		}

		/* If there are any segments left and more than one output segment, there is probably a hole in the mesh */
		if (iseg.h && oseg.h != oseg.t)
			fprintf(stderr, "warning: there is (probably) a hole in the mesh at layer %zd (best_dist = %f)\n", slice_index + 1, sqrt(best_dist));
		goto next_poly;

		add_poly:
		if (oseg.h) {
			ClipperLib::Path poly;
			LIST_FOREACH(&oseg, s)
				poly.push_back(FL_T_TO_INTPOINT(s->x[0], s->y[0]));
			ClipperLib::CleanPolygon(poly, CLEAN_DIST);
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
		oseg.h = oseg.t = NULL;
	}
	free(slice->s);
	ClipperLib::PolyTree tree;
	ClipperLib::Clipper c;
	c.AddPaths(outlines, ClipperLib::ptSubject, true);
	c.Execute(ClipperLib::ctUnion, tree, config.poly_fill_type, config.poly_fill_type);
	generate_islands(slice, &tree);
#if USE_BOUNDING_BOX
	for (struct island &island : slice->islands)
		find_bounding_box(&island);
#endif
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

static void do_offset_square(ClipperLib::Paths &src, ClipperLib::Paths &dest, fl_t dist)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	co.AddPaths(src, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
	co.Execute(dest, FL_T_TO_CINT(dist));
}

static void generate_insets(struct slice *slice)
{
	for (struct island &island : slice->islands) {
		island.insets = (ClipperLib::Paths *) calloc(config.shells, sizeof(ClipperLib::Paths));
		if (!island.insets)
			die(e_nomem, 2);
		if (config.shells > 0) {
			fl_t offset = config.edge_offset + config.extra_offset;
			do_offset(island.outlines, island.insets[0], offset, 1.0 - config.edge_overlap);
			if (config.clean_insets)
				ClipperLib::CleanPolygons(island.insets[0], CLEAN_DIST);
			for (int i = 1; i < config.shells; ++i) {
				offset -= config.extrusion_width;
				do_offset(island.outlines, island.insets[i], offset, 1.0);
				if (config.clean_insets)
					ClipperLib::CleanPolygons(island.insets[i], CLEAN_DIST);
				if (island.insets[i].size() == 0)  /* break if nothing is being generated */
					goto done;
			}
			do_offset(island.outlines, island.infill_insets, offset - config.extrusion_width / 2.0, 0.0);
		}
		else {
			/* The offset distance here is not *technically* correct, but I'm not sure one can expect high dimensional accuracy when only printing infill anyway... */
			do_offset(island.outlines, island.infill_insets, config.edge_offset + config.extra_offset, 0.0);
		}
		ClipperLib::CleanPolygons(island.infill_insets, CLEAN_DIST);

		done:
		do_offset((config.shells > 0) ? island.insets[0] : island.infill_insets, island.boundaries, config.extrusion_width / 8.0, 0.0);
		if (config.solid_infill_clip_offset > 0.0)
			do_offset(island.infill_insets, island.solid_infill_clip, config.solid_infill_clip_offset, 0.0);
		else
			island.solid_infill_clip = island.infill_insets;
		if (config.shells > 1 && config.fill_inset_gaps) {
			ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
			ClipperLib::Paths hole;
			island.inset_gaps = (ClipperLib::Paths *) calloc(config.shells - 1, sizeof(ClipperLib::Paths));
			if (!island.inset_gaps)
				die(e_nomem, 2);
			for (int i = 0; i < config.shells - 1 && island.insets[i].size() > 0; ++i) {
				co.AddPaths(island.insets[i], config.inset_join_type, ClipperLib::etClosedPolygon);
				hole = island.insets[i + 1];
				ClipperLib::ReversePaths(hole);
				co.AddPaths(hole, config.inset_join_type, ClipperLib::etClosedPolygon);
				if (config.fill_threshold > 0.0) {
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT((config.extrusion_width + config.extrusion_width * config.fill_threshold) / -2.0));
					co.Clear();
					co.AddPaths(island.inset_gaps[i], config.outset_join_type, ClipperLib::etClosedPolygon);
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT(config.extrusion_width * config.fill_threshold / 2.0));
				}
				else {
					co.Execute(island.inset_gaps[i], FL_T_TO_CINT(-config.extrusion_width / 2.0));
				}
				co.Clear();
			}
		}
		if (config.align_seams) {
			for (int i = 0; i < config.shells; ++i) {
				for (ClipperLib::Path &p : island.insets[i]) {
					if (p.size() >= 3) {
						fl_t lowest = HUGE_VAL;
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

static void generate_even_line_fill(const ClipperLib::cInt min_x, const ClipperLib::cInt min_y, const ClipperLib::cInt max_x, const ClipperLib::cInt max_y, ClipperLib::Paths &p, fl_t density)
{
	ClipperLib::Path line(2);
	const ClipperLib::cInt move = FL_T_TO_CINT(config.extrusion_width * 1.4142135623 / density);
	line[0].X = min_x;
	line[0].Y = min_y + move;
	line[1].X = min_x + move;
	line[1].Y = min_y;
	while (line[0].X < max_x && line[1].Y < max_y) {
		p.push_back(line);
		if (line[0].Y >= max_y)
			line[0].X += move;
		else
			line[0].Y += move;

		if (line[1].X >= max_x)
			line[1].Y += move;
		else
			line[1].X += move;
	}
}

static void generate_odd_line_fill(const ClipperLib::cInt min_x, const ClipperLib::cInt min_y, const ClipperLib::cInt max_x, const ClipperLib::cInt max_y, ClipperLib::Paths &p, fl_t density)
{
	ClipperLib::Path line(2);
	const ClipperLib::cInt move = FL_T_TO_CINT(config.extrusion_width * 1.4142135623 / density);
	line[0].X = min_x;
	line[0].Y = max_y - move;
	line[1].X = min_x + move;
	line[1].Y = max_y;
	while (line[0].X < max_x && line[1].Y > min_y) {
		p.push_back(line);
		if (line[0].Y <= min_y)
			line[0].X += move;
		else
			line[0].Y -= move;

		if (line[1].X >= max_x)
			line[1].Y -= move;
		else
			line[1].X += move;
	}
}

static void generate_horizontal_line_fill(const ClipperLib::cInt min_x, const ClipperLib::cInt min_y, const ClipperLib::cInt max_x, const ClipperLib::cInt max_y, ClipperLib::Paths &p, fl_t density)
{
	ClipperLib::Path line(2);
	const ClipperLib::cInt move = FL_T_TO_CINT(config.extrusion_width / density);
	line[0].X = min_x;
	line[0].Y = min_y;
	line[1].X = max_x;
	line[1].Y = min_y;
	while (line[0].Y < max_y) {
		p.push_back(line);
		line[0].Y += move;
		line[1].Y += move;
	}
}

static void generate_vertical_line_fill(const ClipperLib::cInt min_x, const ClipperLib::cInt min_y, const ClipperLib::cInt max_x, const ClipperLib::cInt max_y, ClipperLib::Paths &p, fl_t density)
{
	ClipperLib::Path line(2);
	const ClipperLib::cInt move = FL_T_TO_CINT(config.extrusion_width / density);
	line[0].X = min_x;
	line[0].Y = min_y;
	line[1].X = min_x;
	line[1].Y = max_y;
	while (line[0].X < max_x) {
		p.push_back(line);
		line[0].X += move;
		line[1].X += move;
	}
}

static void generate_infill_patterns(struct object *o)
{
	const ClipperLib::cInt min_x = FL_T_TO_CINT(o->c.x - (o->w + config.xy_extra) / 2.0);
	const ClipperLib::cInt min_y = FL_T_TO_CINT(o->c.y - (o->d + config.xy_extra) / 2.0);
	const ClipperLib::cInt max_x = FL_T_TO_CINT(o->c.x + (o->w + config.xy_extra) / 2.0);
	const ClipperLib::cInt max_y = FL_T_TO_CINT(o->c.y + (o->d + config.xy_extra) / 2.0);

	generate_even_line_fill(min_x, min_y, max_x, max_y, o->solid_infill_patterns[0], 1.0);
	generate_odd_line_fill(min_x, min_y, max_x, max_y, o->solid_infill_patterns[1], 1.0);
	if (config.infill_density > 0.0) {
		generate_even_line_fill(min_x, min_y, max_x, max_y, o->sparse_infill_pattern, config.infill_density / 2.0);
		generate_odd_line_fill(min_x, min_y, max_x, max_y, o->sparse_infill_pattern, config.infill_density / 2.0);
	}
	if (config.generate_support) {
		generate_horizontal_line_fill(min_x, min_y, max_x, max_y, o->support_pattern, config.support_density);
		generate_vertical_line_fill(min_x, min_y, max_x, max_y, o->support_interface_pattern, 1.0);
	}
}

#if USE_BOUNDING_BOX
#define BOUNDING_BOX_INTERSECTS(a, b) (!((b).x0 > (a).x1 || (b).x1 < (a).x0 || (b).y0 < (a).y1 || (b).y1 > (a).y0))
#else
#define BOUNDING_BOX_INTERSECTS(a, b) true
#endif

static void generate_infill(struct object *o, ssize_t slice_index)
{
	for (struct island &island : o->slices[slice_index].islands) {
		ClipperLib::Clipper c;
		ClipperLib::PolyTree s;
		ClipperLib::Paths s_tmp;
		if (config.roof_layers > 0) {
			if (slice_index + 1 == o->n_slices)
				island.exposed_surface = island.infill_insets;
			else {
				c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
				for (struct island &clip_island : o->slices[slice_index + 1].islands)
					if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
						c.AddPaths((config.shells > 0) ? clip_island.insets[0] : clip_island.infill_insets, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctDifference, island.exposed_surface, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
			}
			if (island.exposed_surface.size() > 0)
				do_offset(island.exposed_surface, island.exposed_surface, config.extrusion_width / -2.0, 0.0);
		}
		if (config.infill_density == 1.0 || slice_index < config.floor_layers || slice_index + config.roof_layers >= o->n_slices) {
			if (config.fill_threshold > 0.0)
				remove_overlap(island.infill_insets, s_tmp, config.fill_threshold);
			c.AddPaths(o->solid_infill_patterns[slice_index % 2], ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			if (config.fill_inset_gaps)
				for (int i = 0; i < config.shells - 1; ++i)
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
		}
		else if (!config.no_solid && (config.floor_layers > 0 || config.roof_layers > 0)) {
			c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			for (int i = -config.floor_layers; i <= config.roof_layers; ++i) {
				if (i != 0) {
					for (struct island &clip_island : o->slices[slice_index + i].islands)
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
				do_offset_square(s_tmp, s_tmp, config.solid_infill_clip_offset + config.solid_fill_expansion * config.extrusion_width);
				c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
				c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
			}
			c.AddPaths(o->solid_infill_patterns[slice_index % 2], ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			if (config.fill_inset_gaps)
				for (int i = 0; i < config.shells - 1; ++i)
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			c.Clear();
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);

			c.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			c.Clear();
			if (config.fill_threshold > 0.0)
				remove_overlap(s_tmp, s_tmp, config.fill_threshold);
			c.AddPaths(o->sparse_infill_pattern, ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.sparse_infill);
		}
		else {
			if (config.fill_threshold > 0.0)
				remove_overlap(island.infill_insets, s_tmp, config.fill_threshold);
			c.AddPaths(o->sparse_infill_pattern, ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.sparse_infill);

			if (config.fill_inset_gaps) {
				c.Clear();
				c.AddPaths(o->solid_infill_patterns[slice_index % 2], ClipperLib::ptSubject, false);
				for (int i = 0; i < config.shells - 1; ++i)
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
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
	for (struct island &island : o->slices[slice_index - 1].islands)
		co.AddPaths((config.shells > 0) ? island.insets[0] : island.infill_insets, config.outset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(clip_paths, FL_T_TO_CINT(tan(config.support_angle / 180.0 * M_PI) * config.layer_height));
	co.Clear();
	for (struct island &island : o->slices[slice_index].islands)
		c.AddPaths((config.shells > 0) ? island.insets[0] : island.infill_insets, ClipperLib::ptSubject, true);
	c.AddPaths(clip_paths, ClipperLib::ptClip, true);
	c.Execute(ClipperLib::ctDifference, o->slices[slice_index].layer_support_map, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	c.Clear();
	co.AddPaths(o->slices[slice_index].layer_support_map, ClipperLib::jtSquare, ClipperLib::etClosedPolygon);
	co.Execute(o->slices[slice_index].layer_support_map, FL_T_TO_CINT(config.support_xy_expansion + (0.5 + config.support_margin) * config.edge_width - config.edge_offset));
}

static void generate_support_boundaries(struct slice *slice)
{
	ClipperLib::ClipperOffset co(config.offset_miter_limit, config.offset_arc_tolerance);
	for (struct island &island : slice->islands)
		co.AddPaths((config.shells > 0) ? island.insets[0] : island.infill_insets, config.outset_join_type, ClipperLib::etClosedPolygon);
	co.Execute(slice->support_boundaries, FL_T_TO_CINT((0.5 + config.support_margin) * config.edge_width - config.edge_offset));
}

static void generate_support_maps(struct object *o, ssize_t slice_index)
{
	for (ClipperLib::Path &p : o->slices[slice_index].layer_support_map) {
		for (ssize_t k = slice_index; k >= 0; --k) {
			ClipperLib::Clipper c;
			ClipperLib::Paths tmp;
			c.AddPath(p, ClipperLib::ptSubject, true);
			for (int i = (k >= config.support_vert_margin) ? -config.support_vert_margin : -k; k + i <= slice_index && i <= config.support_vert_margin; ++i)
				c.AddPaths(o->slices[k + i].support_boundaries, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctDifference, tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			if (tmp.size() >= 1) {
			#ifdef _OPENMP
				omp_set_lock(&o->slices[k].lock);
			#endif
				o->slices[k].support_map.insert(o->slices[k].support_map.end(), tmp.begin(), tmp.end());
			#ifdef _OPENMP
				omp_unset_lock(&o->slices[k].lock);
			#endif
			}
			else
				break;
		}
	}
}

static void union_support_maps(struct slice *slice)
{
	ClipperLib::SimplifyPolygons(slice->support_map, ClipperLib::pftNonZero);
}

static void remove_supports_not_touching_build_plate(struct object *o)
{
	ClipperLib::Clipper c;
	ClipperLib::Paths clip_paths;
	for (ssize_t i = 0; i < o->n_slices; ++i) {
		clip_paths.insert(clip_paths.end(), o->slices[i].support_boundaries.begin(), o->slices[i].support_boundaries.end());
		ClipperLib::SimplifyPolygons(clip_paths, ClipperLib::pftNonZero);
		c.AddPaths(o->slices[i].support_map, ClipperLib::ptSubject, true);
		c.AddPaths(clip_paths, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, o->slices[i].support_map, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
	}
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
	else if (config.interface_layers > 0) {
		ClipperLib::Paths s_tmp;
		c.AddPaths(slice->support_map, ClipperLib::ptSubject, true);
		for (int i = (slice_index > config.interface_layers) ? -config.interface_layers : -slice_index; slice_index + i < o->n_slices && i <= config.interface_layers; ++i) {
			if (i != 0) {
				c.AddPaths(o->slices[slice_index + i].support_map, ClipperLib::ptClip, true);
				c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
				c.Clear();
				if (i != config.interface_layers)
					c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
			}
		}
		c.AddPaths(slice->support_map, ClipperLib::ptSubject, true);
		c.AddPaths(s_tmp, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
		do_offset_square(s_tmp, s_tmp, config.extrusion_width / config.support_density);
		c.AddPaths(s_tmp, ClipperLib::ptSubject, true);
		c.AddPaths(slice->support_map, ClipperLib::ptClip, true);
		c.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
		c.Clear();
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
	for (int i = 1; i <= config.brim_lines; ++i) {
		ClipperLib::Paths tmp;
		for (struct island &island : o->slices[0].islands) {
			if (config.shells > 0)
				tmp.insert(tmp.end(), island.insets[0].begin(), island.insets[0].end());
			else
				tmp.insert(tmp.end(), island.infill_insets.begin(), island.infill_insets.end());
		}
		if (config.generate_support) {
			tmp.insert(tmp.end(), o->slices[0].support_map.begin(), o->slices[0].support_map.end());
			ClipperLib::SimplifyPolygons(tmp, ClipperLib::pftNonZero);
		}
		do_offset(tmp, tmp, config.extrusion_width * i, 1.0);
		o->brim.insert(o->brim.end(), tmp.begin(), tmp.end());
	}
}

static void slice_object(struct object *o)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> start;
	ssize_t i;
	o->n_slices = (ssize_t) ceil(o->h / config.layer_height);
	o->slices = (struct slice *) calloc(o->n_slices, sizeof(struct slice));
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
			generate_support_maps(o, i);
	#ifdef _OPENMP
		#pragma omp parallel for schedule(dynamic)
	#endif
		for (i = 0; i < o->n_slices; ++i)
			union_support_maps(&o->slices[i]);
		if (!config.support_everywhere)
			remove_supports_not_touching_build_plate(o);
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
#ifdef _OPENMP
	for (i = 0; i < o->n_slices; ++i)
		omp_destroy_lock(&o->slices[i].lock);
#endif

	fprintf(stderr, "sliced in %fs\n",
		(double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0);
}

static void preview_slices(struct object *o)
{
	ssize_t i;
	fprintf(stdout, "set size ratio -1\n");
	fprintf(stdout, "set xrange [%e:%e]\n", o->c.x - (o->w + config.xy_extra) / 2.0, o->c.x + (o->w + config.xy_extra) / 2.0);
	fprintf(stdout, "set yrange [%e:%e]\n", o->c.y - (o->d + config.xy_extra) / 2.0, o->c.y + (o->d + config.xy_extra) / 2.0);
	for (i = 0; i < o->n_slices; ++i) {
		fprintf(stderr, "layer %zd/%zd: intersections = %zd; islands = %zd\n", i + 1, o->n_slices, o->slices[i].n_seg, o->slices[i].islands.size());
		fprintf(stdout, "set title 'Layer %zd/%zd'\n", i + 1, o->n_slices);
		/* Draw outlines */
		fprintf(stdout, "plot \"-\" lc \"red\" dt 3 with lines title \"outline\", \"-\" lc \"blue\" with lines title \"insets\", \"-\" lc \"green\" with lines title \"infill\"\n");
		for (struct island &island : o->slices[i].islands) {
			for (ClipperLib::Path &path : island.outlines) {
				if (path.size() >= 3) {
					for (ClipperLib::IntPoint &p : path)
						fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
					fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / SCALE_CONSTANT, ((double) path[0].Y) / SCALE_CONSTANT);
					putc('\n', stdout);
				}
			}
		#if USE_BOUNDING_BOX
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x0) / SCALE_CONSTANT, ((double) island.box.y0) / SCALE_CONSTANT);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x1) / SCALE_CONSTANT, ((double) island.box.y0) / SCALE_CONSTANT);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x1) / SCALE_CONSTANT, ((double) island.box.y1) / SCALE_CONSTANT);
			fprintf(stdout, "%.4e %.4e\n", ((double) island.box.x0) / SCALE_CONSTANT, ((double) island.box.y1) / SCALE_CONSTANT);
			fprintf(stdout, "%.4e %.4e\n\n", ((double) island.box.x0) / SCALE_CONSTANT, ((double) island.box.y0) / SCALE_CONSTANT);
		#endif
		}
		/* Draw support map */
		for (ClipperLib::Path &path : o->slices[i].support_map) {
			if (path.size() >= 3) {
				for (ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
				fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / SCALE_CONSTANT, ((double) path[0].Y) / SCALE_CONSTANT);
				putc('\n', stdout);
			}
		}
		fprintf(stdout, "e\n");
		/* Draw insets */
		for (struct island &island : o->slices[i].islands) {
			for (int k = 0; k < config.shells; ++k) {
				for (ClipperLib::Path &path : island.insets[k]) {
					if (path.size() >= 3) {
						for (ClipperLib::IntPoint &p : path)
							fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
						fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / SCALE_CONSTANT, ((double) path[0].Y) / SCALE_CONSTANT);
						putc('\n', stdout);
					}
				}
			}
		}
		/* Draw brim */
		if (i == 0) {
			for (ClipperLib::Path &path : o->brim) {
				if (path.size() >= 3) {
					for (ClipperLib::IntPoint &p : path)
						fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
					fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / SCALE_CONSTANT, ((double) path[0].Y) / SCALE_CONSTANT);
					putc('\n', stdout);
				}
			}
		}
		fprintf(stdout, "e\n");
		/* Draw infill */
		for (struct island &island : o->slices[i].islands) {
			for (ClipperLib::Path &path : island.solid_infill) {
				for (ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
				putc('\n', stdout);
			}
			for (ClipperLib::Path &path : island.sparse_infill) {
				for (ClipperLib::IntPoint &p : path)
					fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
				putc('\n', stdout);
			}
		}
		/* Draw support lines */
		for (ClipperLib::Path &path : o->slices[i].support_lines) {
			for (ClipperLib::IntPoint &p : path)
				fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
			putc('\n', stdout);
		}
		for (ClipperLib::Path &path : o->slices[i].support_interface_lines) {
			for (ClipperLib::IntPoint &p : path)
				fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
			putc('\n', stdout);
		}
		fprintf(stdout, "e\n");
		fflush(stdout);
		fprintf(stderr, "press enter for next layer...");
		getchar();
	}
}

/* 0 is colinear, 1 is counter-clockwise and -1 is clockwise */
static int triplet_orientation(ClipperLib::IntPoint &a, ClipperLib::IntPoint &b, ClipperLib::IntPoint &c)
{
	ClipperLib::cInt v = a.X * b.Y - b.X * a.Y + b.X * c.Y - c.X * b.Y + c.X * a.Y - a.X * c.Y;  /* Calculate signed area * 2 */
	return (v == 0) ? 0 : (v > 0) ? 1 : -1;
}

static int is_on_segment(ClipperLib::IntPoint &a, ClipperLib::IntPoint &b, ClipperLib::IntPoint &c)
{
	if (b.X <= MAXIMUM(a.X, c.X) && b.X >= MINIMUM(a.X, c.X) && b.Y <= MAXIMUM(a.Y, c.Y) && b.Y >= MINIMUM(a.Y, c.Y))
		return true;
	return false;
}

static bool intersects(ClipperLib::IntPoint &a, ClipperLib::IntPoint &b, ClipperLib::IntPoint &c, ClipperLib::IntPoint &d)
{
	int o1 = triplet_orientation(a, b, c);
	int o2 = triplet_orientation(a, b, d);
	int o3 = triplet_orientation(c, d, a);
	int o4 = triplet_orientation(c, d, b);
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

static ssize_t get_boundary_crossing(ClipperLib::Path &p, ClipperLib::IntPoint &p0, ClipperLib::IntPoint &p1)
{
	for (size_t i = 1; i < p.size(); ++i) {
		if (intersects(p[i - 1], p[i], p0, p1))
			return (ssize_t) i - 1;
	}
	if (intersects(p[p.size() - 1], p[0], p0, p1))
		return p.size() - 1;
	return -1;
}

static bool crosses_boundary(struct machine *m, struct island *island, ClipperLib::cInt x, ClipperLib::cInt y)
{
	ClipperLib::IntPoint p0(m->x, m->y);
	ClipperLib::IntPoint p1(x, y);
	for (ClipperLib::Path &p : island->boundaries) {
		if (get_boundary_crossing(p, p0, p1) >= 0)
			return true;
	}
	return false;
}

static bool crosses_exposed_surface(struct machine *m, struct island *island, ClipperLib::cInt x, ClipperLib::cInt y)
{
	ClipperLib::IntPoint p0(m->x, m->y);
	ClipperLib::IntPoint p1(x, y);
	for (ClipperLib::Path &p : island->exposed_surface) {
		if (get_boundary_crossing(p, p0, p1) >= 0)
			return true;
		else if (ClipperLib::PointInPolygon(p0, p) || ClipperLib::PointInPolygon(p1, p))
			return true;
	}
	return false;
}

static size_t find_nearest_path(ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist, size_t *r_start)
{
	size_t best = 0, start = 0;
	fl_t best_dist = HUGE_VAL;
	fl_t x0 = CINT_TO_FL_T(x), y0 = CINT_TO_FL_T(y);
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
		*r_dist = best_dist;
	if (r_start)
		*r_start = start;
	return best;
}

static size_t find_nearest_aligned_path(ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist)
{
	size_t best = 0;
	fl_t best_dist = HUGE_VAL;
	fl_t x0 = CINT_TO_FL_T(x), y0 = CINT_TO_FL_T(y);
	for (size_t i = 0; i < p.size(); ++i) {
		fl_t x1 = CINT_TO_FL_T(p[i][0].X), y1 = CINT_TO_FL_T(p[i][0].Y);
		fl_t dist = (x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0);
		if (dist < best_dist) {
			best_dist = dist;
			best = i;
		}
	}
	if (r_dist)
		*r_dist = best_dist;
	return best;
}

static size_t find_nearest_segment(ClipperLib::Paths &p, ClipperLib::cInt x, ClipperLib::cInt y, fl_t *r_dist, bool *r_flip)
{
	bool flip = false;
	size_t best = 0;
	fl_t best_dist = HUGE_VAL;
	fl_t x0 = CINT_TO_FL_T(x), y0 = CINT_TO_FL_T(y);
	for (size_t i = 0; i < p.size(); ++i) {
		fl_t x1_0 = CINT_TO_FL_T(p[i][0].X), y1_0 = CINT_TO_FL_T(p[i][0].Y);
		fl_t x1_1 = CINT_TO_FL_T(p[i][1].X), y1_1 = CINT_TO_FL_T(p[i][1].Y);
		fl_t dist0 = (x1_0 - x0) * (x1_0 - x0) + (y1_0 - y0) * (y1_0 - y0);
		fl_t dist1 = (x1_1 - x0) * (x1_1 - x0) + (y1_1 - y0) * (y1_1 - y0);
		if (dist0 < best_dist) {
			best_dist = dist0;
			best = i;
			flip = false;
		}
		if (dist1 < best_dist) {
			best_dist = dist1;
			best = i;
			flip = true;
		}
	}
	if (r_dist)
		*r_dist = best_dist;
	if (r_flip)
		*r_flip = flip;
	return best;
}

static void append_g_move(struct slice *slice, struct g_move &move, fl_t len)
{
	/* FIXME: should probably take acceleration into account... */
	slice->layer_time += len / move.feed_rate;
	slice->moves.push_back(move);
}

static void linear_move(struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt x, ClipperLib::cInt y, ClipperLib::cInt z, fl_t extra_e_len, fl_t feed_rate, fl_t flow_adjust, bool scalable, bool is_travel, bool doing_infill)
{
	fl_t f_x = CINT_TO_FL_T(x), f_y = CINT_TO_FL_T(y), f_z = CINT_TO_FL_T(z);
	fl_t f_mx = CINT_TO_FL_T(m->x), f_my = CINT_TO_FL_T(m->y), f_mz = CINT_TO_FL_T(m->z);
	struct g_move move = { x, y, z, 0.0, feed_rate, scalable, is_travel, false };
	fl_t len = sqrt((f_mx - f_x) * (f_mx - f_x) + (f_my - f_y) * (f_my - f_y) + (f_mz - f_z) * (f_mz - f_z));
	if (is_travel) {
		if (!m->is_retracted && config.retract_len > 0.0
			&& (m->force_retract
				|| len > ((doing_infill) ? config.retract_threshold : config.retract_min_travel)
				|| (config.retract_within_island && len > config.retract_min_travel)
				|| (island && crosses_boundary(m, island, x, y))
				|| (island && len > config.extrusion_width * 2.0 && crosses_exposed_surface(m, island, x, y))
			)
		) {
			struct g_move retract_move = { m->x, m->y, m->z, -config.retract_len, config.retract_speed, false, false, false };
			append_g_move(slice, retract_move, config.retract_len);
			m->is_retracted = true;
		}
		m->force_retract = false;
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
}

#if 0
/* This function is unused at the moment, but may be useful in the future */
static fl_t get_partial_path_len(ClipperLib::Path &p, size_t start, size_t end, bool reverse)
{
	fl_t l = 0.0;
	fl_t x0 = CINT_TO_FL_T(p[start].X), y0 = CINT_TO_FL_T(p[start].Y);
	size_t i = start;
	do {
		if (reverse)
			i = (i < 1) ? p.size() - 1 : i - 1;
		else
			i = (i >= p.size() - 1) ? 0 : i + 1;
		fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		x0 = x1;
		y0 = y1;
	} while (i != end);
	return l;
}
#endif

static bool path_len_is_greater_than(ClipperLib::Path &p, fl_t len)
{
	fl_t l = 0.0;
	fl_t x0 = CINT_TO_FL_T(p[0].X), y0 = CINT_TO_FL_T(p[0].Y);
	for (size_t i = 1; i < p.size(); ++i) {
		fl_t x1 = CINT_TO_FL_T(p[i].X), y1 = CINT_TO_FL_T(p[i].Y);
		l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		if (l > len)
			return true;
		x0 = x1;
		y0 = y1;
	}
	fl_t x1 = CINT_TO_FL_T(p[0].X), y1 = CINT_TO_FL_T(p[0].Y);
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
		fl_t x1 = CINT_TO_FL_T((p.end() - 1)->X), y1 = CINT_TO_FL_T((p.end() - 1)->Y);
		fl_t xv = x1 - x0, yv = y1 - y0;
		fl_t norm = sqrt(xv * xv + yv * yv);
		l += norm;
		if (l == clip)
			break;
		else if (l > clip) {
			fl_t new_x = x1 - (l - clip) * (xv / norm), new_y = y1 - (l - clip) * (yv / norm);
			p.push_back(FL_T_TO_INTPOINT(new_x, new_y));
			break;
		}
		x0 = x1;
		y0 = y1;
	}
	if (clipped_points)
		std::reverse(clipped_points->begin(), clipped_points->end());
}

static void generate_closed_path_moves(ClipperLib::Path &p, size_t start_idx, struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, fl_t feed_rate)
{
	if (p.size() < 3)
		return;
	bool first_point = true;
	bool do_clip = (config.shell_clip > 0.0 && path_len_is_greater_than(p, config.shell_clip * 4.0));
	bool do_coast = (config.coast_len > 0.0 && path_len_is_greater_than(p, config.coast_len * 2.0 + ((do_clip) ? config.shell_clip : 0.0)));
	if (start_idx != 0)
		std::rotate(p.begin(), p.begin() + start_idx, p.end());
	p.push_back(p[0]);
	if (do_clip)
		clip_path_from_end(p, NULL, config.shell_clip);
	ClipperLib::Path coast_path;
	if (do_coast)
		clip_path_from_end(p, &coast_path, config.coast_len);
	for (ClipperLib::IntPoint &point : p) {
		if (first_point) {
			linear_move(slice, island, m, point.X, point.Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, false);
			first_point = false;
		}
		else
			linear_move(slice, island, m, point.X, point.Y, z, 0.0, feed_rate, 1.0, true, false, false);
	}
	m->is_retracted = true;  /* Make sure we don't retract */
	for (ClipperLib::IntPoint &point : coast_path)
		linear_move(slice, island, m, point.X, point.Y, z, 0.0, feed_rate, 1.0, true, true ,false);
	m->is_retracted = false;
}

static void plan_brim(struct object *o, struct machine *m, ClipperLib::cInt z)
{
	while (!o->brim.empty()) {
		size_t best = 0, start = 0;
		find_nearest_path(o->brim, m->x, m->y, NULL, &start);
		generate_closed_path_moves(o->brim[best], start, &o->slices[0], NULL, m, z, config.perimeter_feed_rate);
		o->brim.erase(o->brim.begin() + best);
	}
}

static void plan_support(struct slice *slice, ClipperLib::Paths &lines, struct machine *m, ClipperLib::cInt z, ssize_t layer_num, fl_t min_len, fl_t connect_threshold)
{
	ClipperLib::Path last_line(2);
	bool first = true;
	fl_t flow_adjust = (layer_num > 0) ? config.support_flow_mult : 1.0;
	fl_t feed_rate = (layer_num > 0) ? config.support_feed_rate : config.perimeter_feed_rate;
	while (!lines.empty()) {
		bool flip_points;
		fl_t best_dist;
		size_t best = find_nearest_segment(lines, m->x, m->y, &best_dist, &flip_points);
		ClipperLib::Path &p = lines[best];
		fl_t x0 = CINT_TO_FL_T(p[0].X), y0 = CINT_TO_FL_T(p[0].Y);
		fl_t x1 = CINT_TO_FL_T(p[1].X), y1 = CINT_TO_FL_T(p[1].Y);
		fl_t len = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		if (len > min_len) {
			bool crosses_boundary = false;
			if (!first) {
				ClipperLib::IntPoint p0(m->x, m->y);
				for (ClipperLib::Path &bound : slice->support_boundaries) {
					if (get_boundary_crossing(bound, last_line[0], last_line[1]) >= 0 || get_boundary_crossing(bound, p0, (flip_points) ? p[1] : p[0]) >= 0) {
						crosses_boundary = true;
						m->force_retract = true;
						break;
					}
				}
			}
			bool connect = (!first && !crosses_boundary && best_dist < connect_threshold);
			if (flip_points) {
				if (connect)
					linear_move(slice, NULL, m, p[1].X, p[1].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
				else
					linear_move(slice, NULL, m, p[1].X, p[1].Y, z, 0.0, config.travel_feed_rate, flow_adjust, false, true, true);
				linear_move(slice, NULL, m, p[0].X, p[0].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
			}
			else {
				if (connect)
					linear_move(slice, NULL, m, p[0].X, p[0].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
				else
					linear_move(slice, NULL, m, p[0].X, p[0].Y, z, 0.0, config.travel_feed_rate, flow_adjust, false, true, true);
				linear_move(slice, NULL, m, p[1].X, p[1].Y, z, 0.0, feed_rate, flow_adjust, true, false, true);
			}
			last_line = p;
			first = false;
		}
		lines.erase(lines.begin() + best);
	}
}

static void plan_insets_weighted(struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z, bool outside_first)
{
	for (;;) {
		bool done = true;
		fl_t best_dist = HUGE_VAL;
		size_t best = 0, inset = 0, start = 0;
		for (int i = 0; i < config.shells; ++i) {
			if (!island->insets[i].empty()) {
				fl_t dist;
				size_t r, start_tmp = 0;
				if (config.align_seams)
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
			i = (outside_first) ? i + i : i - 1;
			continue;
		}
		size_t best, start = 0;
		if (config.align_seams)
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
}

static void plan_infill(ClipperLib::Paths &lines, struct slice *slice, struct island *island, struct machine *m, ClipperLib::cInt z)
{
	while (!lines.empty()) {
		bool flip_points;
		size_t best = find_nearest_segment(lines, m->x, m->y, NULL, &flip_points);
		ClipperLib::Path &p = lines[best];
		if (flip_points) {
			linear_move(slice, island, m, p[1].X, p[1].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
			linear_move(slice, island, m, p[0].X, p[0].Y, z, 0.0, config.infill_feed_rate, 1.0, true, false, true);
		}
		else {
			linear_move(slice, island, m, p[0].X, p[0].Y, z, 0.0, config.travel_feed_rate, 1.0, false, true, true);
			linear_move(slice, island, m, p[1].X, p[1].Y, z, 0.0, config.infill_feed_rate, 1.0, true, false, true);
		}
		lines.erase(lines.begin() + best);
	}
}

static void plan_moves(struct object *o, struct slice *slice, ssize_t layer_num, struct machine *m)
{
	ClipperLib::cInt z = FL_T_TO_CINT(((fl_t) layer_num) * config.layer_height + config.layer_height);
	if (layer_num == 0 && config.brim_lines > 0)
		plan_brim(o, m, z);
	if (config.generate_support) {
		plan_support(slice, slice->support_interface_lines, m, z, layer_num, config.extrusion_width, (config.connect_support_lines) ? config.extrusion_width * 1.9 : 0.0);
		plan_support(slice, slice->support_lines, m, z, layer_num, config.extrusion_width * 2.0, (config.connect_support_lines) ? (layer_num == 0 && config.solid_support_base) ? config.extrusion_width * 1.9 : config.extrusion_width / config.support_density * 10.0 : 0.0);
	}
	while (slice->islands.size() > 0) {
		size_t best = 0;
		fl_t best_dist = HUGE_VAL;
		for (size_t i = 0; i < slice->islands.size(); ++i) {
			fl_t dist;
			if (config.align_seams && config.shells > 0)
				find_nearest_aligned_path(slice->islands[i].insets[0], m->x, m->y, &dist);
			else
				find_nearest_path((config.shells > 0) ? slice->islands[i].insets[0] : slice->islands[i].infill_insets, m->x, m->y, &dist, NULL);
			if (dist < best_dist) {
				best = i;
				best_dist = dist;
			}
		}
		struct island &island = slice->islands[best];
		if (!config.solid_infill_first) {
			/* Concat sparse infill into solid infill vector so the infill pathing is planned together */
			island.solid_infill.insert(island.solid_infill.end(), island.sparse_infill.begin(), island.sparse_infill.end());
			island.sparse_infill.clear();
		}
		if (config.infill_first && layer_num != 0) {
			plan_infill(island.solid_infill, slice, &island, m, z);
			if (config.solid_infill_first)
				plan_infill(island.sparse_infill, slice, &island, m, z);
			plan_insets(slice, &island, m, z, config.outside_first || layer_num == 0);
		}
		else {
			plan_insets(slice, &island, m, z, config.outside_first || layer_num == 0);
			plan_infill(island.solid_infill, slice, &island, m, z);
			if (config.solid_infill_first)
				plan_infill(island.sparse_infill, slice, &island, m, z);
		}
		free(island.insets);
		free(island.inset_gaps);
		slice->islands.erase(slice->islands.begin() + best);
	}
	m->force_retract = true;  /* force retract on layer change */
}

static void write_gcode_string(const char *s, FILE *f)
{
	bool line_start = true;
	if (!s || strlen(s) == 0)
		return;
	while (*s) {
		if (*s == '%') {
			switch (*(++s)) {
			case 't':
				fprintf(f, "%.1f", config.temp);
				break;
			case 'b':
				fprintf(f, "%.1f", config.bed_temp);
				break;
			case 'R':
				fprintf(f, "%.5f", config.retract_len);
				break;
			case '%':
				fputc('%', f);
				break;
			}
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
	fputc('\n', f);
}

static void write_gcode_move(FILE *f, struct g_move *move, struct machine *m, fl_t feed_rate_mult, bool force_xyz)
{
	fl_t feed_rate = move->feed_rate;
	if (move->scalable) {
		feed_rate *= feed_rate_mult;
		fl_t min_feed_rate = (move->is_restart) ? config.min_feed_rate * config.extrusion_area / config.material_area : config.min_feed_rate;
		if (feed_rate < min_feed_rate)
			feed_rate = min_feed_rate;
	}
	if (move->is_travel && move->z != m->z && config.separate_z_travel) {
		fprintf(f, "G1 Z%.3f", CINT_TO_FL_T(move->z));
		if (feed_rate != m->feed_rate)
			fprintf(f, " F%ld", lround(feed_rate * 60.0));
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
	if (feed_rate != m->feed_rate) {
		long int feed_rate_min = lround(feed_rate * 60.0);
		fprintf(f, " F%ld", (feed_rate_min < 1) ? 1 : feed_rate_min);
	}
	fputc('\n', f);
	m->x = move->x;
	m->y = move->y;
	m->z = move->z;
	m->e += move->e;
	m->feed_rate = feed_rate;
}

static int write_gcode(const char *path, struct object *o)
{
	FILE *f;
	if (strcmp(path, "-") == 0)
		f = stdout;
	else
		f = fopen(path, "w");
	if (!f)
		return 1;
	struct machine plan_m = {};
	plan_m.x = FL_T_TO_CINT(o->c.x - (o->w + config.xy_extra) / 2.0);
	plan_m.y = FL_T_TO_CINT(o->c.y - (o->d + config.xy_extra) / 2.0);
	plan_m.is_retracted = true;
	bool is_first_move = true;
	struct machine export_m = {};
	fl_t total_e = 0.0;
	fl_t feed_rate_mult = config.first_layer_mult;
	fprintf(stderr, "plan moves and write gcode to %s...", path);
	write_gcode_string(config.start_gcode, f);
	for (ssize_t i = 0; i < o->n_slices; ++i) {
		struct slice *slice = &o->slices[i];
		plan_moves(o, slice, i, &plan_m);
		fprintf(f, "; layer %zd (z = %f)\n", i, ((fl_t) i) * config.layer_height + config.layer_height);
		if (i == config.cool_layer)
			write_gcode_string(config.cool_on_gcode, f);
		fl_t average_layer_time = slice->layer_time;
		for (int k = 1; k < config.layer_time_samples; ++k)
			average_layer_time += (k < i) ? o->slices[i - k].layer_time : o->slices[0].layer_time;
		average_layer_time /= config.layer_time_samples;
		if (average_layer_time < config.min_layer_time)
			feed_rate_mult = average_layer_time / config.min_layer_time;
		for (struct g_move &move : o->slices[i].moves) {
			write_gcode_move(f, &move, &export_m, feed_rate_mult, is_first_move);
			is_first_move = false;
		}
		slice->moves.clear();
		feed_rate_mult = 1.0;
		total_e += export_m.e;
		export_m.e = 0.0;
		fputs("G92 E0\n", f);
	}
	write_gcode_string(config.cool_off_gcode, f);
	write_gcode_string(config.end_gcode, f);
	fputs(" done\n", stderr);
	fl_t mass = config.material_area * total_e * config.material_density / config.flow_multiplier;
	fprintf(f, "; material length = %.4f\n", total_e / config.flow_multiplier);
	fprintf(f, "; material mass   = %.4f\n", mass);
	fprintf(f, "; material cost   = %.4f\n", mass * config.material_cost);
	fprintf(stderr, "material length = %.4f\n", total_e / config.flow_multiplier);
	fprintf(stderr, "material mass   = %.4f\n", mass);
	fprintf(stderr, "material cost   = %.4f\n", mass * config.material_cost);
	long int bytes = ftell(f);
	if (bytes >= 2048 * 1024)
		fprintf(stderr, "wrote %.2fMiB\n", bytes / 1024.0 / 1024.0);
	else if (bytes >= 2048)
		fprintf(stderr, "wrote %.2fKiB\n", bytes / 1024.0);
	else
		fprintf(stderr, "wrote %ldB\n", bytes);
	fclose(f);
	return 0;
}

static const char * get_poly_fill_type_string(ClipperLib::PolyFillType pft)
{
	switch (pft) {
	case ClipperLib::pftEvenOdd:  return "even_odd";
	case ClipperLib::pftNonZero:  return "non_zero";
	case ClipperLib::pftPositive: return "positive";
	case ClipperLib::pftNegative: return "negative";
	}
	return NULL;
}

static const char * get_join_type_string(ClipperLib::JoinType jt)
{
	switch (jt) {
	case ClipperLib::jtMiter:  return "miter";
	case ClipperLib::jtSquare: return "square";
	case ClipperLib::jtRound:  return "round";
	}
	return NULL;
}

#define GET_FEED_RATE(x, m) (((x) >= 0.0) ? (x) : (m) * -(x))

int main(int argc, char *argv[])
{
	int opt;
	char *path, *output_path = NULL;
	struct object o;
	fl_t scale_factor = 1.0, x_translate = 0.0, y_translate = 0.0, z_chop = 0.0;
	bool do_preview = false;

	/* Parse options */
	while ((opt = getopt(argc, argv, ":hpo:c:O:l:w:t:s:d:n:r:f:b:C:x:y:z:")) != -1) {
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
			output_path = strdup(optarg);
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
			key = strdup(optarg);
			value = isolate(key, '=');
			if (set_config_option(key, value, 0, NULL))
				return 1;
			free(key);
			break;
		case 'l':
			if (set_config_option("layer_height", optarg, 0, NULL))
				return 1;
			break;
		case 'w':
			if (set_config_option("extrusion_width", optarg, 0, NULL))
				return 1;
			break;
		case 't':
			if (set_config_option("tolerance", optarg, 0, NULL))
				return 1;
			break;
		case 's':
			scale_factor = atof(optarg);
			CHECK_VALUE(scale_factor != 0.0, "scale factor", "!= 0");
			break;
		case 'd':
			if (set_config_option("infill_density", optarg, 0, NULL))
				return 1;
			break;
		case 'n':
			if (set_config_option("shells", optarg, 0, NULL))
				return 1;
			break;
		case 'r':
			if (set_config_option("roof_thickness", optarg, 0, NULL))
				return 1;
			break;
		case 'f':
			if (set_config_option("floor_thickness", optarg, 0, NULL))
				return 1;
			break;
		case 'b':
			if (set_config_option("brim_width", optarg, 0, NULL))
				return 1;
			break;
		case 'C':
			if (set_config_option("coarseness", optarg, 0, NULL))
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

	/* config.tolerance needs to be squared */
	config.tolerance *= config.tolerance;
	config.roof_layers = lround(config.roof_thickness / config.layer_height);
	config.floor_layers = lround(config.floor_thickness / config.layer_height);
	config.extrusion_area = config.extrusion_width * config.layer_height * config.packing_density;
	config.edge_width = (config.extrusion_area - (config.layer_height * config.layer_height * M_PI / 4.0)) / config.layer_height + config.layer_height;
	config.edge_offset = config.edge_width / -2.0 - (config.extrusion_area * (1.0 - config.edge_packing_density)) / config.layer_height;
	config.shell_clip = (config.extrusion_width * config.packing_density) * M_PI / 4.0 * (1.0 - config.seam_packing_density);
	config.material_area = M_PI * config.material_diameter * config.material_diameter / 4.0;
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
	/* set feed rates */
	config.perimeter_feed_rate = GET_FEED_RATE(config.perimeter_feed_rate, config.feed_rate);
	config.loop_feed_rate = GET_FEED_RATE(config.loop_feed_rate, config.feed_rate);
	config.infill_feed_rate = GET_FEED_RATE(config.infill_feed_rate, config.feed_rate);
	config.support_feed_rate = GET_FEED_RATE(config.support_feed_rate, config.feed_rate);
	config.travel_feed_rate = GET_FEED_RATE(config.travel_feed_rate, config.feed_rate);
	config.restart_speed = GET_FEED_RATE(config.restart_speed, config.retract_speed);

	fprintf(stderr, "configuration:\n");
	fprintf(stderr, "  layer_height          = %f\n", config.layer_height);
	fprintf(stderr, "  tolerance             = %f\n", sqrt(config.tolerance));
	fprintf(stderr, "  coarseness            = %f (%f units)\n", config.coarseness, config.coarseness / SCALE_CONSTANT);
	fprintf(stderr, "  extrusion_width       = %f\n", config.extrusion_width);
	fprintf(stderr, " *edge_width            = %f\n", config.edge_width);
	fprintf(stderr, " *extrusion_area        = %f\n", config.extrusion_area);
	fprintf(stderr, "  scale_factor (-s)     = %f\n", scale_factor);
	fprintf(stderr, "  xy_scale_factor       = %f\n", config.xy_scale_factor);
	fprintf(stderr, "  z_scale_factor        = %f\n", config.z_scale_factor);
	fprintf(stderr, "  x_center              = %f\n", config.x_center);
	fprintf(stderr, "  y_center              = %f\n", config.y_center);
	fprintf(stderr, "  packing_density       = %f\n", config.packing_density);
	fprintf(stderr, "  edge_packing_density  = %f\n", config.edge_packing_density);
	fprintf(stderr, "  seam_packing_density  = %f\n", config.seam_packing_density);
	fprintf(stderr, "  extra_offset          = %f\n", config.extra_offset);
	fprintf(stderr, " *edge_offset           = %f\n", config.edge_offset);
	fprintf(stderr, " *shell_clip            = %f\n", config.shell_clip);
	fprintf(stderr, "  infill_density        = %f\n", config.infill_density);
	fprintf(stderr, "  shells                = %d\n", config.shells);
	fprintf(stderr, "  roof_thickness        = %f\n", config.roof_thickness);
	fprintf(stderr, " *roof_layers           = %d\n", config.roof_layers);
	fprintf(stderr, "  floor_thickness       = %f\n", config.floor_thickness);
	fprintf(stderr, " *floor_layers          = %d\n", config.floor_layers);
	fprintf(stderr, "  min_shell_contact     = %f\n", config.min_shell_contact);
	fprintf(stderr, "  material_diameter     = %f\n", config.material_diameter);
	fprintf(stderr, " *material_area         = %f\n", config.material_area);
	fprintf(stderr, "  flow_multiplier       = %f\n", config.flow_multiplier);
	fprintf(stderr, "  feed_rate             = %f\n", config.feed_rate);
	fprintf(stderr, "  perimeter_feed_rate   = %f\n", config.perimeter_feed_rate);
	fprintf(stderr, "  loop_feed_rate        = %f\n", config.loop_feed_rate);
	fprintf(stderr, "  infill_feed_rate      = %f\n", config.infill_feed_rate);
	fprintf(stderr, "  support_feed_rate     = %f\n", config.support_feed_rate);
	fprintf(stderr, "  travel_feed_rate      = %f\n", config.travel_feed_rate);
	fprintf(stderr, "  first_layer_mult      = %f\n", config.first_layer_mult);
	fprintf(stderr, "  coast_len             = %f\n", config.coast_len);
	fprintf(stderr, "  retract_len           = %f\n", config.retract_len);
	fprintf(stderr, "  retract_speed         = %f\n", config.retract_speed);
	fprintf(stderr, "  restart_speed         = %f\n", config.restart_speed);
	fprintf(stderr, "  retract_min_travel    = %f\n", config.retract_min_travel);
	fprintf(stderr, "  retract_threshold     = %f\n", config.retract_threshold);
	fprintf(stderr, "  retract_within_island = %s\n", (config.retract_within_island) ? "true" : "false");
	fprintf(stderr, "  extra_restart_len     = %f\n", config.extra_restart_len);
	fprintf(stderr, "  cool_layer            = %d\n", config.cool_layer);
	fprintf(stderr, "  temp                  = %f\n", config.temp);
	fprintf(stderr, "  bed_temp              = %f\n", config.bed_temp);
	fprintf(stderr, "  edge_overlap          = %f\n", config.edge_overlap);
	fprintf(stderr, "  strict_shell_order    = %s\n", (config.strict_shell_order) ? "true" : "false");
	fprintf(stderr, "  infill_first          = %s\n", (config.infill_first) ? "true" : "false");
	fprintf(stderr, "  align_seams           = %s\n", (config.align_seams) ? "true" : "false");
	fprintf(stderr, "  clean_insets          = %s\n", (config.clean_insets) ? "true" : "false");
	fprintf(stderr, "  fill_inset_gaps       = %s\n", (config.fill_inset_gaps) ? "true" : "false");
	fprintf(stderr, "  no_solid              = %s\n", (config.no_solid) ? "true" : "false");
	fprintf(stderr, "  outside_first         = %s\n", (config.outside_first) ? "true" : "false");
	fprintf(stderr, "  solid_infill_first    = %s\n", (config.solid_infill_first) ? "true" : "false");
	fprintf(stderr, "  separate_z_travel     = %s\n", (config.separate_z_travel) ? "true" : "false");
	fprintf(stderr, "  combine_all           = %s\n", (config.combine_all) ? "true" : "false");
	fprintf(stderr, "  generate_support      = %s\n", (config.generate_support) ? "true" : "false");
	fprintf(stderr, "  support_everywhere    = %s\n", (config.support_everywhere) ? "true" : "false");
	fprintf(stderr, "  solid_support_base    = %s\n", (config.solid_support_base) ? "true" : "false");
	fprintf(stderr, "  connect_support_lines = %s\n", (config.connect_support_lines) ? "true" : "false");
	fprintf(stderr, "  poly_fill_type        = %s\n", get_poly_fill_type_string(config.poly_fill_type));
	fprintf(stderr, "  inset_join_type       = %s\n", get_join_type_string(config.inset_join_type));
	fprintf(stderr, "  outset_join_type      = %s\n", get_join_type_string(config.outset_join_type));
	fprintf(stderr, "  offset_miter_limit    = %f\n", config.offset_miter_limit);
	fprintf(stderr, "  offset_arc_tolerance  = %f\n", config.offset_arc_tolerance);
	fprintf(stderr, "  fill_threshold        = %f\n", config.fill_threshold);
	fprintf(stderr, "  support_angle         = %f\n", config.support_angle);
	fprintf(stderr, "  support_margin        = %f\n", config.support_margin);
	fprintf(stderr, "  support_vert_margin   = %d\n", config.support_vert_margin);
	fprintf(stderr, "  interface_layers      = %d\n", config.interface_layers);
	fprintf(stderr, "  support_xy_expansion  = %f\n", config.support_xy_expansion);
	fprintf(stderr, "  support_density       = %f\n", config.support_density);
	fprintf(stderr, "  support_flow_mult     = %f\n", config.support_flow_mult);
	fprintf(stderr, "  min_layer_time        = %f\n", config.min_layer_time);
	fprintf(stderr, "  layer_time_samples    = %d\n", config.layer_time_samples);
	fprintf(stderr, "  min_feed_rate         = %f\n", config.min_feed_rate);
	fprintf(stderr, "  brim_lines            = %d\n", config.brim_lines);
	fprintf(stderr, "  material_density      = %f\n", config.material_density);
	fprintf(stderr, "  material_cost         = %f\n", config.material_cost);
#ifdef _OPENMP
	fprintf(stderr, "OpenMP enabled\n");
#endif

	fprintf(stderr, "load object...\n");
	if (read_binary_stl(&o, path)) {
		fprintf(stderr, "error: failed to read stl: %s: %s\n", path, strerror(errno));
		return 1;
	}

	fprintf(stderr, "  polygons = %zd\n", o.n);
	fprintf(stderr, "  center   = (%f, %f, %f)\n", o.c.x, o.c.y, o.c.z);
	fprintf(stderr, "  height   = %f\n", o.h);
	fprintf(stderr, "  width    = %f\n", o.w);
	fprintf(stderr, "  depth    = %f\n", o.d);

	fprintf(stderr, "scale and translate object...\n");
	scale_object(&o, config.xy_scale_factor * scale_factor, config.xy_scale_factor * scale_factor, config.z_scale_factor * scale_factor);
	translate_object(&o, -o.c.x + config.x_center, -o.c.y + config.y_center, o.h / 2 - o.c.z - z_chop);
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
