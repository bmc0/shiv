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
extern "C" {
	#include <getopt.h>
}
#ifdef _OPENMP
#include <omp.h>
#endif
#include "list.h"
#include "clipper.hpp"

#define SCALE_CONSTANT             1000  /* Clipper uses integers, so we need to scale floating point values. Precision is 1 / SCALE_CONSTANT units. */
#define CLEAN_DIST                 (1.41421356 * config.coarseness)
#define REMOVE_OVERLAPPING_SHELLS  true
#define CLIPPER_MITER_LIMIT        2.0
#define CLIPPER_ARC_TOLERANCE      10.0
#define CLIPPER_JOIN_TYPE          ClipperLib::jtSquare
#define USE_BOUNDING_BOX           1
#define SHIV_DEBUG                 1
typedef double fl_t;

#define MINIMUM(a, b) ((a < b) ? a : b)
#define MAXIMUM(a, b) ((a > b) ? a : b)
#define FL_T_TO_INTPOINT(x, y)  ClipperLib::IntPoint((ClipperLib::cInt) lround((x) * SCALE_CONSTANT), (ClipperLib::cInt) lround((y) * SCALE_CONSTANT))
#define INTPOINT_TO_FL_T(point) ((fl_t) point.X) / SCALE_CONSTANT, ((fl_t) point.Y) / SCALE_CONSTANT

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
	"            [-r roof_thickness] [-f floor_thickness] [-C coarseness]\n"
	"            [-x x_translate] [-y y_translate] binary_stl_file\n"
	"\n"
	"flags:\n"
	"  -h                    show this help\n"
	"  -p                    preview slices (pipe stdout to gnuplot)\n"
	"  -o output_path        output gcode path\n"
	"  -c config_path        configuration file path\n"
	"  -O option=value       set option to value\n"
	"  -l layer_height       layer height\n"
	"  -w extrusion_width    constrained (infill) extrusion width\n"
	"  -t tolerance          segment connection tolerance\n"
	"  -s scale_factor       object scale factor\n"
	"  -d infill_density     sparse infill density\n"
	"  -n shells             number of shells\n"
	"  -r roof_thickness     solid roof thickness\n"
	"  -f floor_thickness    solid floor thickness\n"
	"  -C coarseness         output coarseness\n"
	"  -x x_translate        translate object in the x-axis\n"
	"  -y y_translate        translate object in the y-axis\n";

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
	fl_t packing_density       = 0.99;     /* Solid packing density (should be slightly less than 1; 0.99 seems to work well for PLA) */
	fl_t infill_density        = 0.2;      /* Sparse infill density */
	int shells                 = 2;        /* Number of loops/perimeters/shells (whatever you want to call them) */
	fl_t roof_thickness        = 0.8;      /* Solid surface thickness when looking upwards */
	int roof_layers;                       /* Calculated from roof_thickness */
	fl_t floor_thickness       = 0.8;      /* Solid surface thickness when looking downwards */
	int floor_layers;                      /* Calculated from floor_thickness */
	fl_t material_diameter     = 1.75;     /* Diameter of material */
	fl_t material_area;                    /* Cross-sectional area of filament (calculated from material_diameter) */
	fl_t flow_multiplier       = 1.0;      /* Flow rate adjustment */
	fl_t perimeter_feed_rate   = 25.0;     /* Outer shell feed rate */
	fl_t loop_feed_rate        = 40.0;     /* Inner shell feed rate */
	fl_t infill_feed_rate      = 50.0;
	fl_t travel_feed_rate      = 120.0;
	fl_t first_layer_mult      = 0.5;      /* First layer feed rates (except travel) are multiplied by this value */
	fl_t retract_len           = 1.0;
	fl_t retract_speed         = 20.0;
	fl_t restart_speed         = -1.0;     /* -1 means same as retract_speed */
	fl_t retract_min_travel    = 1.6;
	fl_t retract_threshold     = 8.0;      /* Unconditional retraction threshold */
	fl_t retract_within_island = false;
	int cool_layer             = 1;        /* Turn on part cooling at this layer */
	char *start_gcode          = NULL;
	char *end_gcode            = NULL;
	char *cool_on_gcode        = NULL;
	char *cool_off_gcode       = NULL;
	fl_t temp                  = 220.0;    /* Hotend temperature */
	fl_t bed_temp              = 65.0;
	bool infill_first          = false;    /* Do infill before shells */
	bool align_seams           = true;     /* Align seams to the lower left corner */
	bool clean_insets          = true;     /* Do ClipperLib::CleanPolygon operation on all insets (only the initial outline is cleaned if this is false) */
	bool fill_inset_gaps       = true;     /* Fill gaps between shells */
	bool anchor                = true;     /* Clip and anchor inset paths */
	bool outside_first         = false;    /* Prefer exterior shells */
	fl_t fill_threshold        = 0.2;      /* Remove infill or inset gap fill when it would be narrower than extrusion_width * fill_threshold */
	fl_t min_layer_time        = 8.0;      /* Slow down if the estimated layer time is less than this value */
	int layer_time_samples     = 5;        /* Number of samples in the layer time moving average */
	fl_t min_feed_rate         = 10.0;
	fl_t material_density      = 0.00125;  /* Material density in <arbitrary mass unit> / <input/output unit>^3. The default is correct for PLA and millimeter input/output units */
	fl_t material_cost         = 0.01499;  /* Material cost in <arbitrary currency> / <arbitrary mass unit>. The arbitrary mass unit must be the same as used in material_density */
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
	struct cint_rect box;  /* bounding box */
};

struct g_move {
	fl_t x, y, z, e, feed_rate;
	bool scalable;
};

struct slice {
#ifdef _OPENMP
	volatile int lock;
	volatile ssize_t n_seg, s_len;
#else
	ssize_t n_seg, s_len;
#endif
	struct segment *s;
	std::vector<struct island> islands;
	std::vector<struct g_move> moves;
	fl_t layer_time;
};

struct machine {
	fl_t x, y, z, e;
	fl_t feed_rate;
	bool is_retracted, new_island;
};

#ifdef _OPENMP
static void spin_lock(volatile int *x)
{
	while (__sync_lock_test_and_set(x, 1))
		while (*x);
}

static void spin_unlock(volatile int *x)
{
	__sync_lock_release(x);
}
#endif

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
	else if (strcmp(key, "material_diameter") == 0) {
		config.material_diameter = atof(value);
		CHECK_VALUE(config.material_diameter > 0.0, "material diameter", "> 0");
	}
	else if (strcmp(key, "flow_multiplier") == 0) {
		config.flow_multiplier = atof(value);
		CHECK_VALUE(config.flow_multiplier >= 0.0, "flow multiplier", ">= 0");
	}
	else if (strcmp(key, "perimeter_feed_rate") == 0) {
		config.perimeter_feed_rate = atof(value);
		CHECK_VALUE(config.perimeter_feed_rate > 0.0, "perimeter feed rate", "> 0");
	}
	else if (strcmp(key, "loop_feed_rate") == 0) {
		config.loop_feed_rate = atof(value);
		CHECK_VALUE(config.loop_feed_rate > 0.0, "loop feed rate", "> 0");
	}
	else if (strcmp(key, "infill_feed_rate") == 0) {
		config.infill_feed_rate = atof(value);
		CHECK_VALUE(config.infill_feed_rate > 0.0, "infill feed rate", "> 0");
	}
	else if (strcmp(key, "travel_feed_rate") == 0) {
		config.travel_feed_rate = atof(value);
		CHECK_VALUE(config.travel_feed_rate > 0.0, "travel feed rate", "> 0");
	}
	else if (strcmp(key, "first_layer_mult") == 0) {
		config.first_layer_mult = atof(value);
		CHECK_VALUE(config.first_layer_mult > 0.0, "first layer multiplier", "> 0");
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
		CHECK_VALUE(config.restart_speed > 0.0 || config.restart_speed == -1.0, "restart_speed", "> 0 or -1");
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
	else if (strcmp(key, "anchor") == 0) {
		config.anchor = PARSE_BOOL(value);
	}
	else if (strcmp(key, "outside_first") == 0) {
		config.outside_first = PARSE_BOOL(value);
	}
	else if (strcmp(key, "fill_threshold") == 0) {
		config.fill_threshold = atof(value);
		CHECK_VALUE(config.fill_threshold >= 0.0, "fill threshold", ">= 0");
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
			while (next && (*(next + 1) == ' ' || *(next + 1) == '\t'))
				next = strchr(next + 1, '\n');
			next = isolate(next, '\n');
		}
		if (*key != '\n' && *key != '#') {
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
		spin_lock(&slices[i].lock);
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
		spin_unlock(&slices[i].lock);
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
			/* If we flipped more than half of the segments, the first segment in the outline was probably oriented incorrectly */
			if (flip_count > segment_count / 2) {
				DEBUG("reversed outline order at layer %zd\n", slice_index + 1);
				ClipperLib::ReversePath(poly);
			}
			outlines.push_back(poly);
			/* DEBUG("inexact_count = %zd / %zd\n", inexact_count, segment_count); */
		}
		next_poly:
		oseg.h = oseg.t = NULL;
	}
	free(slice->s);
	ClipperLib::PolyTree tree;
	ClipperLib::Clipper c;
	c.AddPaths(outlines, ClipperLib::ptSubject, true);
	c.Execute(ClipperLib::ctUnion, tree, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
	generate_islands(slice, &tree);
#if USE_BOUNDING_BOX
	for (struct island &island : slice->islands)
		find_bounding_box(&island);
#endif
}

static void do_offset(ClipperLib::Paths &src, ClipperLib::Paths &dest, fl_t dist, bool remove_overlap)
{
	ClipperLib::ClipperOffset co(CLIPPER_MITER_LIMIT, CLIPPER_ARC_TOLERANCE);
	co.AddPaths(src, CLIPPER_JOIN_TYPE, ClipperLib::etClosedPolygon);
	if (remove_overlap) {
		fl_t extra = (dist > 0) ? (config.extrusion_width / 2.0) : (config.extrusion_width / -2.0);
		co.Execute(dest, (dist + extra) * SCALE_CONSTANT);
		co.Clear();
		co.AddPaths(dest, CLIPPER_JOIN_TYPE, ClipperLib::etClosedPolygon);
		co.Execute(dest, -extra * SCALE_CONSTANT);
	}
	else
		co.Execute(dest, dist * SCALE_CONSTANT);
}

static void generate_insets(struct slice *slice)
{
	for (struct island &island : slice->islands) {
		island.insets = (ClipperLib::Paths *) calloc(config.shells, sizeof(ClipperLib::Paths));
		if (!island.insets)
			die(e_nomem, 2);
		fl_t edge_offset = -config.edge_width / 2.0;
		if (config.shells > 0) {
			do_offset(island.outlines, island.insets[0], edge_offset, false);
			if (config.clean_insets)
				ClipperLib::CleanPolygons(island.insets[0], CLEAN_DIST);
			for (int i = 1; i < config.shells; ++i) {
				do_offset(island.insets[i - 1], island.insets[i], -config.extrusion_width, REMOVE_OVERLAPPING_SHELLS);
				if (config.clean_insets)
					ClipperLib::CleanPolygons(island.insets[i], CLEAN_DIST);
				if (island.insets[i].size() == 0)  /* break if nothing is being generated */
					goto done;
			}
			do_offset(island.insets[config.shells - 1], island.infill_insets, -config.extrusion_width / 2.0, false);
		}
		else {
			do_offset(island.outlines, island.infill_insets, edge_offset, false);
		}
		ClipperLib::CleanPolygons(island.infill_insets, CLEAN_DIST * 4.0);

		done:
		if (config.shells > 1 && config.fill_inset_gaps) {
			ClipperLib::ClipperOffset co(CLIPPER_MITER_LIMIT, CLIPPER_ARC_TOLERANCE);
			ClipperLib::Paths hole;
			island.inset_gaps = (ClipperLib::Paths *) calloc(config.shells - 1, sizeof(ClipperLib::Paths));
			if (!island.inset_gaps)
				die(e_nomem, 2);
			for (int i = 0; i < config.shells - 1 && island.insets[i].size() > 0; ++i) {
				co.AddPaths(island.insets[i], CLIPPER_JOIN_TYPE, ClipperLib::etClosedPolygon);
				hole = island.insets[i + 1];
				ClipperLib::ReversePaths(hole);
				co.AddPaths(hole, CLIPPER_JOIN_TYPE, ClipperLib::etClosedPolygon);
				if (config.fill_threshold > 0.0) {
					co.Execute(island.inset_gaps[i], (config.extrusion_width + config.extrusion_width * config.fill_threshold) / -2.0 * SCALE_CONSTANT);
					co.Clear();
					co.AddPaths(island.inset_gaps[i], CLIPPER_JOIN_TYPE, ClipperLib::etClosedPolygon);
					co.Execute(island.inset_gaps[i], config.extrusion_width * config.fill_threshold / 2.0 * SCALE_CONSTANT);
				}
				else {
					co.Execute(island.inset_gaps[i], -config.extrusion_width / 2.0 * SCALE_CONSTANT);
				}
				co.Clear();
			}
		}
		if (config.align_seams) {
			for (int i = 0; i < config.shells; ++i) {
				for (ClipperLib::Path &p : island.insets[i]) {
					if (p.size() > 1) {
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

static void generate_infill_patterns(struct object *o)
{
	ClipperLib::Path line(2);
	const ClipperLib::cInt move = lround(config.extrusion_width * sqrt(2.0) * SCALE_CONSTANT);
	const ClipperLib::cInt min_x = lround((-o->w / 2.0 + o->c.x) * SCALE_CONSTANT);
	const ClipperLib::cInt min_y = lround((-o->d / 2.0 + o->c.y) * SCALE_CONSTANT);
	const ClipperLib::cInt max_x = lround((o->w / 2.0 + o->c.x) * SCALE_CONSTANT);
	const ClipperLib::cInt max_y = lround((o->d / 2.0 + o->c.y) * SCALE_CONSTANT);

	/* generate even solid fill */
	line[0].X = min_x;
	line[0].Y = min_y + move;
	line[1].X = min_x + move;
	line[1].Y = min_y;
	while (line[0].X < max_x && line[1].Y < max_y) {
		o->solid_infill_patterns[0].push_back(line);
		if (line[0].Y >= max_y)
			line[0].X += move;
		else
			line[0].Y += move;

		if (line[1].X >= max_x)
			line[1].Y += move;
		else
			line[1].X += move;
	}

	/* generate odd solid fill */
	line[0].X = min_x;
	line[0].Y = max_y - move;
	line[1].X = min_x + move;
	line[1].Y = max_y;
	while (line[0].X < max_x && line[1].Y > min_y) {
		o->solid_infill_patterns[1].push_back(line);
		if (line[0].Y <= min_y)
			line[0].X += move;
		else
			line[0].Y -= move;

		if (line[1].X >= max_x)
			line[1].Y -= move;
		else
			line[1].X += move;
	}

	if (config.infill_density > 0.0) {
		const ClipperLib::cInt sparse_move = lround(config.extrusion_width * sqrt(2.0) / config.infill_density * 2.0 * SCALE_CONSTANT);
		line[0].X = min_x;
		line[0].Y = min_y + sparse_move;
		line[1].X = min_x + sparse_move;
		line[1].Y = min_y;
		while (line[0].X < max_x && line[1].Y < max_y) {
			o->sparse_infill_pattern.push_back(line);
			if (line[0].Y >= max_y)
				line[0].X += sparse_move;
			else
				line[0].Y += sparse_move;

			if (line[1].X >= max_x)
				line[1].Y += sparse_move;
			else
				line[1].X += sparse_move;
		}
		line[0].X = min_x;
		line[0].Y = max_y - sparse_move;
		line[1].X = min_x + sparse_move;
		line[1].Y = max_y;
		while (line[0].X < max_x && line[1].Y > min_y) {
			o->sparse_infill_pattern.push_back(line);
			if (line[0].Y <= min_y)
				line[0].X += sparse_move;
			else
				line[0].Y -= sparse_move;

			if (line[1].X >= max_x)
				line[1].Y -= sparse_move;
			else
				line[1].X += sparse_move;
		}
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
		if (config.infill_density == 1.0 || slice_index < config.floor_layers || slice_index + config.roof_layers >= o->n_slices) {
			if (config.fill_threshold > 0.0) {
				do_offset(island.infill_insets, island.infill_insets, -config.extrusion_width * config.fill_threshold / 2.0, false);
				do_offset(island.infill_insets, island.infill_insets, config.extrusion_width * config.fill_threshold / 2.0, false);
			}
			c.AddPaths(o->solid_infill_patterns[slice_index % 2], ClipperLib::ptSubject, false);
			c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
			if (config.fill_inset_gaps)
				for (int i = 0; i < config.shells - 1; ++i)
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
		}
		else if (config.floor_layers > 0 || config.roof_layers > 0) {
			ClipperLib::Clipper sc;
			ClipperLib::Paths s_tmp;
			
			sc.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			for (int i = -config.floor_layers; i <= config.roof_layers; ++i) {
				if (i != 0) {
					if (config.shells > 0) {
						for (struct island &clip_island : o->slices[slice_index + i].islands)
							if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
								sc.AddPaths(clip_island.insets[config.shells - 1], ClipperLib::ptClip, true);
					}
					else {
						for (struct island &clip_island : o->slices[slice_index + i].islands)
							if (BOUNDING_BOX_INTERSECTS(island.box, clip_island.box))
								sc.AddPaths(clip_island.infill_insets, ClipperLib::ptClip, true);
					}
					sc.Execute(ClipperLib::ctIntersection, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
					sc.Clear();
					if (i != config.roof_layers)
						sc.AddPaths(s_tmp, ClipperLib::ptSubject, true);
				}
			}
			sc.AddPaths(s_tmp, ClipperLib::ptClip, true);
			sc.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			sc.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			sc.Clear();
			if (config.fill_threshold > 0.0) {
				do_offset(s_tmp, s_tmp, -config.extrusion_width * config.fill_threshold / 2.0, false);
				do_offset(s_tmp, s_tmp, config.extrusion_width * config.fill_threshold / 2.0, false);
			}
			c.AddPaths(o->solid_infill_patterns[slice_index % 2], ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			if (config.fill_inset_gaps)
				for (int i = 0; i < config.shells - 1; ++i)
					c.AddPaths(island.inset_gaps[i], ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.solid_infill);
			c.Clear();

			sc.AddPaths(s_tmp, ClipperLib::ptClip, true);
			sc.AddPaths(island.infill_insets, ClipperLib::ptSubject, true);
			sc.Execute(ClipperLib::ctDifference, s_tmp, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			if (config.fill_threshold > 0.0) {
				do_offset(s_tmp, s_tmp, -config.extrusion_width * config.fill_threshold / 2.0, false);
				do_offset(s_tmp, s_tmp, config.extrusion_width * config.fill_threshold / 2.0, false);
			}
			c.AddPaths(o->sparse_infill_pattern, ClipperLib::ptSubject, false);
			c.AddPaths(s_tmp, ClipperLib::ptClip, true);
			c.Execute(ClipperLib::ctIntersection, s, ClipperLib::pftNonZero, ClipperLib::pftNonZero);
			ClipperLib::OpenPathsFromPolyTree(s, island.sparse_infill);
		}
		else {
			if (config.fill_threshold > 0.0) {
				do_offset(island.infill_insets, island.infill_insets, -config.extrusion_width * config.fill_threshold / 2.0, false);
				do_offset(island.infill_insets, island.infill_insets, config.extrusion_width * config.fill_threshold / 2.0, false);
			}
			c.AddPaths(o->sparse_infill_pattern, ClipperLib::ptSubject, false);
			c.AddPaths(island.infill_insets, ClipperLib::ptClip, true);
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

	fprintf(stderr, "sliced in %fs\n",
		(double) std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now() - start).count() / 1000000.0);
}

static void preview_slices(struct object *o)
{
	ssize_t i;
	fprintf(stdout, "set size ratio -1\n");
	fprintf(stdout, "set xrange [%e:%e]\n", o->c.x - o->w / 2.0, o->c.x + o->w / 2.0);
	fprintf(stdout, "set yrange [%e:%e]\n", o->c.y - o->d / 2.0, o->c.y + o->d / 2.0);
	for (i = 0; i < o->n_slices; ++i) {
		fprintf(stderr, "layer %zd/%zd: intersections = %zd; islands = %zd\n", i + 1, o->n_slices, o->slices[i].n_seg, o->slices[i].islands.size());
		fprintf(stdout, "set title 'Layer %zd/%zd'\n", i + 1, o->n_slices);
		/* Draw outlines */
		fprintf(stdout, "plot \"-\" lc \"red\" dt 3 with lines title \"outline\", \"-\" lc \"blue\" with lines title \"insets\", \"-\" lc \"green\" with lines title \"infill\"\n");
		for (struct island &island : o->slices[i].islands) {
			for (ClipperLib::Path &path : island.outlines) {
				if (path.size() > 0) {
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
		fprintf(stdout, "e\n");
		/* Draw insets */
		for (struct island &island : o->slices[i].islands) {
			for (int k = 0; k < config.shells; ++k) {
				for (ClipperLib::Path &path : island.insets[k]) {
					if (path.size() > 0) {
						for (ClipperLib::IntPoint &p : path)
							fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
						fprintf(stdout, "%.4e %.4e\n", ((double) path[0].X) / SCALE_CONSTANT, ((double) path[0].Y) / SCALE_CONSTANT);
						putc('\n', stdout);
					}
				}
			}
		}
		fprintf(stdout, "e\n");
		/* Draw infill */
		for (struct island &island : o->slices[i].islands) {
			for (ClipperLib::Path &path : island.solid_infill) {
				if (path.size() > 0) {
					for (ClipperLib::IntPoint &p : path)
						fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
					putc('\n', stdout);
				}
			}
			for (ClipperLib::Path &path : island.sparse_infill) {
				if (path.size() > 0) {
					for (ClipperLib::IntPoint &p : path)
						fprintf(stdout, "%.4e %.4e\n", ((double) p.X) / SCALE_CONSTANT, ((double) p.Y) / SCALE_CONSTANT);
					putc('\n', stdout);
				}
			}
		}
		fprintf(stdout, "e\n");
		fflush(stdout);
		fprintf(stderr, "press enter for next layer...");
		getchar();
	}
}

static bool ccw(ClipperLib::IntPoint &a, ClipperLib::IntPoint &b, ClipperLib::IntPoint &c)
{
	return (c.Y - a.Y) * (b.X - a.X) > (b.Y - a.Y) * (c.X - a.X);
}

static bool intersects(ClipperLib::IntPoint &a, ClipperLib::IntPoint &b, ClipperLib::IntPoint &c, ClipperLib::IntPoint &d)
{
	return ccw(a, c, d) != ccw(b, c, d) && ccw(a, b, c) != ccw(a, b, d);
}

static bool crosses_boundary(struct machine *m, struct island *island, fl_t x, fl_t y)
{
	ClipperLib::IntPoint c((ClipperLib::cInt) (m->x * SCALE_CONSTANT), (ClipperLib::cInt) (m->y * SCALE_CONSTANT));
	ClipperLib::IntPoint d((ClipperLib::cInt) (x * SCALE_CONSTANT), (ClipperLib::cInt) (y * SCALE_CONSTANT));
	for (ClipperLib::Path &p : island->outlines) {
		ClipperLib::IntPoint a = p[0];
		for (size_t i = 1; i < p.size(); ++i) {
			ClipperLib::IntPoint b = p[i];
			if (intersects(a, b, c, d))
				return true;
			a = b;
		}
	}
	return false;
}

static void append_g_move(struct slice *slice, struct g_move &move, fl_t len)
{
	/* FIXME: should probably take acceleration into account... */
	slice->layer_time += len / move.feed_rate;
	slice->moves.push_back(move);
}

static void linear_move(struct slice *slice, struct island *island, struct machine *m, fl_t x, fl_t y, fl_t z, fl_t extra_e_len, fl_t feed_rate, bool is_travel)
{
	struct g_move move = { x, y, z, 0.0, feed_rate, !is_travel };
	fl_t e_len = 0.0;
	fl_t len = sqrt((m->x - x) * (m->x - x) + (m->y - y) * (m->y - y) + (m->z - z) * (m->z - z));
	if (is_travel) {
		if (!m->is_retracted && config.retract_len > 0.0 && (m->new_island || len > config.retract_threshold || (config.retract_within_island && len > config.retract_min_travel) || crosses_boundary(m, island, x, y))) {
			struct g_move retract_move = { m->x, m->y, m->z, -config.retract_len, config.retract_speed, false };
			append_g_move(slice, retract_move, config.retract_len);
			m->e -= config.retract_len;
			m->is_retracted = true;
		}
		m->new_island = false;
	}
	else {
		if (m->is_retracted) {
			struct g_move restart_move = { m->x, m->y, m->z, config.retract_len, config.restart_speed, false };
			append_g_move(slice, restart_move, config.retract_len);
			m->e += config.retract_len;
			m->is_retracted = false;
		}
		move.e = len * config.extrusion_area * config.flow_multiplier / config.material_area;
	}
	if (extra_e_len != 0.0) {
		struct g_move restart_move = { m->x, m->y, m->z, extra_e_len, feed_rate * config.extrusion_area / config.material_area, false };
		append_g_move(slice, restart_move, extra_e_len);
		m->e += extra_e_len;
	}
	if (x != m->x || y != m->y || z != m->z || e_len != 0) {
		append_g_move(slice, move, len);
		m->x = x;
		m->y = y;
		m->z = z;
		m->e += e_len;
	}
}

static bool path_len_is_greater_than(ClipperLib::Path &p, fl_t len)
{
	fl_t l = 0.0;
	fl_t x0 = ((fl_t) p[0].X) / SCALE_CONSTANT;
	fl_t y0 = ((fl_t) p[0].Y) / SCALE_CONSTANT;
	for (size_t i = 1; i < p.size(); ++i) {
		fl_t x1 = ((fl_t) p[i].X) / SCALE_CONSTANT;
		fl_t y1 = ((fl_t) p[i].Y) / SCALE_CONSTANT;
		l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
		if (l > len)
			return true;
		x0 = x1;
		y0 = y1;
	}
	fl_t x1 = ((fl_t) p[0].X) / SCALE_CONSTANT;
	fl_t y1 = ((fl_t) p[0].Y) / SCALE_CONSTANT;
	l += sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
	if (l > len)
		return true;
	return false;
}

/* Returns the new start index. Do not call this function on a path that has a total length of less than clip. */
static size_t clip_path(ClipperLib::Path &p, size_t start_idx, fl_t clip)
{
	fl_t l = 0.0;
	fl_t x0 = ((fl_t) p[start_idx].X) / SCALE_CONSTANT;
	fl_t y0 = ((fl_t) p[start_idx].Y) / SCALE_CONSTANT;
	/* We want to leave the original start point alone because it becomes the end point */
	start_idx = (start_idx + 1 == p.size()) ? 0 : start_idx + 1;
	for (;;) {
		fl_t x1 = ((fl_t) p[start_idx].X) / SCALE_CONSTANT;
		fl_t y1 = ((fl_t) p[start_idx].Y) / SCALE_CONSTANT;
		fl_t xv = x1 - x0;
		fl_t yv = y1 - y0;
		fl_t norm = sqrt(xv * xv + yv * yv);
		l += norm;
		if (l == clip)
			break;
		else if (l > clip) {
			fl_t new_x = x1 - (l - clip) * (xv / norm);
			fl_t new_y = y1 - (l - clip) * (yv / norm);
			p.insert(p.begin() + start_idx, FL_T_TO_INTPOINT(new_x, new_y));
			break;
		}
		x0 = x1;
		y0 = y1;
		p.erase(p.begin() + start_idx);
		if (start_idx == p.size())
			start_idx = 0;
	}
	return start_idx;
}

static void write_insets(struct slice *slice, struct island *island, struct machine *m, fl_t z, bool outside_first)
{
	if (config.shells > 0) {
		int inset_num;
		do {
			auto best = island->insets[0].begin();
			fl_t best_dist = HUGE_VAL;
			inset_num = -1;
			size_t start_idx = 0;
			for (int i = 0; i < config.shells; ++i) {
				if (!island->insets[i].empty()) {
					if (config.align_seams) {
						for (auto it = island->insets[i].begin(); it != island->insets[i].end(); ++it) {
							fl_t x = ((fl_t) (*it)[0].X) / SCALE_CONSTANT;
							fl_t y = ((fl_t) (*it)[0].Y) / SCALE_CONSTANT;
							fl_t dist = (x - m->x) * (x - m->x) + (y - m->y) * (y - m->y);
							if (outside_first) {
								if (i != 0)
									dist = dist * 2.0 * (i + 1) + 10.0;  /* prefer exterior */
							}
							else {
								if (i != config.shells - 1)
									dist = dist * 2.0 * (config.shells - i) + 10.0;  /* prefer interior */
							}
							if (dist < best_dist) {
								best_dist = dist;
								best = it;
								start_idx = 0;
								inset_num = i;
							}
						}
					}
					else {
						for (auto it = island->insets[i].begin(); it != island->insets[i].end(); ++it) {
							for (size_t k = 0; k < (*it).size(); ++k) {
								fl_t x = ((fl_t) (*it)[k].X) / SCALE_CONSTANT;
								fl_t y = ((fl_t) (*it)[k].Y) / SCALE_CONSTANT;
								fl_t dist = (x - m->x) * (x - m->x) + (y - m->y) * (y - m->y);
								if (outside_first) {
									if (i != 0)
										dist = dist * 2.0 * (i + 1) + 10.0;  /* prefer exterior */
								}
								else {
									if (i != config.shells - 1)
										dist = dist * 2.0 * (config.shells - i) + 10.0;  /* prefer interior */
								}
								if (dist < best_dist) {
									best_dist = dist;
									best = it;
									start_idx = k;
									inset_num = i;
								}
							}
						}
					}
				}
			}
			if (inset_num != -1) {
				ClipperLib::Path &p = *best;
				if (p.size() > 1) {
					bool first_point = true;
					fl_t feed_rate = (inset_num == 0) ? config.perimeter_feed_rate : config.loop_feed_rate;
					bool do_anchor = false, did_anchor = false;
					if (config.anchor && path_len_is_greater_than(p, config.extrusion_width * 3.0)) {
						start_idx = clip_path(p, start_idx, config.extrusion_width / 2.0);
						do_anchor = true;
					}
					size_t k = start_idx;
					do {
						if (first_point) {
							linear_move(slice, island, m, INTPOINT_TO_FL_T(p[k]), z, 0.0, config.travel_feed_rate, true);
						}
						else {
							fl_t anchor_e_len = 0.0;
							if (do_anchor && !did_anchor) {
								anchor_e_len = config.extrusion_width / 2.0 * config.extrusion_area * config.flow_multiplier / config.material_area;
								did_anchor = true;
							}
							linear_move(slice, island, m, INTPOINT_TO_FL_T(p[k]), z, anchor_e_len, feed_rate, false);
						}
						first_point = false;
						k = (k + 1 == p.size()) ? 0 : k + 1;
					} while (k != start_idx);
					if (!do_anchor)  /* If the path was not clipped, the start and end points are the same, so we need to do one more extrusion move back to the starting point. */
						linear_move(slice, island, m, INTPOINT_TO_FL_T(p[start_idx]), z, 0.0, feed_rate, false);
				}
				island->insets[inset_num].erase(best);
			}
		} while (inset_num != -1);
	}
}

static void write_infill(ClipperLib::Paths &paths, struct slice *slice, struct island *island, struct machine *m, fl_t z)
{
	while (!paths.empty()) {
		auto best = paths.end();
		fl_t best_dist = HUGE_VAL;
		bool flip_points = false;
		for (auto it = paths.begin(); it != paths.end(); ++it) {
			fl_t x0 = ((fl_t) (*it)[0].X) / SCALE_CONSTANT;
			fl_t y0 = ((fl_t) (*it)[0].Y) / SCALE_CONSTANT;
			fl_t x1 = ((fl_t) (*it)[1].X) / SCALE_CONSTANT;
			fl_t y1 = ((fl_t) (*it)[1].Y) / SCALE_CONSTANT;
			fl_t dist0 = (x0 - m->x) * (x0 - m->x) + (y0 - m->y) * (y0 - m->y);
			fl_t dist1 = (x1 - m->x) * (x1 - m->x) + (y1 - m->y) * (y1 - m->y);
			fl_t dist = MINIMUM(dist0, dist1);
			if (dist < best_dist) {
				flip_points = (dist1 < dist0);
				best_dist = dist;
				best = it;
			}
		}
		if (best != paths.end()) {
			ClipperLib::Path &p = *best;
			if (flip_points) {
				linear_move(slice, island, m, INTPOINT_TO_FL_T(p[1]), z, 0.0, config.travel_feed_rate, true);
				linear_move(slice, island, m, INTPOINT_TO_FL_T(p[0]), z, 0.0, config.infill_feed_rate, false);
			}
			else {
				linear_move(slice, island, m, INTPOINT_TO_FL_T(p[0]), z, 0.0, config.travel_feed_rate, true);
				linear_move(slice, island, m, INTPOINT_TO_FL_T(p[1]), z, 0.0, config.infill_feed_rate, false);
			}
		}
		if (best != paths.end())
			paths.erase(best);
	}
}

static void plan_moves(struct object *o, struct slice *slice, ssize_t layer_num)
{
	struct machine m = {};
	m.x = o->c.x - o->w / 2.0;
	m.y = o->c.y - o->d / 2.0;
	m.z = ((fl_t) layer_num) * config.layer_height + config.layer_height;
	m.is_retracted = true;  /* NOTE: This will cause unconditional retraction on layer change */
	while (slice->islands.size() > 0) {
		auto best = slice->islands.begin();
		fl_t best_dist = HUGE_VAL;
		for (auto it = slice->islands.begin(); it != slice->islands.end(); ++it) {
			if (config.align_seams) {
				for (ClipperLib::Path &p : it->outlines) {
					fl_t x = ((fl_t) p[0].X) / SCALE_CONSTANT;
					fl_t y = ((fl_t) p[0].Y) / SCALE_CONSTANT;
					fl_t dist = (x - m.x) * (x - m.x) + (y - m.y) * (y - m.y);
					if (dist < best_dist) {
						best_dist = dist;
						best = it;
					}
				}
			}
			else {
				for (ClipperLib::Path &p : it->outlines) {
					for (size_t i = 0; i < p.size(); ++i) {
						fl_t x = ((fl_t) p[i].X) / SCALE_CONSTANT;
						fl_t y = ((fl_t) p[i].Y) / SCALE_CONSTANT;
						fl_t dist = (x - m.x) * (x - m.x) + (y - m.y) * (y - m.y);
						if (dist < best_dist) {
							best_dist = dist;
							best = it;
						}
					}
				}
			}
		}
		struct island &island = *best;
		m.new_island = true;
		/* Concat sparse infill into solid infill vector so the infill pathing is planned together */
		island.solid_infill.insert(island.solid_infill.end(), island.sparse_infill.begin(), island.sparse_infill.end());
		island.sparse_infill.clear();
		if (config.infill_first && layer_num != 0) {
			write_infill(island.solid_infill, slice, &island, &m, m.z);
			write_insets(slice, &island, &m, m.z, config.outside_first || layer_num == 0);
		}
		else {
			write_insets(slice, &island, &m, m.z, config.outside_first || layer_num == 0);
			write_infill(island.solid_infill, slice, &island, &m, m.z);
		}
		free(island.insets);
		free(island.inset_gaps);
		slice->islands.erase(best);
	}
	/* We need to retract at the end of each layer */
	if (!m.is_retracted) {
		struct g_move retract_move = { m.x, m.y, m.z, -config.retract_len, config.retract_speed, false };
		append_g_move(slice, retract_move, config.retract_len);
	}
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

static void write_gcode_move(FILE *f, struct g_move *move, struct machine *m, fl_t feed_rate_mult)
{
	fl_t feed_rate = move->feed_rate;
	if (move->scalable) {
		feed_rate *= feed_rate_mult;
		if (feed_rate < config.min_feed_rate)
			feed_rate = config.min_feed_rate;
	}
	fputs("G1", f);
	if (move->x != m->x)
		fprintf(f, " X%.3f", move->x);
	if (move->y != m->y)
		fprintf(f, " Y%.3f", move->y);
	if (move->z != m->z)
		fprintf(f, " Z%.3f", move->z);
	if (move->e != 0.0)
		fprintf(f, " E%.5f", m->e + move->e);
	if (feed_rate != m->feed_rate)
		fprintf(f, " F%ld", lround(feed_rate * 60.0));
	fputc('\n', f);
	m->x = move->x;
	m->y = move->y;
	m->z = move->z;
	m->e += move->e;
	m->feed_rate = feed_rate;
}

static int write_gcode(const char *path, struct object *o)
{
	fputs("plan moves...", stderr);
#ifdef _OPENMP
	#pragma omp parallel for schedule(dynamic)
#endif
	for (ssize_t i = 0; i < o->n_slices; ++i)
		plan_moves(o, &o->slices[i], i);
	fputs(" done\n", stderr);

	FILE *f;
	if (strcmp(path, "-") == 0)
		f = stdout;
	else
		f = fopen(path, "w");
	if (!f)
		return 1;
	struct machine m = {};
	m.x = o->c.x - o->w / 2.0;
	m.y = o->c.y - o->d / 2.0;
	fl_t feed_rate_mult = config.first_layer_mult;
	fprintf(stderr, "write gcode to %s...", path);
	write_gcode_string(config.start_gcode, f);
	for (ssize_t i = 0; i < o->n_slices; ++i) {
		struct slice *slice = &o->slices[i];
		fprintf(f, "; layer %zd (z = %f)\n", i, ((fl_t) i) * config.layer_height + config.layer_height);
		if (i == config.cool_layer)
			write_gcode_string(config.cool_on_gcode, f);
		fl_t average_layer_time = slice->layer_time;
		for (int k = 1; k < config.layer_time_samples; ++k)
			average_layer_time += (k < i) ? o->slices[i - k].layer_time : o->slices[0].layer_time;
		average_layer_time /= config.layer_time_samples;
		if (average_layer_time < config.min_layer_time)
			feed_rate_mult = average_layer_time / config.min_layer_time;
		for (struct g_move &move : o->slices[i].moves)
			write_gcode_move(f, &move, &m, feed_rate_mult);
		feed_rate_mult = 1.0;
	}
	write_gcode_string(config.cool_off_gcode, f);
	write_gcode_string(config.end_gcode, f);
	fputs(" done\n", stderr);
	fl_t mass = config.material_area * m.e * config.material_density / config.flow_multiplier;
	fprintf(f, "; material length = %.4f\n", m.e / config.flow_multiplier);
	fprintf(f, "; material mass   = %.4f\n", mass);
	fprintf(f, "; material cost   = %.4f\n", mass * config.material_cost);
	fprintf(stderr, "material length = %.4f\n", m.e / config.flow_multiplier);
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

int main(int argc, char *argv[])
{
	int opt;
	char *path, *output_path = NULL;
	struct object o;
	fl_t scale_factor = 1.0, x_translate = 0.0, y_translate = 0.0;
	bool do_preview = false;

	/* Parse options */
	while ((opt = getopt(argc, argv, ":hpo:c:O:l:w:t:s:d:n:r:f:C:x:y:")) != -1) {
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
	config.edge_width = (config.extrusion_area - (config.layer_height * config.layer_height * M_PI / 4)) / config.layer_height + config.layer_height;
	config.material_area = M_PI * config.material_diameter * config.material_diameter / 4.0;
	if (config.restart_speed == -1.0)
		config.restart_speed = config.retract_speed;
	if (config.cool_on_gcode == NULL)
		config.cool_on_gcode = strdup(DEFAULT_COOL_ON_STR);
	if (config.cool_off_gcode == NULL)
		config.cool_off_gcode = strdup(DEFAULT_COOL_OFF_STR);
	config.x_center += x_translate;
	config.y_center += y_translate;

	fprintf(stderr, "configuration:\n");
	fprintf(stderr, "  layer height          = %f\n", config.layer_height);
	fprintf(stderr, "  tolerance             = %f\n", sqrt(config.tolerance));
	fprintf(stderr, "  output coarseness     = %f (%f units)\n", config.coarseness, config.coarseness / SCALE_CONSTANT);
	fprintf(stderr, "  extrusion width       = %f\n", config.extrusion_width);
	fprintf(stderr, "  edge width            = %f\n", config.edge_width);
	fprintf(stderr, "  extrusion area        = %f\n", config.extrusion_area);
	fprintf(stderr, "  object scale factor   = %f\n", scale_factor);
	fprintf(stderr, "  xy scale factor       = %f\n", config.xy_scale_factor);
	fprintf(stderr, "  z scale factor        = %f\n", config.z_scale_factor);
	fprintf(stderr, "  x center              = %f\n", config.x_center);
	fprintf(stderr, "  y center              = %f\n", config.x_center);
	fprintf(stderr, "  packing density       = %f\n", config.packing_density);
	fprintf(stderr, "  infill density        = %f\n", config.infill_density);
	fprintf(stderr, "  shells                = %d\n", config.shells);
	fprintf(stderr, "  roof layers           = %d\n", config.roof_layers);
	fprintf(stderr, "  floor layers          = %d\n", config.floor_layers);
	fprintf(stderr, "  material diameter     = %f\n", config.material_diameter);
	fprintf(stderr, "  material area         = %f\n", config.material_area);
	fprintf(stderr, "  flow multiplier       = %f\n", config.flow_multiplier);
	fprintf(stderr, "  perimeter feed rate   = %f\n", config.perimeter_feed_rate);
	fprintf(stderr, "  loop feed rate        = %f\n", config.loop_feed_rate);
	fprintf(stderr, "  infill feed rate      = %f\n", config.infill_feed_rate);
	fprintf(stderr, "  travel feed rate      = %f\n", config.travel_feed_rate);
	fprintf(stderr, "  first layer mult      = %f\n", config.first_layer_mult);
	fprintf(stderr, "  retract len           = %f\n", config.retract_len);
	fprintf(stderr, "  retract speed         = %f\n", config.retract_speed);
	fprintf(stderr, "  restart speed         = %f\n", config.restart_speed);
	fprintf(stderr, "  retract min travel    = %f\n", config.retract_min_travel);
	fprintf(stderr, "  retract threshold     = %f\n", config.retract_threshold);
	fprintf(stderr, "  retract within island = %s\n", (config.retract_within_island) ? "true" : "false");
	fprintf(stderr, "  cool layer            = %d\n", config.cool_layer);
	fprintf(stderr, "  temp                  = %f\n", config.temp);
	fprintf(stderr, "  bed temp              = %f\n", config.bed_temp);
	fprintf(stderr, "  infill first          = %s\n", (config.infill_first) ? "true" : "false");
	fprintf(stderr, "  align seams           = %s\n", (config.align_seams) ? "true" : "false");
	fprintf(stderr, "  clean insets          = %s\n", (config.clean_insets) ? "true" : "false");
	fprintf(stderr, "  fill inset gaps       = %s\n", (config.fill_inset_gaps) ? "true" : "false");
	fprintf(stderr, "  anchor                = %s\n", (config.anchor) ? "true" : "false");
	fprintf(stderr, "  outside first         = %s\n", (config.outside_first) ? "true" : "false");
	fprintf(stderr, "  fill threshold        = %f\n", config.fill_threshold);
	fprintf(stderr, "  min layer time        = %f\n", config.min_layer_time);
	fprintf(stderr, "  layer time samples    = %d\n", config.layer_time_samples);
	fprintf(stderr, "  min feed rate         = %f\n", config.min_feed_rate);
	fprintf(stderr, "  material density      = %f\n", config.material_density);
	fprintf(stderr, "  material cost         = %f\n", config.material_cost);
#ifdef _OPENMP
	fprintf(stderr, "  openmp                = yes\n");
#else
	fprintf(stderr, "  openmp                = no\n");
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
	translate_object(&o, -o.c.x + config.x_center, -o.c.y + config.y_center, o.h / 2 - o.c.z);
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
