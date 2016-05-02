#!/bin/sh

MACHINE=${MACHINE:-"ultra3d"}
EXTRUDER=${EXTRUDER:-"left"}
MATERIAL=${MATERIAL:-"inland_pla"}

if [ -z "$1" ]; then
	echo "usage: $0 binary_stl_file [shiv_options]"
	exit 1
fi

BASE_DIR="$(dirname "$0")"
CONFIG_DIR="$BASE_DIR/configs"
INFILE="$1"
OUTFILE="${INFILE%.stl}"
OUTFILE="${OUTFILE%.STL}"
OUTFILE="${OUTFILE}_s.gcode"
shift 1

time "$BASE_DIR/shiv" -c "$CONFIG_DIR/global" -c "$CONFIG_DIR/$MACHINE/$MACHINE" -c "$CONFIG_DIR/$MACHINE/$EXTRUDER" -c "$CONFIG_DIR/$MACHINE/$MATERIAL" -o "$OUTFILE" "$@" "$INFILE"
