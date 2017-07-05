#!/bin/sh

MACHINE=${MACHINE:-"ultra3d"}
EXTRUDER=${EXTRUDER:-"left"}
MATERIAL=${MATERIAL:-"inland_pla"}

if [ -z "$1" ]; then
	echo "usage: $0 binary_stl[.gz|.bz2|.xz|.lz4] [shiv_options]"
	exit 1
fi

BASE_DIR="$HOME/src/shiv"
CONFIG_DIR="$BASE_DIR/configs"
INFILE="$1"
INFILE_BASE="$(basename "$INFILE")"
INFILE_DIR="$(dirname "$INFILE")"
EXT="${INFILE_BASE##*.}"
OUTFILE="$INFILE_DIR/${INFILE_BASE%%.*}.gcode"
shift 1

run_shiv() {
	time "$BASE_DIR/shiv" \
		-S gcode_variable=date="$(date)" \
		-S gcode_variable=machine="$MACHINE" \
		-S gcode_variable=extruder="$EXTRUDER" \
		-S gcode_variable=material="$MATERIAL" \
		-c "$CONFIG_DIR/global" \
		-c "$CONFIG_DIR/$MACHINE/$MACHINE" \
		-c "$CONFIG_DIR/$MACHINE/$EXTRUDER" \
		-c "$CONFIG_DIR/$MACHINE/$MATERIAL" \
		"$@"
}

case "$EXT" in
	gz)  gzip  -dc "$INFILE" | run_shiv -o "$OUTFILE" "$@" - ;;
	bz2) bzip2 -dc "$INFILE" | run_shiv -o "$OUTFILE" "$@" - ;;
	xz)  xz    -dc "$INFILE" | run_shiv -o "$OUTFILE" "$@" - ;;
	lz4) lz4   -dc "$INFILE" | run_shiv -o "$OUTFILE" "$@" - ;;
	*)   run_shiv -o "$OUTFILE" "$@" "$INFILE" ;;
esac
