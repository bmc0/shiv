#!/bin/sh

fail() {
	echo ">> Failed"
	exit $1
}

[ -z "$CXX" ] && CXX="g++"
CXXFLAGS="-O2 -std=gnu++11 -Wall -Wextra -fopenmp $CXXFLAGS"
LDFLAGS="-lm $LDFLAGS"

$CXX -o shiv $CXXFLAGS $LDFLAGS shiv.cpp clipper.cpp || fail $?

echo ">> Build successful"
