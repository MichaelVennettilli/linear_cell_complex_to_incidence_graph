#!/bin/sh
cgal_create_CMakeLists -c Qt5
if [ ! -d build ]; then
	echo DNE
	mkdir build;
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
rm ../CMakeLists.txt.bak
