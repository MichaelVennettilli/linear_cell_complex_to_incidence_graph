#!/bin/sh
cgal_create_CMakeLists -c Qt5
sed -i '1 s/^/set(CMAKE_CXX_STANDARD 17)\n/' CMakeLists.txt
#echo 'target_compile_features("main" PRIVATE cxx_std_17)' >>  CMakeLists.txt
if [ ! -d build ]; then
	echo DNE
	mkdir build;
fi
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make
rm ../CMakeLists.txt.bak
