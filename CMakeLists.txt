cmake_minimum_required(VERSION 2.4.6)
Project("SDFGen")

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries

# SET(CMAKE_BUILD_TYPE Coverage)
SET(CMAKE_BUILD_TYPE Release)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#These flags might not work on every system, especially the release flags, comment out as needed
set(CMAKE_CXX_FLAGS "-O3")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -msse4.2 -mfpmath=sse -mtune=native -march=native")


add_executable(${PROJECT_NAME}
main.cpp  makelevelset3.cpp array1.h  array2.h  array3.h  hashgrid.h  hashtable.h makelevelset3.h  util.h  vec.h
)
