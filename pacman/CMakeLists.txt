cmake_minimum_required(VERSION 3.9.0)

#set (CMAKE_CXX_STANDARD 11)		# needed for root 6
set( CMAKE_EXPORT_COMPILE_COMMANDS ON ) # generate database for vim autocompletion 

project(xephyr)
#########  IMPORTANT  ######################
# make sure the compiler supports c++11 libraries
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.8.1")
		message(FATAL_ERROR "GNU-gcc version less than 4.8.1, this does not support c++11 library needed for ROOT 6")
	endif()
endif()
##########################################

################ ROOT setup #####################
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED )
include(${ROOT_USE_FILE})
################################################

# dir of headers (which is the same as source)
include_directories(Xephyr/src)

# Get all source file and saves it in a list called SOURCES
file(GLOB SOURCES "Xephyr/src/*.cxx")

# Producing a shared library with all xephyr stuff
add_library(xelib SHARED ${SOURCES})
target_link_libraries(xelib ${ROOT_LIBRARIES}) # link the xephyr lib to the ROOT libraries


############  USER LIBRARIES AND EXE  #####################
