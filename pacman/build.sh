#!/bin/bash

cp Xephyr/pacman/CMakeLists.txt .
{ #try
	python3.4 Xephyr/pacman/check.py
} || {
	python3 Xephyr/pacman/check.py
} || {
	python Xephyr/pacman/check.py
}

mkdir -p build
cd build

#cmake options 
if [ "$1" = "midway" ] 
then
	cmake -DCMAKE_C_COMPILER='/usr/bin/gcc' -DCMAKE_CXX_COMPILER='/usr/bin/g++' ..
	make
elif [ "$1" = "stop" ] 
then
	echo 'STOPPING, you need to run your build command yourself'
elif [ "$1" = "" ] 
then
	cmake ..
	make
fi
cd ..
