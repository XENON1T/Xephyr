#!/bin/bash

#colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
PRPL='\033[0;35m'

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

build_success=123

#cmake options 
if [ "$1" = "midway" ] 
then
	cmake -DCMAKE_C_COMPILER='/usr/bin/gcc' -DCMAKE_CXX_COMPILER='/usr/bin/g++' ..
	make
	build_success=$?
elif [ "$1" = "stop" ] 
then
	echo 'STOPPING, you need to run your build command yourself'
	#	build_success=0
elif [ "$1" = "" ] 
then
	cmake ..
	make
	build_success=$?
fi
cd ..

echo
if [ $build_success -eq 0 ]
then
	echo
	echo
	echo -e "                  ${GREEN}SUCCESS!!!!! :) you are good to go!${NC}"
	echo -e "Congrats! Now all your compiled code, libraries and executables are in the ./build directory"
	xe_location=$(pwd)
	echo
	echo -e   "${GREEN}But remeber to set this environment variable (and add to your .bashrc), do: "
	echo -e   "--------------->     ${PRPL}export XEPHYR_DIR=$xe_location ${NC}"
	echo
	echo
fi
