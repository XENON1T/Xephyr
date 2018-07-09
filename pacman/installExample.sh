#!/bin/bash

#colors
RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
PRPL='\033[0;35m'

dir='Xephyr/examples/'$1

if [ -d "$dir" ]
then
	mkdir -p xephyr_examples
	cp -r "$dir" xephyr_examples/.
	mv xephyr_examples/$1/_pkg.json xephyr_examples/$1/info.json
	echo  -e "\n\nExample ${RED}${1} ${NC} installed succesfully.   Now run:"
	echo -e   "${GREEN} source Xephyr/pacman/build.sh ${NC}\n\n"
else
	echo -e "\n\nExample ${RED}${1} ${NC} does not exist\n\n"
fi

