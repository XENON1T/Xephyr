# Xephyr
# Userfriendly documentation, tutorials etc. at: https://xenon1t.github.io/Xephyr/docs/
# Doxigen style code reference, classes and description at: https://xenon1t.github.io/Xephyr/class_reference/

# Requirements:
- g++ >= 4.8.1
- gcc >= 4.8.1
- cmake >= 3.0
- ROOT >= 6.0.0 but reccomended is >= 6.10.04 to get all latest jupyter notebook features.
- python 3  (used for the package manager) 

# Quick Start:
### To load the library simply do: 
<code> root -l loadXephyr.C </code>

# Setting Up on Midway2 
### I couldn't manage on Midway1, the following works only on Midway2

- Setting up:
```bash
	mkdir XEPHYR
	cd XEPHYR
	git clone https://github.com/XENON1T/Xephyr.git Xephyr
	module load cmake			#version 3.6.2  
	module unload ROOT			
	module load  ROOT/6.06.08		#This is the only version that worked for me because of linking with sys library
	export XEPHYR_DIR=Your_dir
#	module unload python                    #this is not neded for the moment
#	module load  python/3.5.2+gcc-4.8       #this is not neded for the moment
```

- Cloning your code to compile against the library, in this case is SR1 (here we are only pulling the necessary directories):
```bash
	mkdir SR1
	cd SR1
	git init
	git remote add -f origin https://github.com/XENON1T/SR1Results.git
	git config core.sparseCheckout true
	echo "StatisticalAnalyses/xephyr_sr1_likelihood/" >> .git/info/sparse-checkout
	echo "StatisticalAnalyses/inputs_for_likelihood/" >> .git/info/sparse-checkout
	git pull origin master 
```

- Compiling Xephyr
```bash
	cd ..
	cp Xephyr/pacman/CMakeLists_copyMe.txt CMakeLists.txt   #pacman has a bug for now
	mkdir build
	cd build
	cmake -DCMAKE_C_COMPILER="/usr/bin/gcc" -DCMAKE_CXX_COMPILER="/usr/bin/g++" ..
	make
#	source Xephyr/pacman/build.sh        #all the above will in the future be done with this command but give problems now.
```
- 
