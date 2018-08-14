# Xephyr
# Userfriendly documentation, tutorials etc. at: https://xenon1t.github.io/Xephyr/docs/
# Doxigen style code reference, classes and description at: https://xenon1t.github.io/Xephyr/class_reference/

# Requirements:
- g++ >= 4.8.1
- gcc >= 4.8.1
- cmake >= 3.9
- ROOT >= 6.0.0 but reccomended is >= 6.10.04 to get all latest jupyter notebook features.
- python >=3  (used for the package manager) 
- if you use Anaconda ROOT, you will pay for your sins :fire:.


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
	export XEPHYR_DIR=Your_dir/    # the final slash is important
	module unload python                  
	module unload  mkl
	module load  python/3.5.2
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
	source Xephyr/pacman/build.sh

	# OR if you are on midway
	source Xephyr/pacman/build.sh midway
	
	# OR in case you need your custom cmake command and you just want the makefile to be build do:
	source Xephyr/pacman/build.sh stop
```
