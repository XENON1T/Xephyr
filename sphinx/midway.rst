Setting Up on Midway2 
=======================

I couldn't manage on Midway1, the following works only on Midway2

Setting up:
-----------

.. code-block:: bash

        mkdir XEPHYR
	cd XEPHYR
	git clone https://github.com/XENON1T/Xephyr.git Xephyr
	module load cmake			#version 3.6.2  
	module unload ROOT			
	module load  ROOT/6.06.08		#This is the only version that worked for me because of linking with sys library
	export XEPHYR_DIR=Your_Work_dir/    # the final slash is important
	module unload python                  
	module unload  mkl
	module load  python/3.5.2



Compiling Xephyr
----------------

.. code-block:: bash

	source Xephyr/pacman/build.sh midway
	
	# OR in case you need your custom cmake command and you just want the makefile to be build do:
	source Xephyr/pacman/build.sh stop

