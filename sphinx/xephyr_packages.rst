.. _packages:

Xephyr Packages
===============

Xephyr strongly encourages to keep your code organized and easy sharable via a simple package structure. In this section 
we describe the XEPHYR package manager **"pacman"**, which basically consist of a python script and some conventions.
The main purpose of this package manager is to make stupid simple to compile your private code and link it to the 
xephyr library, but also helps to keep your code organized.

Overview
--------

We use cmake to build the set of dependency for your code, so that you don't have to write a complicated makefile, everything is automatic.
But still, one need to tell to cmake where the code is and what you want to compile, this is what "pacman" is doing, so that following a 
few convention one can compile code written by others, without knowing anything about it.

The main idea is to split your code into 3 classes, let's fix the naming conventions: 

 - **Scripts**, these are all your ROOT scripts that will work only interactively with a ROOT CLI, they will not be compiled.
 - **Executable**, are all the files that must be compiled and for each one you expect an executable. These files contain "int main(){...}".
 - **Libraries**, most of the time you want to split your code into re-usable building blocks, for example classes, functions, so that you don't have to rewrite them 
   for each of your scripts or executables. At compiling time these libraries will be linked (togheter with XEPHYR libraries) to your executables. They can also be dinamically loaded from scripts (obviously).


Conventions
-----------

To build a XEPHYR package you must follow a number of conventions, mind that these are case sensitive:

  - First define a master dir where all the code is supposed to be. This must be refered with an env variable $XEPHYR_DIR. See the :ref:`Quick Start reference <quickstart>` for detail on how to do it.
  - Package definition: the package manager scans the whole master directory in search of a JSON file called **info.json**, each directory that contains that file is considered a XEPHYR package.
  - Within the package directory (but not in its subdirectories) all the files that end with **"_main.cxx"** are considered executables, and scheduled for being build.
  - Within each package it will look for a directory called **"src"**, where all the user library are. All the files in this directory that have a **".cxx"** extension are scheduled to be build as library.
  - You can put your scripts anywhere, as long as they are not in "src" nor their name end with "_main.cxx", the package manager will just ignore them.


Setup
-----

In the master dir make a new directory for your package and initialize it by creating the json file:

.. code-block:: bash
        
        cd $XEPHYR_DIR
        mkdir myAwsomePackage
        echo '{ "pkg_name": "package_name", 
		"pkg_version" : "0.0", 
		"dependencies": ["ROOT==6.1.0"] }'  > myAwsomePackage/info.json 

The "myAwsomePackage" directory can be also a subdirectory nested in whatever complicated directory tree and it will work anyway, 
the important thing is that it is inside $XEPHYR_DIR. Now create the content of your package:

.. code-block:: bash

	cd myAwsomePackage
	mkdir src

	echo  ' // Xephyr includes examples 
		#include "XePdfObjects.h"
		#include "XeLikelihoods.h"

	       //some fancy code here... 
	      '                               > src/fancyLib.cxx

	echo  ' // Include your library
		#include "fancyLib.cxx"
		
	       int main()
		{
		// your code here
		}
		'                              > myExe_main.cxx
 

**Did you notice that all the includes are global?**  You don't need to specify the directory of Xephyr files or library files you want to include!!!
Note that this feature is not available with scripts, so if you want to load your "fancyLib.cxx" from a script you would add to your script:

.. code-block:: c++

	// in case you want to load the library inside a script
	 gROOT->ProcessLine(".L PATH_TO_PAKAGE/src/fancyLib.cxx");
	
	//======================================================//
	
	// Or in case you want to use the include syntax in a macro you must run the 
        // following from the CLI (or from a sript)
	> gInterpreter->AddIncludePath("PATH_TO_PAKAGE/src")
	// Now you can use global includes even iside your macro.
	



Build
------

Now that you created your pakage is time to compile it and link it. Go back to the master directory and use the pakage manager script:

.. code-block:: bash

        cd $XEPHYR_DIR
        source Xephyr/pacman/build.sh


That's it! (if it worked) this should have created a new directory "build" where you can find all the executables of all your packages.
This step scans the directory tree and produces a makefile, so it is only needed at the beginning and when you add a new file, 
for all other small changes you can do:

.. code-block:: bash

       $ cd build
       $ make




        
