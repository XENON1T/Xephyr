.. _quickstart:

Quick Start 
===========
This section is for all of those who **"I've got no time, just want to run the stupid thing..."**

There are 3 ways to use XEPHYR library:
 - Compiling the code as library and link to your own code.
 - load dynamically xephyr as root library, when using interactive root client.
 - Using xephyr in a jupyter notebook (only for brave hearts).


System Requirements
-------------------

Before you even start, you need to make sure you've got the following:


 - g++ >= 4.8.1
 - gcc >= 4.8.1
 - cmake >= 3.9
 - ROOT >= 6.0.0 but reccomended is >= 6.10.04 to get all latest jupyter notebook features.
 - python >=3 (used by the package manager)

**:::Warning:::** don't use ROOT installed from **Anaconda**, high are the chances that wont work.

**:::Warning2:::** make sure that your root installation has minuit, you can check it simply by running the following command in a ROOT CLI:

.. code-block:: bash

	$ root -l
	> ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad") == 0 ? "ERROR" : "You've got minuit!" 

 
 
Initial Setup
---------------

First create a Master directory (in this case called "XEPHYR_PKG") where to store Xephyr and your code, then clone Xephyr form github.

.. code-block:: bash

        mkdir XEPHYR_PKG
        cd XEPHYR_PKG
        export XEPHYR_DIR=$(pwd)
        echo "export XEPHYR_DIR=$(pwd)" >> ~/.bashrc
        git clone git@github.com:XENON1T/Xephyr.git .

You maybe noticed above the "export" of enviroment variable **XEPHYR_DIR** that points to the master dir XEPHYR_PKG, this tells Xephyr where its code is located, it is very important and so I reccomend to add it to your .bashrc (as the snipplet above is doing).

Compile Your Code With Xephyr
-----------------------------

XEPHYR strongly encourage package structure, so that you can easily share your code and compile it (hopefully) out of the box.
So, we strongly recommend you to :ref:`create your own package <packages>`, but since this is a "quick start for people in a hurry"  we'll skip that for now, 
you can always do it later. 

Xephyr provides a simple build command to compile any xephyr package automatically, **without writing any makefile** yourself. 
First step: the package you want to compile must lie inside the Master xephyr dir (**$XEPHYR_DIR**), which as a reminder is the dir that contains the Xephyr package.  Then run the following from the master dir:
 
.. code-block:: bash

       source Xephyr/pacman/build.sh

That's it! (if it worked) this should have created a new directory "build" where you can find all the executables of all your packages.
This step scans the directory tree and produces a makefile, so it is only needed at the beginning and when you add a new file, 
for all other small changes you can do:

.. code-block:: bash

       $ cd build
       $ make
        


Usage with Scripts 
--------------------

Xephyr provides a simple loader that can be used within a ROOT CLI session (or inside a macro) to load the library dinamically.
First, keep a structure (don't write your |pig| |pig|  scripts inside the Xephyr code dir), create a new directory where to store your scripts  
(inside the Master dir "XEPHYR_PKG") and copy the xephyr loader in it. 
Consider later to make a package out of it.

.. |pig| image:: static/pig.svg
        :width: 14pt       

.. code-block:: bash

        mkdir myAmazingCode
        cd myAmazingCode
        cp ../Xephyr/loadXephyr.C .

Now you can load xephyr libraries for **interactive use**, meaning you can run your scripts that use Xephyr classes out of the box (without any include statement).
To load the library for interactive use do:

.. code-block:: javascript 

        $ root -l loadXephyr.C
        // now I can use Xephyr classes in ROOT CLI :)
        > pdfComponent p("pippo","pluto.root")

Or if for some reason you prefer to call it within a script, then add to your script the following:

.. code-block:: c++

        {       
                // This is the top of your ROOT script 
                gROOT->ProcessLine(".x loadXephyr.C");

                // Here your code //
                // here I can use all Xephyr classes for example
                // pdfLikelihood * p(....)
        }


Jupyter Notebooks
-----------------

You brave hero! I pray for your soul... `FIXME: add at least a small comment here`_




.. raw:: html

        <div>Icons made by <a href="https://www.flaticon.com/authors/freepik" title="Pig">Pig</a> from <a href="https://www.flaticon.com/"     title="Flaticon">www.flaticon.com</a> is licensed by <a href="http://creativecommons.org/licenses/by/3.0/"     title="Creative Commons BY 3.0" target="_blank">CC 3.0 BY</a></div>


