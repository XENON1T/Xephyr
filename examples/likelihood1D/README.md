# 1D likelihood Example

Xephyr currently supports only a two dimensional discrimination space (for example S1 VS S2). This means that only a 2D likelihood is implemented.
We are working on implementing an N-dimensional likelihood, but for now it is not possible to use more that 2 dimensions.

However it is always possible to reduce the dimensionality from two to one (the purists wont like it, but hey, we are physicist after all).
This is an example to demonstrate how to generate some sample 1D data and 1D models  and feed them to a 2D likelihood.

### Setup 

The example comes as a Xephyr pachage (see the docs to know more about it), there are a couple of executable that need to be compiled.
To install the example do:

```bash
cd $XEPHYR_DIR
source Xephyr/pacman/installExample.sh likelihood1D # this just copies the content
source Xephyr/pacman/build.sh       # this compiles the executables of the example
``` 

This will create a directory "xephyr\_examples" in "XEPHYR\_DIR" and copy the code  content in it.  
 
### Generated the data and the models 

A script is made available to you to get inspiration, where an example of the production of mock data is shown. You can have a look at 
`produceDataAndModels.C`. To produce data simply do:

```bash
cd xephyr_examples/likelihood1D
root -l produceDataAndModels.C
``` 

This will create an additionl directory "data" and save the models in it.



### Likelihood, Scripts and Executables
| file name | Description |
|-----------|-------------|
| src/likelihoodDef.cxx | A simple likelihood setup is provided in this library |
| likelihoodShow.C | A simple example script to show a few helpfull feature, maximize the 1D-Likelihood and report the values of the parameters|
