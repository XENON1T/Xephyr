# 1D likelihood Example

Xephyr currently supports only a two dimensional discrimination space (for example S1 VS S2). This means that only a 2D likelihood is implemented.
We are working on implementing an N-dimensional likelihood, but for now it is not possible to use more that 2 dimensions.

However it is always possible to reduce the dimensionality from two to one (the purists wont like it, but hey, we are physicist after all).
This is an example to demonstrate how to generate some sample 1D data and 1D models  and feed them to a 2D likelihood.

### Setup 

The example comes as a Xephyr pachage (see the docs to know more about it), there are a couple of executable that need to be compiled.
To install the example do:

This will copy the content of the example in the "XEPHYR\_DIR". 
 
### Generated the data and the models 

A script is made available to you to get inspiration, where an example of the production of mock data is shown. You can have a look at 
`produceDataAndModels.C`. To produce data simply do:

```bash
root -l produceDataAndModels.C
``` 


 
