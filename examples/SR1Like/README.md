# The SR1 likelihood example

This example is supposed to be as close as possible to a template for SR1 like analysis, 
where you basically just want to swap the signal model with your own and run the rest as the same.

This is a lighter and a bit more self-consistent version of the real [SR1 repository](https://github.com/XENON1T/SR1Results/tree/master/StatisticalAnalyses/xephyr_sr1_likelihood).

### What it needs

SR1 data repository [that you find here](https://github.com/XENON1T/SR1Results/tree/master/StatisticalAnalyses/inputs_for_likelihood) 
is needed I'd suggest to pull only the directory you need for the repository by doing the following, **warning this downloads 2.8Gb**:

```bash
cd $XEPHYR_DIR
mkdir SR1
cd SR1/
git init
git config core.sparseCheckout true
echo "StatisticalAnalyses/inputs_for_likelihood/lax_1.5.1_egg/"  >> .git/info/sparse-checkout
git remote add -f origin git@github.com:XENON1T/SR1Results.git
git pull origin master
```

Then you need to "install" this example by doing:

```bash
cd $XEPHYR_DIR
source Xephyr/pacman/installExample.sh SR1Like
```
This will copy the code in your $XEPHYR\_DIR, so that is easy for you to make your changes and compile 
the code, without being afraid of compromising Xephyr.

### What it does

The Xephyr SR1 likelihood is broken into 3 sub-likelihood, corresponding to:

- Inner Egg volume, R <= 34.5904 cm, mass 646.157 kg. 
- The so-called U-volume, 34.59 < R <= 36.48 cm,  mass 360.478 kg. 
- Wall volume, 36.48 < R <= 41.2396 cm, mass 238.908 kg.

**Total mass 1245.54 kg.**

Each of these sub-volume can be loaded separately, meaning that you can run inference using just one of them,
but the default configuration in this example is to combine the 3 volumes in a final likelihood.

### What shall I modify

There are basically just a couple of functions that you need to modify to run the SR1 analysis on your signal model 
and you can find it in the **"src/signalDef.cxx"**. Everything is explained in the code itself.


### More info

If you want to learn more how all this works I suggest you have a look at the [notebook tutorials](https://xenon1t.github.io/Xephyr/docs/tutorials.html).


