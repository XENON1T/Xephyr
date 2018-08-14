.. _tutorials:

Tutorials and Examples 
=======================

Notebook Tutorials
-------------------

You can find the HTML version of the following ROOT jupyter notebooks here:

- `tutorial likelihood`_: tutorial on how to setup a likelihood and run limits with asymptotics.
- `tutorial toy generation`_: First tutorial on Neymann construction, the generation of toys.
- `tutorial toy fitting and limits`_: Second tutorial on Neymann construction, generating the limits distributions.
- `simple script on limit plotting`_: Putting all togheter and plotting limits! (with nice colors).


.. _`tutorial likelihood`: https://github.com/XENON1T/Xephyr/blob/master/notebooks/asymptotic_inference.ipynb 
.. _`tutorial toy generation`: https://github.com/XENON1T/Xephyr/blob/master/notebooks/ToyGenerator_example.ipynb 
.. _`tutorial toy fitting and limits`: https://github.com/XENON1T/Xephyr/blob/master/notebooks/ToyFitter.ipynb
.. _`simple script on limit plotting`: https://github.com/XENON1T/Xephyr/blob/master/notebooks/Simple_limit_plotting_script.ipynb


Examples of Packages
---------------------

Togheter with tutorials, XEPHYR offers a list of examples to show you how actually implement things in real life.
The examples comes as XEPHYR packages and are very easy to install, once you have XEPHYR setup (see :ref:`quick start <quickstart>`) do:

.. code-block:: bash
        
        source Xephyr/pacman/installExample.sh name_of_example




**Xephyr supported Examples**

You find a list of valid examples in the xephyr `package itself`_, here a summary:

- `1D likelihood example:`_ how to "hack" Xephyr to feed 1D data and models to a 2D likelihood. 
- `SR1 like example:`_ this is a plug and play version of the SR1 analysis, you can easily plug in your signal model and 
  get a result using SR1 bkg models.


**Noticeable code examples (from other authors):**

- `OFFICIAL SR1 Xephyr likelihood package`_
- `Inelastic Dark Matter`_

.. _`OFFICIAL SR1 Xephyr likelihood package`: https://github.com/XENON1T/SR1Results/tree/master/StatisticalAnalyses/xephyr_sr1_likelihood
.. _`1D likelihood example:`: https://github.com/XENON1T/Xephyr/tree/master/examples/likelihood1D
.. _`Inelastic Dark Matter`: https://github.com/XENON1T/iDM
.. _`SR1 like example:`: https://github.com/XENON1T/Xephyr/tree/master/examples/SR1Like
.. _`package itself`:   https://github.com/XENON1T/Xephyr/tree/master/examples

If you write a package for your analysis we'd like to post here a link to your code. 
More example are under construction and will be posted on this page. 
