.. _basics:

Xephyr basic concepts
=====================

In this section we describe the core component of Xephyr, we wont go trough the implementation details which is 
described in the :ref:`Tutorials <tutorials>`, but we'll go over the concepts.

Xephyr core libraries can be divided into 4 different main areas, for each of them there is a c++ class or a group of classes that handle it:
 

 - :ref:`Input Model/Data Handling <input_handling>`
 - :ref:`Nuissance Parameter Handling <nuissance>`, including interpolation.
 - :ref:`Likelihood building <likelihood>`
 - :ref:`Limit setting <limit>`: with Neymann construction and asymptotics.
 

.. _input_handling:

Input handling 
--------------

Currently Xephyr supports only 2D inputs (it is possible tough to use Xephyr for 1D likelihood, check out how to do it in **this example** `FIXME`_).
The input handling is done via two classes: `dataHandler`_, to handle science data, toy MC and calibration and a class that handles
inputs from the models as a `pdfComponent`_.

The **dataHandler** is basically just a helper class to handle ROOT Trees and is used to "fill the likelihood" with data points.
You can create a dataHandler in this way:

.. code-block:: c++

        dataHandler data(TString name, TString fileName, TString dmTree);
                       // just a name  // ROOT file name and Tree name inside the file

Currently the Tree must contain at least two variables (braches), the name of the branches can be set but the default is "cs1", "cs2". 
For neymann construction, when one needs to loop over many toy datasets (where each of them is a separated Tree), one can use the dataHandler to loop 
over all Tree at once using the convenient method **dataHandler::setTreeIndex**.

The **pdfComponent** handels instead the loading of PDF models through histograms. Again, currently this supports only **TH2F** as input, so you can
do 1D and 2D discriminating space with it but not more. 
It is an object representing a normalized PDF, meaning that is normalized to the total expected number of events for that background or signal,
This component is designed to interplay with :ref:`nuissance parameter <nuissance>`, so that the scale (total number of events) and the shape  of the histogram can change 
as a function of the parameters. The pdfComponent implements linear interpolation between any number of Shape parameters. To load models in the 
form of histogram you have to follow a naming convention for the histogram, here an example:

.. code-block:: c++

         prefixName_parameterNameA_x.yz_parameterNameB_x.yz......_parameterNameZ_x.yz

where the "prefixName" is whatever name for your bkg model, "parameterName" is the name of one of your parameters, and "x.y.z" are numbers.

You can find a detailed **tutorial** on how to use these classes for the most common cases in  :ref:`Tutorials <tutorials>` and you can find the official class documentation in
`pdfComponent`_ and `dataHandler`_ 

.. _`pdfComponent`: https://xenon1t.github.io/Xephyr/class_reference/classpdfComponent.html
.. _`dataHandler`: https://xenon1t.github.io/Xephyr/class_reference/classdataHandler.html



.. _nuissance:

Nuissance Parameter Handling
-----------------------------

Xephyr comes with two types of nuissance parameter (which in short is a fancy name to say uncertainties), a scale uncertainty, **scaleSys**, that takes
care of the absolute rate variation for that component, and a shape uncertainty, **shapeSys**, for its shape variation. We use the "natural" unit scale 
for the uncertainties: mean zero and sigma one. Each parameter can have three different behaviours: can be just a constant **fixed** parameter, or 
can be a **free** parameter, which is able to vary withouth any contraint, or can be a **nuissance** parameter, which can vary but is constrained by a gaussian
constrain.

For a **shapeSys** to be usable we need to load the PDF variations as a function of the parameter value, this is done loading different histogram 
corresponding to specific parameters value, so one need to build a kind of grid in the parameter space, populated with histograms, which aftewards
is going to be interpolated.

You can find the class documentation for the paramter handling in `shapeSys`_ and `scaleSys`_ reference, while a detailed tutorial on how 
to use them can be found in :ref:`Tutorials <tutorials>`. 

.. _`shapeSys`: https://xenon1t.github.io/Xephyr/class_reference/classshapeSys.html
.. _`scaleSys`: https://xenon1t.github.io/Xephyr/class_reference/classscaleSys.html

.. _likelihood:

Building the Likelihood
-------------------------

bla bla bla


.. _limit:

Limit setting
--------------

bla bla bla


