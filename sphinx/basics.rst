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


        **Advanced:**
        Parameter that are used in combination of two likelihoods **combinedParameter** class. Exaplain a bit `FIXME`_


.. _likelihood:

Building the Likelihood
-------------------------

The likelihood classes are intended as a scaffold that collects your inputs and has the ability to run methods on them. The general idea is that a likelihood class
would be initialized with your inputs and  later on it can return the likelihood value for any parameter set. Note that the likelihood structure is predefined and
cannot be changed. There are a few types of likelihood classes and they differ based on the implementation of the likelihood function (see below),
but they all share common infrastructure like: printing parameter values and parameter access, a method called **maximize** that will actually minimize the -log(L),
methods computing likelihood scans, etc., you can find details on all of this in the :ref:`Tutorials <tutorials>` section.

Xephyr provides three types of likelihood classes: 
 
 - **pdfLikelihood**: implements the standard 2D (S1,S2) unbinned likelihood, it is composed of a poisson term that accounts for the total number of events times an 
   extended term that accounts for the shape, additionally each (non-free) nuissance parameter is constrained by a gaussian term. 
   This likelihood can be used also for a 1D PDF.

 .. FIXME formula

 - **binnedLikelihood**: implements a binned likelihood that can be either 1D or 2D, it returns just a product of poissons computed for each bin. TH2F errors
   on the bins are assigned as scale uncertainty to the bin itself. 
   Nowadays the trend is to go for a PDF likelihood because is "more sensitive", but remember this comes with high costs, since is harder to introduce uncertainties in a 
   PDF shape rather than consider them as scale independent factors on bins. So this type of likelihood is handy for all those search with low rate and large unknow
   uncertainties.

 .. FIXME formula
 
 - **myLikelihood**: this is just an empty scheleton, it provides all the functionalities but does not implement the method **computeTheLikelihood**, 
   so the developer is free to invent his own. Note that since the input infrastructure accepts 2D only histograms this is not an easy proxy for adding
   a dimension.

 .. FIXME formula


.. _safeguard:

The Safeguard
--------------

The idea of the safeguard is to try to address the difficulty of assigning 
uncertainty to the background model, especially when using an unbinned likelihood. This procedure needs a calibration dataset, the **assumption** 
is that the calibration represents exactly the background component of science data that one is trying to model. The science data and the calibration data
are then simultaneously fitted, so one need to provide the calibration data to the likelihood as well, a signal component (the safeguard) is injected then 
in shape of the background model and constrained by the fit on calibration. To turn On/Off the safeguard computation is quite easy, 
use the **setWithSafeGuard** flag. Currentl Xephyr implement only a positive injection of signal (which for a few reason can be considered more conservative), 
the possibility for a negative safeguard is under scrutiny and may be added in future releases.

        **WARNING:::** The safeguard is not the holy grail and it comes with a few subtle problems. The main problem is the contamination of your calibration sample
        by other background components (different from the ones you want to model). Xephyr allows you to introduce an additional background component to be 
        used only in the calibration fit to mimic the expected "known" contamination level. In case there is reason to think that the calibration sample, within the signal region,
        might  suffer of uknown or largely uncertain contamination the safeguard may lead to unwanted feature, hiding a real existing signal (for positive safeguard) or  
        even creating a fake signal (for negative safeguard).
        


.. _limit:

Limit setting
--------------

Xephyr uses the profiled likelihood approach



