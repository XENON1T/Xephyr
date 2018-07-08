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

Currently Xephyr supports only 2D inputs (it is possible tough to use Xephyr for 1D likelihood, check out how to do it in `this example`_).
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
.. _`this example`: https://github.com/XENON1T/Xephyr/tree/master/examples/likelihood1D



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

        **IMPORTANT NOTES:** the likelihood expects pdfComponents to return models normalized to number of expected events. The likelihood
        rescales internally the parameter of interest (mu) to correspond to the number of "observed" signal events, for example, best fit mu=1.5 means
        that you are "fitting" 1.5 signal events. 

Xephyr provides four types of likelihood classes: 
 
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

 - **CombinedProfileLikelihood**: this class contains a list of likelihoods that need to be combined. Just add likelihood istances to its list and
   use same common methods. The combination is done just by multiplying likelihood togheter, sharing the parameter of interest, it takes into account 
   the relative exposure using the total signal integral for each likelihood. Nuissance parameter by default are added being considered independent from 
   each other, if you want to correlate them then you need to use the **CombinedParameter** class.

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

Xephyr uses the profiled likelihood approach to get rid of nuissance parameter depencency and the likelihood ratio as test statistic.
The parameter of interest (mu) is bound to be positive and the test statistic can be difined as two-sided or one-sided. 
If you are unfamiliar with this we recommend to have a look for example at `K. Cranmer`_ review.

 .. FIXME: add formula

.. _`k. Cranmer`: https://arxiv.org/abs/1503.07622

Sensitivity and limits may be computed in two different way: using asymptotics formulae or Neymann construction (following Feldmann-Cousins approach). 
For the asymptotics approach one can use the class `AsymptoticExclusion`_, the syntax is simple (either for limit and sensitivity) and you can refer to the tutorials for this (it is basically two lines of code), while for the construction of the test statistic distribution via MC it is worth to outline the procedure here.

.. _`AsymptoticExclusion`: https://xenon1t.github.io/Xephyr/class_reference/classAsymptoticExclusion.html

**Test statistic distro with MC:** The goal is to compute the graph of the test statistic's 90% quantiles as a funtion of 
the parameter of interest, this can be done in a following a few steps.

 #. Choose a list of parameter of interest (POI) values. Given that Xephyr is re-scaling internally the POI this list can be the same for all masses, the
    recommended one is mu = [ 0.5, 1, 2, 2.5, 3, 4, 5, 10, 15 ].
 #. Generate for each of the hypotheses in the list (and for each mass) MC toys datasets with nuber of signal injected events according to the hypothesis.
    For a final result with reasonable precision one should target 10'000 datasets per hypothesis. The datasets can be generated using the `ToyGenerator`_ class,
    see tutorials.
 #. Each of the generated dataset must be fitted with a conditional fit mu fixed to the true generated hypothesis. Don't worry, there is a class also for this,
    `ToyFitterExclusion`_ that does it semi-automaically, you just need to feed it with a preloaded likelihood. This will loop over a set of datasets 
    and store their post-fit parameter values, test statistic and other in a output ROOT Tree.
 #. Once You've got the post fits test statistic for all the samples (in the form of a ROOT Tree), you can extract their distribution with one of 
    the **plotHelpers::makeQuantiles**.

.. _`ToyGenerator`: https://xenon1t.github.io/Xephyr/class_reference/classToyGenerator.html
.. _`ToyFitterExclusion`: https://xenon1t.github.io/Xephyr/class_reference/classToyFitterExclusion.html
.. _`SR1`: https://github.com/XENON1T/SR1Results/tree/master/StatisticalAnalyses/xephyr_sr1_likelihood


**Computing limits with MC:** to compute limits for a single dataset you can use again the class `ToyFitterExclusion`_ and its method **spitTheLimit**, that will take as input
the 90% quantile graph generated before. In case you want to compute sensitivity, you need compute limit for a set of null hypothesys datatets, then compute the median and
1 and 2 sigma quantiles of their distribution. For this, again, there is an helper in **plotHelpers**. Have a look at the :ref:`tutorials <tutorials>` for more info,
also you can have a look to a real life example in the  `SR1`_ repository.



