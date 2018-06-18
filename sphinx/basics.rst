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
This component is designed to integrate with nuissance parameter, so that the Scale (total number of events) and the shape either shape or scale (see next section).
You can immagine it as an 

You can find a detailed tutorial on how to use these classes for the most common cases in `FIXME`_

**(link class)** `FIXME`_
**(linking tutorial needed)**


.. _nuissance:

Nuissance Parameter Handling
-----------------------------

bla bla bla

.. _likelihood:

Building the Likelihood
-------------------------

bla bla bla


.. _limit:

Limit setting
--------------

bla bla bla


