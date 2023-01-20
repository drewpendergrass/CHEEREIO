.. _Observations:

Observations
==========

This page is under construction, check back later!


The observation operator toolkit
-------------

CHEEREIO handles observations by using objects inheriting from the Observation_Translator class, a low-level translator which loads observations from file and compares them to GEOS-Chem output. In object-oriented programming, inheritance can be thought of as a sophisticated form of templating. Indeed, the Observation_Translator class itself is mostly empty, and contains instructions to the user on how to write two standardized methods to (1) read observations from file and process them into a Python dictionary formatted for CHEEREIO, and (2) generate simulated observations :math:`y_i^b` from GEOS-Chem output. Users can easily write their own class inheriting from Observation_Translator for a specific use case (like a particular surface or satellite instrument) by implementing these two methods, optionally employing a provided observation toolkit. Any class written with this strict template will then plug in automatically to the rest of the CHEEREIO workflow and can be activated from the main configuration file. CHEEREIO also comes with some pre-written observation operators (such as for the TROPOMI and OMI satellite instruments). Many different observation operators can be used simultaneously, making it natural to perform multispecies data assimilation or assimilation using both surface and satellite data within the CHEEREIO framework. Again, because Observation_Translators handle the details of interpreting a specific observation type, the rest of CHEEREIO can remain ignorant of specifics and operate in a fully abstract environment that can be reused for all simulations.

.. _ObsData:

The ObsData class
-------------

TKTKTK

Existing observation toolkits
-------------

TKTKTK

.. _TROPOMI tools:

TROPOMI tools
~~~~~~~~~~~~~

TKTKTK

.. _OMI tools:

OMI tools
~~~~~~~~~~~~~

TKTKTK


Supplementing an existing observation type
~~~~~~~~~~~~~

TKTKTK

.. _New observation:

Workflow to add a new observation operator
-------------

TKTKTKT

(1) Create a class inheriting from Observation_Translator 
~~~~~~~~~~~~~

TKTKTK.

(2) Implement getObservations() function 
~~~~~~~~~~~~~

TKTKTK

(3) Implement gcCompare() function 
~~~~~~~~~~~~~

TKTKTK

(4) Update operators.json
~~~~~~~~~~~~~

TKTKTK

(5) Link observational files from ens_config.json
~~~~~~~~~~~~~

TKTKTK

(6) [optional] Add observation filters via an extension
~~~~~~~~~~~~~

TKTKTK

.. _New superobservation:

(7) [optional] Add a new super observation error function
~~~~~~~~~~~~~

TKTKTK

