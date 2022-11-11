.. _LETKF modules:

The CHEEREIO LETKF implementation
==========

The 4D localized ensemble transform Kalman filter (4D-LETKF) algorithm lies at the heart of CHEEREIO. Everything else merely supports the assimilation operations conducted by this workflow.

Overview of the 4D-LETKF implementation
-------------

LETKF assimilation is implemented in CHEEREIO using a structure of nested Python objects, designed primarily to ensure that new observation operators can immediately plug into CHEEREIO and work automatically, without requiring users to have deep knowledge of the CHEEREIO code structure. We use Python because of its familiarity to a broader user base and because the object-oriented structure of the language makes it ideally suited to the modular design of CHEEREIO. Here we sketch the objects involved in the CHEEREIO implementation of the LETKF.
State vectors x_i^b are handled by an array of GC_Translator objects, one per ensemble member, which wrap around the GEOS-Chem restart file and scaling factors. These objects translate between gridded NetCDF files used by GEOS-Chem and the one-dimensional vectors needed by the LETKF process. The Assimilator object contains an array of GC_Translator objects and uses them to form the background perturbation matrix X_i^b, a key input into the LETKF algorithm.
To facilitate the 4D-LETKF calculation outlined in section 2.2, the HIST_Translator object (one per ensemble member) wraps around GEOS-Chem species concentration output files along with other files necessary for calculating simulated observations, such as meteorological information. Because the 4D version of LETKF requires aligning observations with the model state in time, the HIST_Translator object accesses model data from throughout the assimilation window. HIST_Translator objects are stored in an array within a HIST_Ens object (one per Assimilator object), which is responsible for aligning observations with simulated observations obtained by applying an observation operator to model output. The HIST_Ens object handles observations by using objects inheriting from the Observation_Translator class. In object-oriented programming, inheritance can be thought of as a sophisticated form of templating. As shown in Figure 4, the Observation_Translator class is a template that contains instructions to the user on how to write two standardized methods to (1) read observations from file and process them into a Python dictionary formatted for CHEEREIO, and (2) generate simulated observations y_i^b from GEOS-Chem output. Users can easily write their own class inheriting from Observation_Translator for a specific use case (like a particular surface or satellite instrument) by implementing these two methods, optionally employing a provided observation toolkit. Any class written with this strict template will then plug in automatically to the rest of the CHEEREIO workflow, and can be activated like usual from the main configuration file. CHEEREIO also comes with some pre-written observation operators (such as one for the TROPOMI satellite). Many different observation operators can be used simultaneously, making it natural to perform multispecies data assimilation or assimilation using both surface and satellite data within the CHEEREIO framework. The HIST_Ens object processes this information into a simulated observation perturbation matrix Y_i^b, which it passes to the Assimilator object.


The LETKF classes
-------------

TKTKTK

.. _GC Translator:

The GC Translator class
~~~~~~~~~~~~~

TKTKTK

.. _HIST Translator:

The HIST Translator class
~~~~~~~~~~~~~

TKTKTK

.. _HIST Ensemble:

The HIST Ensemble class
~~~~~~~~~~~~~

TKTKTK

The Observation Translator class type
~~~~~~~~~~~~~

TKTKTK

.. _GT Container:

The GT Container class
~~~~~~~~~~~~~

TKTKTK

.. _Assimilator:

The Assimilator class
~~~~~~~~~~~~~

TKTKTK