.. _Deploying the Ensemble:

Deploying the ensemble
==========

Because CHEEREIO requires many of the ensemble settings to be globally available to many different types of scripts and programs, written in different languages and stored in in different directories, it expects a very particular kind of installation process and for user settings to be made available in a strict format (so that they can be read by scripts in multiple languages). The overall install process is not too different from the way GEOS-Chem run directories are created as of version 13.0.0, but it does require some care on the part of the user. This page explains how to install the ensemble and customize it to meet your needs. Important science goes into this step!

.. _Configuration:

The Ensemble configuration file
-------------

All ensemble configuration is set by the ``ens_config.json`` file in the main CHEEREIO code directory. This file contains all the information that the user would normally input in the process of making a GEOS-Chem run directory, in addition to other settings like cluster configuration, assimilation global variables (like the localization radius), the composition of the state and control vectors, links to observation operators, and details about which emissions can be assimilated. This file is in JSON format, which is a file type that is human readable but expects a fairly strict syntax. 

**Important:** CHEEREIO constantly references the ``ens_config.json`` file throughout runtime. Changes to this file during runtime will change the settings of the ongoing assimilation calculation. However, there are some settings that are applied at the time of template run directory creation (e.g. cluster settings like memory usage) and cannot be adjusted at runtime via the ``ens_config.json`` file. **To be safe, you should redeploy the ensemble when you make changes to the configuration file and avoid making any changes during runtime.**

Let's walk through the configuration file section-by-section, after first covering some basic formatting requirements.

JSON formatting
~~~~~~~~~~~~~

JSON stores data in a form very similar to a Python dictionary; more details on the format are available on the JSON `website <https://www.json.org>`__. 


Basic GEOS-Chem and ensemble settings
~~~~~~~~~~~~~


Cluster settings
~~~~~~~~~~~~~


Species in state/control/observation vectors
~~~~~~~~~~~~~


Miscellaneous LETKF settings
~~~~~~~~~~~~~


The Setup Ensemble script
-------------

TKTKTK

