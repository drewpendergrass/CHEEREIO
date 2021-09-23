.. _Deploying the Ensemble:

Deploying the ensemble
==========

Because CHEEREIO requires many of the ensemble settings to be globally available to many different types of scripts and programs, written in different languages and stored in in different directories, it expects a very particular kind of installation process and for user settings to be made available in a strict format. This is not too different from the way GEOS-Chem run directories are created as of version 13.0.0, but it does require some care on the part of the user. This page explains how to install the ensemble and customize it to meet your needs. Important science goes into this step!

.. _Configuration:

The Ensemble configuration file
-------------

If

**Important:** CHEEREIO constantly references the ``ens_config.json`` file throughout runtime. Changes to this file during runtime will change the settings of the ongoing assimilation calculation. However, there are some settings that are applied at the time of template run directory creation (e.g. cluster settings like memory usage) and cannot be adjusted at runtime via the ``ens_config.json`` file. **To be safe, you should redeploy the ensemble when you make changes to the configuration file and avoid making any changes during runtime.**

The Setup Ensemble script
-------------

TKTKTK

