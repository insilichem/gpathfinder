.. GPathFinder: Identification of ligand binding pathways 
.. by a multi-objective genetic algorithm

   https://github.com/insilichem/gpathfinder

   Copyright 2019 José-Emilio Sánchez Aparicio, Giuseppe Sciortino,
   Daniel Villadrich Herrmannsdoerfer, Pablo Orenes Chueca, 
   Jaime Rodríguez-Guerra Pedregal and Jean-Didier Maréchal
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

==============
How to install
==============

Recommended steps:

1 - Download the `latest stable copy of UCSF Chimera <http://www.cgl.ucsf.edu/chimera/download.html>`_ and install it with the `.dmg` file (macOS) or with the following command (Linux):

::

  chmod +x chimera-*.bin && sudo ./chimera-*.bin
  
.. tip:: 

   When Chimera installer asks about `Install symbolic link to chimera executable?`, we recommend to choose option `/usr/bin`, so GaudiMM will find Chimera installation without problem.

2 - Download `Miniconda Python 2.7 Distribution <http://conda.pydata.org/miniconda.html>`_ for your platform and install it with:

::

  bash Miniconda2*.sh

3 - Install ``gpathfinder`` with ``conda`` in a new environment called ``gpathfinder`` (or whatever name you prefer after the ``-n`` flag), using these custom channels (``-c`` flags):

::

  conda create -n gpathfinder -c omnia -c conda-forge -c bioconda -c josan_bcn gpathfinder


4 - Activate the new environment as proposed:

::

  conda activate gpathfinder

or

::

  source activate gpathfinder
 

5 - Run it!

::

  gpath


Check everything is OK
======================

If everything went OK, you will get the usage screen:

.. code-block:: console

    Usage: gpath [OPTIONS] COMMAND [ARGS]...

      GPathFinder: Indentification of ligand pathways by a multi-objective
      genetic algorithm

      (C) 2019, InsiliChem
      https://github.com/insilichem/gpathfinder

    Options:
      --version   Show the version and exit.
      -h, --help  Show this message and exit.

    Commands:
      prepare  Create or edit a GPATH input file.
      run      Launch a GPATH input file.
      view     Analyze the results in a GPATH output file.

OS Compatibility
================

GPathFinder is compatible with Linux and macOS. This installation procedure has been checked with Chimera v.1.13.1 and the following distributions:

- macOS Mojave 10.14

- Mint 19.1
- Debian 9.9.0
- Ubuntu 16.04 and 18.04
- OpenSUSE Leap 15.1
- Manjaro 18.0.4

If you find some difficulties when installing it in a concrete distribution, please use the `issues page <https://github.com/insilichem/gpathfinder/issues>`_ to report them.
