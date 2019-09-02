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

===========
Quick usage
===========

If you created a conda environment to use GPathFinder (as proposed in installation), first you need to activate it with:

::

  conda activate name_of_the_environment

or

::

  source activate name_of_the_environment

Running GAUDI jobs is quite easy with :mod:`gaudi.cli.gaudi_run`. Put in your terminal:

::

    gpath run /path/to/input_file.yaml

You will need at least three input files. A `.yaml` file with the configuration of the job and two `.mol2` files for the ligand and the receptor molecules. To learn how to create input files, go to :ref:`input`. You can also check the tutorials :ref:`tutorial-input` and :ref:`tutorial-mol2`.

After the job is completed, you can use our home-made scripts to analyze them. A complete description of the output files is provided in :ref:`output`, and some tutorials on how to perform different analysis are available in :ref:`tutorial-output`.

To understand better the complete process of a GPathFinder calculation, you have also available the tutorial :ref:`tutorial-first`.
