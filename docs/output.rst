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

.. _output:

============
Output files
============

This section aims to provide the user with a complete description of all the files returned as output of a GPathFinder calculation. To learn how to perform different analysis over this data, you can check the tutorial :ref:`tutorial-output`.

The results of the calculation will be saved inside a folder on the path indicated in the ``output.path`` parameter of the `.yaml` input file. This folder will contain:

- **A `.gaudi-log` file**: contains the log information of the calculation, with data about the time employed in the execution, the evolution of the scores along the calculation and possible errors/warnings.
- **A `.yaml` file**: contains a replica of the `.yaml` input file employed in the calculation.
- **A `.gpath-output` file**: contains a summary of the obtained results, with data about their scores and name of the file that actually has the pathway. Can be used to get an overview of the quality of the solutions.
- **A `summary.csv` file**: contains a summary of the scoring of all obtained results, detailed by frame. Also contains the coordinates of every solution (considering the center of mass of the ligand at every frame of the pathway). Can be used to see at a frame detail the quality of a concrete solution and differences between solutions.
- **A set of `.zip` files**: contain the different pathways obtained from the calculation. The number of solutions will be equal to the population size if ``output.pareto`` was set to `False` or equal to the size of the pareto frontier (i.e. dominant solutions) otherwise.

**Optional files**

- **A `.nmd` file**: contains the information about the prody modes calculated in Normal Mode Analysis of the receptor molecule. This file only will be present if ``gaudi.genes.path_normalmodes.write_modes`` was set to `True`.
- **A `.samples` file**: contains the information about all the samples generated during Normal Mode Analysis of the receptor molecule. This file only will be present if ``gaudi.genes.path_normalmodes.write_samples`` was set to `True`.

Contents of the `.zip` corresponding to each solution
=====================================================

Each `.zip` file contains all the necessary information to reproduce at atomic level the pathway proposed by GPathFinder as a solution for your problem. Inside the file there are the following contents:

- **Two `.mol2` files** with the original 3D structures of the ligand and the receptor.
- **A `.gaudi` file** with the summary of the files.
- **Another `.zip` file** which contains a set of `.pdb` files with all the frames of the pathway (each one has a conformation of the ligand and receptor molecules). You can open these `.pdb` files in any visualization tool like a MD-movie and work with them and analyze each frame. An example of a pathway that represents a GPathFinder solution can be found `here <https://raw.githubusercontent.com/insilichem/gpathfinder/master/examples/output_files/example_pathway.zip>`_. 

Moreover, an `allele.txt` and a `scores.txt` files are present with all the necessary information to reconstruct the pathway from the original structures and the score information for every frame. These additional files can be used by your own scripts to analyze the results in further detail.

Finally, a `trajectory.pdb` file which contains the route that follows the ligand, as a set of points that represent the center of the ligand at each step of the trajectory. This can be used to easily compare the trajectory of different solutions obtained from GPathFinder calculations.
