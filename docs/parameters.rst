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

.. _parameters:

==================
List of parameters
==================


Parameters for version 1.2.0 
============================

**Genetic Algorithm parameters**

- **population** (int): size of the starting population, in number of individuals. **Default: 12**
- **generations** (int): number of generations to simulate. **Default: 500**
- **mu** (float): the number of children to select at each generation, expressed as a multiplier of ``ga.population``. **Default: 1**
- **lambda_** (float): the number of children to produce at each generation, expressed as a multiplier of ``ga.population``. **Default: 1.0**
- **mut_eta** (float): crowding degree of the mutation. A high eta will produce a mutant resembling its parent, while a small eta will produce a solution much more different. **Default: 5**
- **mut_pb** (float): the probability that an offspring is produced by mutation. **Default: 0.8**
- **mut_indpb** (float): independent probability for each gene to be mutated. **Default: 1.0**
- **cx_eta** (float): crowding degree of the crossover. A high eta will produce children resembling to their parents, while a small eta will produce solutions much more different. **Default: 5**
- **cx_pb** (float): the probability that an offspring is produced by crossover. **Default: 0.2**

**Path gene**

- **radius_rotamers** (float): Maximum distance (in Angstroms) from any point of the ligand in every frame that is searched for possible rotamers of the protein side-chain. **Default: 3.0**
- **max_step_separation** (float): Maximum distance (in Angstroms) from one point of the ligand in the pathway to the next one. If not set by the user, GPathFinder calculates the value from the size of the ligand. **Default: None**
- **min_step_increment** (float): Minimum distance increment (in Angstroms) from the ligand's origin  that has to be the ligand in one frame of the pathway with respect of the ligand's distance from the origin of the previous frame of the pathway. If not set by the user, GPathFinder calculates the value as 2/5 by max_step_separation. **Default: None**
- **mut_pos_pb** (float): When a mutation occurs, this value is the probability of such mutation to be of the type `positions`, that is, the mutation changes the actual trajectory of the pathway. Warning: parameter for advanced users, usually the default value is correct for the vast majority of the systems. **Default: 0.10**

**Path_torsion gene**

- **flexibility** (int or float): Maximum number of degrees a bond can rotate. **Default: 360.0**
- **anchor** (str): Molecule/atom_serial_number of reference atom for torsions. If not set, the nearest atom to the molecule geometric center is selected. **Default: None**
- **rotatable_atom_types** (list of str): Which type of atom types (as in chimera.Atom.idatmType) should rotate. **Default: ('C3', 'N3', 'C2', 'N2', 'P')**
- **rotatable_atom_names** (list of str): Which type of atom names (as in chimera.Atom.name) should rotate. **Default: ()**
- **non_rotatable_bonds** (list of str): Which bonds (identified by [Molecule/SerialNumber1, Molecule/SerialNumber2]) are not allowed to rotate. **Default: ()**
- **non_rotatable_selection** (str): Which bonds (identified by a Chimera selection query) are not allowed to rotate. **Default: ()**

**Path_rotamers gene**

- **library** {'Dunbrack', 'Dynameomics'}: The rotamer library to use. **Default: Dunbrack**

**Path_normalmodes gene**

- **method** {'prody', 'gaussian', 'pca'}: `prody` to calculate normal modes using prody algorithms, `gaussian` to read normal modes from a gaussian output file and `pca` to perform a PCA analysis over a MD trajectory (.dcd file). **Default: prody**
- **modes** (list of int): Modes to be used to move the molecule. **Default: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]**
- **pca_atoms** {'calpha', 'backbone', 'all'}: which atoms of the receptor are selected in the PCA analysis. **Default: calpha**
- **group_by** {'residues', 'mass', 'calpha', ''}: method to group atoms when using a coarse grain model of the receptor. **Default: residues**
- **group_lambda** (int): Either number of residues per group **(default=15)**, or total mass per group **(default=100)**
- **n_samples** (int): number of conformations to generate. **Default: 100**
- **rmsd** (float): average RMSD, in Angstroms, that conformations will have with respect to the initial conformation. **Default: 2.0**
- **minimize** (bool): whether to minimize the resulting samples or not. **Default: False**
- **forcefields** (list of str):  Used when `minimize` is `True` to indicate which forcefields to use. It can be a .prmtop file containing the parametrization and topology of the Protein. **Default: ('amber99sbildn.xml',)**
- **minimization_tolerance** (float): used when `minimize` is `True`. Convergence criteria for energy minimization. In kJ/mol. **Default: 10.0**
- **minimization_iterations** (int): used when `minimize` is `True`. Max attempts to converge at minimization. **Default: 1000**

**Path_scoring objective**

- **radius** (float): maximum distance, in Angstroms, from any point of the ligand in every frame that is searched for possible interactions. **Default: 5.0**
- **method** {'sum', 'average', 'max'}: Method used to calculate the score (i.e. sum, average or maximum of the scores of all frames). **Default: average**
- **clash_threshold** (float): used when the scoring is `clashes`. Maximum overlap of van-der-Waals spheres. If the overlap is greater, it's considered a clash. **Default: 0.6**
- **bond_separation** (int): used when the scoring is `clashes`. Ignore clashes between atoms within n bonds. **Default: 4**
- **same_residue** (bool): used when the scoring is `clashes`. Include intra-molecular clashes. **Default: True**
- **smoothness_threshold** (float): used when method is `smoothness`. RMSD between ligands on two consecutive frames that is permitted considering a perfect score of smoothness. **Default: 0.0**
- **smina_scoring** (str): used when the scoring is `smina` to specify alternative builtin scoring function (e.g. vinardo). **Default: None**
- **smina_custom_scoring** (str): used when the scoring is `smina` to specify a custom scoring function file. **Default: None**
- **smina_custom_atoms** (str): used when the scoring is `smina` to specify a custom atom type parameters file. **Default: None**

Parameters for version 1.1.0 
============================

**Genetic Algorithm parameters**

- **population** (int): size of the starting population, in number of individuals. **Default: 12**
- **generations** (int): number of generations to simulate. **Default: 500**
- **mu** (float): the number of children to select at each generation, expressed as a multiplier of ``ga.population``. **Default: 1**
- **lambda_** (float): the number of children to produce at each generation, expressed as a multiplier of ``ga.population``. **Default: 1.0**
- **mut_eta** (float): crowding degree of the mutation. A high eta will produce a mutant resembling its parent, while a small eta will produce a solution much more different. **Default: 5**
- **mut_pb** (float): the probability that an offspring is produced by mutation. **Default: 0.8**
- **mut_indpb** (float): independent probability for each gene to be mutated. **Default: 1.0**
- **cx_eta** (float): crowding degree of the crossover. A high eta will produce children resembling to their parents, while a small eta will produce solutions much more different. **Default: 5**
- **cx_pb** (float): the probability that an offspring is produced by crossover. **Default: 0.2**

**Path gene**

- **radius_rotamers** (float): Maximum distance (in Angstroms) from any point of the ligand in every frame that is searched for possible rotamers of the protein side-chain. **Default: 3.0**
- **max_step_separation** (float): Maximum distance (in Angstroms) from one point of the ligand in the pathway to the next one. If not set by the user, GPathFinder calculates the value from the size of the ligand. **Default: None**
- **min_step_increment** (float): Minimum distance increment (in Angstroms) from the ligand's origin  that has to be the ligand in one frame of the pathway with respect of the ligand's distance from the origin of the previous frame of the pathway. If not set by the user, GPathFinder calculates the value as 2/5 by max_step_separation. **Default: None**
- **mut_pos_pb** (float): When a mutation occurs, this value is the probability of such mutation to be of the type `positions`, that is, the mutation changes the actual trajectory of the pathway. Warning: parameter for advanced users, usually the default value is correct for the vast majority of the systems. **Default: 0.10**

**Path_torsion gene**

- **flexibility** (int or float): Maximum number of degrees a bond can rotate. **Default: 360.0**
- **anchor** (str): Molecule/atom_serial_number of reference atom for torsions. If not set, the nearest atom to the molecule geometric center is selected. **Default: None**
- **rotatable_atom_types** (list of str): Which type of atom types (as in chimera.Atom.idatmType) should rotate. **Default: ('C3', 'N3', 'C2', 'N2', 'P')**
- **rotatable_atom_names** (list of str): Which type of atom names (as in chimera.Atom.name) should rotate. **Default: ()**
- **non_rotatable_bonds** (list of str): Which bonds (identified by [Molecule/SerialNumber1, Molecule/SerialNumber2]) are not allowed to rotate. **Default: ()**
- **non_rotatable_selection** (str): Which bonds (identified by a Chimera selection query) are not allowed to rotate. **Default: ()**

**Path_rotamers gene**

- **library** {'Dunbrack', 'Dynameomics'}: The rotamer library to use. **Default: Dunbrack**

**Path_normalmodes gene**

- **method** {'prody', 'gaussian'}: `prody` to calculate normal modes using prody algorithms and `gaussian` to read normal modes from a gaussian output file. **Default: prody**
- **modes** (list of int): Modes to be used to move the molecule. **Default: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]**
- **group_by** {'residues', 'mass', 'calpha', ''}: method to group atoms when using a coarse grain model of the receptor. **Default: residues**
- **group_lambda** (int): Either number of residues per group **(default=15)**, or total mass per group **(default=100)**
- **n_samples** (int): number of conformations to generate. **Default: 100**
- **rmsd** (float): average RMSD, in Angstroms, that conformations will have with respect to the initial conformation. **Default: 2.0**
- **minimize** (bool): whether to minimize the resulting samples or not. **Default: False**
- **forcefields** (list of str):  Used when `minimize` is `True` to indicate which forcefields to use. It can be a .prmtop file containing the parametrization and topology of the Protein. **Default: ('amber99sbildn.xml',)**
- **minimization_tolerance** (float): used when `minimize` is `True`. Convergence criteria for energy minimization. In kJ/mol. **Default: 10.0**
- **minimization_iterations** (int): used when `minimize` is `True`. Max attempts to converge at minimization. **Default: 1000**

**Path_scoring objective**

- **radius** (float): maximum distance, in Angstroms, from any point of the ligand in every frame that is searched for possible interactions. **Default: 5.0**
- **method** {'sum', 'average', 'max'}: Method used to calculate the score (i.e. sum, average or maximum of the scores of all frames). **Default: average**
- **clash_threshold** (float): used when the scoring is `clashes`. Maximum overlap of van-der-Waals spheres. If the overlap is greater, it's considered a clash. **Default: 0.6**
- **bond_separation** (int): used when the scoring is `clashes`. Ignore clashes between atoms within n bonds. **Default: 4**
- **same_residue** (bool): used when the scoring is `clashes`. Include intra-molecular clashes. **Default: True**
- **smoothness_threshold** (float): used when method is `smoothness`. RMSD between ligands on two consecutive frames that is permitted considering a perfect score of smoothness. **Default: 0.0**
- **smina_scoring** (str): used when the scoring is `smina` to specify alternative builtin scoring function (e.g. vinardo). **Default: None**
- **smina_custom_scoring** (str): used when the scoring is `smina` to specify a custom scoring function file. **Default: None**
- **smina_custom_atoms** (str): used when the scoring is `smina` to specify a custom atom type parameters file. **Default: None**

Parameters for versions 1.0.x 
=============================

**Genetic Algorithm parameters**

- **population** (int): size of the starting population, in number of individuals. **Default: 12**
- **generations** (int): number of generations to simulate. **Default: 500**
- **mu** (float): the number of children to select at each generation, expressed as a multiplier of ``ga.population``. **Default: 1**
- **lambda_** (float): the number of children to produce at each generation, expressed as a multiplier of ``ga.population``. **Default: 1.0**
- **mut_eta** (float): crowding degree of the mutation. A high eta will produce a mutant resembling its parent, while a small eta will produce a solution much more different. **Default: 5**
- **mut_pb** (float): the probability that an offspring is produced by mutation. **Default: 0.8**
- **mut_indpb** (float): independent probability for each gene to be mutated. **Default: 1.0**
- **cx_eta** (float): crowding degree of the crossover. A high eta will produce children resembling to their parents, while a small eta will produce solutions much more different. **Default: 5**
- **cx_pb** (float): the probability that an offspring is produced by crossover. **Default: 0.2**

**Path gene**

- **radius_rotamers** (float): Maximum distance (in Angstroms) from any point of the ligand in every frame that is searched for possible rotamers of the protein side-chain. **Default: 3.0**
- **max_step_separation** (float): Maximum distance (in Angstroms) from one point of the ligand in the pathway to the next one. If not set by the user, GPathFinder calculates the value from the size of the ligand. **Default: None**
- **min_step_increment** (float): Minimum distance increment (in Angstroms) from the ligand's origin  that has to be the ligand in one frame of the pathway with respect of the ligand's distance from the origin of the previous frame of the pathway. If not set by the user, GPathFinder calculates the value as 2/5 by max_step_separation. **Default: None**
- **mut_pos_pb** (float): When a mutation occurs, this value is the probability of such mutation to be of the type `positions`, that is, the mutation changes the actual trajectory of the pathway. Warning: parameter for advanced users, usually the default value is correct for the vast majority of the systems. **Default: 0.10**

**Path_torsion gene**

- **flexibility** (int or float): Maximum number of degrees a bond can rotate. **Default: 360.0**
- **anchor** (str): Molecule/atom_serial_number of reference atom for torsions. If not set, the nearest atom to the molecule geometric center is selected. **Default: None**
- **rotatable_atom_types** (list of str): Which type of atom types (as in chimera.Atom.idatmType) should rotate. **Default: ('C3', 'N3', 'C2', 'N2', 'P')**
- **rotatable_atom_names** (list of str): Which type of atom names (as in chimera.Atom.name) should rotate. **Default: ()**
- **non_rotatable_bonds** (list of str): Which bonds (identified by [Molecule/SerialNumber1, Molecule/SerialNumber2]) are not allowed to rotate. **Default: ()**
- **non_rotatable_selection** (str): Which bonds (identified by a Chimera selection query) are not allowed to rotate. **Default: ()**

**Path_rotamers gene**

- **library** {'Dunbrack', 'Dynameomics'}: The rotamer library to use. **Default: Dunbrack**

**Path_normalmodes gene**

- **method** {'prody', 'gaussian'}: `prody` to calculate normal modes using prody algorithms and `gaussian` to read normal modes from a gaussian output file. **Default: prody**
- **modes** (list of int): Modes to be used to move the molecule. **Default: [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]**
- **group_by** {'residues', 'mass', 'calpha', ''}: method to group atoms when using a coarse grain model of the receptor. **Default: residues**
- **group_lambda** (int): Either number of residues per group **(default=15)**, or total mass per group **(default=100)**
- **n_samples** (int): number of conformations to generate. **Default: 100**
- **rmsd** (float): average RMSD, in Angstroms, that conformations will have with respect to the initial conformation. **Default: 2.0**
- **minimize** (bool): whether to minimize the resulting samples or not. **Default: False**
- **forcefields** (list of str):  Used when `minimize` is `True` to indicate which forcefields to use. It can be a .prmtop file containing the parametrization and topology of the Protein. **Default: ('amber99sbildn.xml',)**
- **minimization_tolerance** (float): used when `minimize` is `True`. Convergence criteria for energy minimization. In kJ/mol. **Default: 10.0**
- **minimization_iterations** (int): used when `minimize` is `True`. Max attempts to converge at minimization. **Default: 1000**

**Path_scoring objective**

- **radius** (float): maximum distance, in Angstroms, from any point of the ligand in every frame that is searched for possible interactions. **Default: 5.0**
- **method** {'sum', 'average', 'max'}: Method used to calculate the score (i.e. sum, average or maximum of the scores of all frames). **Default: average**
- **clash_threshold** (float): used when the scoring is `clashes`. Maximum overlap of van-der-Waals spheres. If the overlap is greater, it's considered a clash. **Default: 0.6**
- **bond_separation** (int): used when the scoring is `clashes`. Ignore clashes between atoms within n bonds. **Default: 4**
- **same_residue** (bool): used when the scoring is `clashes`. Include intra-molecular clashes. **Default: True**
- **smoothness_threshold** (float): used when method is `smoothness`. RMSD between ligands on two consecutive frames that is permitted considering a perfect score of smoothness. **Default: 0.0**
