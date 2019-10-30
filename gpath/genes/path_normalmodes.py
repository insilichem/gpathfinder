#!/usr/bin/env python
# -*- coding: utf-8 -*-

##############
# GPathFinder: Identification of ligand pathways by a multi-objective
# genetic algorithm
# 
# https://github.com/insilichem/gpathfinder
#
# Copyright 2019 José-Emilio Sánchez Aparicio, Giuseppe Sciortino,
# Daniel Villadrich Herrmannsdoerfer, Pablo Orenes Chueca, 
# Jaime Rodríguez-Guerra Pedregal, Laura Tiessler-Sala and 
# Jean-Didier Maréchal
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############

"""
This module allows to explore molecular folding of the receptor
through normal modes analysis.
It works by calculating normal modes for the input molecule and moving 
along a combination of normal modes.
It needs at least a GPathFinder `path` gene and a GPathFinder `molecule` gene.
"""

# Python
from __future__ import print_function, division
import random
import logging
import numpy
import os
import pickle
# Chimera
import chimera
from _multiscale import get_atom_coordinates
# 3rd party
import simtk.openmm.app as openmm_app
import prody
from boltons.cacheutils import LRU
from cclib.parser import Gaussian
# GPATH
from gpath.genes import GeneProvider
from gpath import parse
from gpath.objectives import path_energy as energy

logger = logging.getLogger(__name__)
_openmm_builtin_forcefields = os.listdir(os.path.join(openmm_app.__path__[0], 'data'))

def enable(**kwargs):
    kwargs = NormalModes.validate(kwargs)
    return NormalModes(**kwargs)

class NormalModes(GeneProvider):

    """
    NormalModes class
    Parameters
    ----------
    method : str
        Either:
        - prody : calculate normal modes using prody algorithms
        - gaussian : read normal modes from a gaussian output file
        - pca: calculate normal modes using pca analysis over a trajectory
    target : str
        Name of the Gene containing the actual molecule
    modes : list, optional, default=range(20)
        Modes to be used to move the molecule
    pca_atoms : str, optional
        Either:
        - calpha : use alpha carbon to calculate normal modes
        - backbone : use all backbone atoms to calculate normal modes 
        - all : use all atoms of the receptor to calculate normal modes
        Default to calpha
    group_by : str or callable, optional, default=None
        group_by_*: algorithm name or callable
        coarseGrain(prm) which makes ``mol.select().setBetas(i)``,
        where ``i`` is the index Coarse Grain group,
        and ``prm`` is ``prody.AtomGroup``
    group_lambda : int, optional
        Either: number of residues per group (default=15), or
        total mass per group (default=100)
    path : str
        Gaussian or prody modes output path. Required if ``method`` is
        ``gaussian``.
    samples_path : str, optional
        Path of the samples information file. If not None, the gene
        will use the samples contained in the file instead of generating
        them from scratch.
    trajectory : str, optional
        Path to a .dcd file containing a trajectory for pca analysis.
        Default to None.
    write_modes: bool, optional
        write a ``modes.nmd`` file with the ProDy modes
    write_samples: bool, optional
        write a ``samples.samples`` file with the samples information
    n_samples : int, optional, default=100
        number of conformations to generate
    rmsd : float, optional, default=2.0
        average RMSD that the conformations will have with respect
        to the initial conformation
        forcefields : list of str, optional, default=('amber99sbildn.xml',)
        Used when method is ``energy`` to indicate which forcefields to 
        use.
        It can be a .prmtop file containing the parametrization and 
        topology of the whole system Protein + Ligand (in this order).
    minimize : bool, optional, defaults to False
        Whether to minimize the resulting samples or not.
    minimization_tolerance : float, optional, defaults to 10
        Used when ``minimize`` is True. Convergence criteria for energy 
        minimization. In kJ/mol.
    minimization_iterations : int, optional, defaults to 100
        Used when ``minimize`` is True. Max attempts to converge at 
        minimization.
    forcefields : list of str, optional, default=('amber99sbildn.xml',)
        Used when ``minimize`` is True to indicate which forcefields to 
        use.
        It can be a .prmtop file containing the parametrization and 
        topology of the Protein.
    auto_parametrize: list of str, optional, default=None
        Used when ``minimize`` is True. List of ``molecule`` GPathFinder gene 
        instances that will be tried to auto parametrize with antechamber.
    parameters : list of 2-item list of str, optional, default=None
        Used when ``minimize`` is True. List of (gaff.mol2, .frcmod) 
        files to use as parametrization source.
    frozen_atoms : str
        Used when ``minimize`` is True. MDTraj DSL query that will 
        select atoms to be froze; e.g., ``backbone``
    system_options : dict
        Used when ``minimize`` is True. Key-value options that will be 
        passed to ``.createSystem()`` calls. Useful for implicit solvent 
        declaration, custom non bonded methods, and so on.
    Attributes
    ----------
    allele : slice of prody.ensemble
        Randomly picked coordinates from NORMAL_MODE_SAMPLES
    NORMAL_MODES : prody.modes
        normal modes calculated for the molecule or readed
        from the gaussian frequencies output file stored
        in a prody modes class (ANM or RTB)
    NORMAL_MODE_SAMPLES : prody.ensemble
        configurations applying modes to molecule
    _original_coords : numpy.array
        Parent coordinates
    _chimera2prody : dict
        _chimera2prody[chimera_index] = prody_index
    """

    _validate = {
        'method': parse.In(['prody', 'gaussian','pca']),
        'path': parse.RelPathToInputFile(),
        'samples_path': parse.RelPathToInputFile(), 
        'trajectory': parse.RelPathToInputFile(),
        'write_modes': bool,
        'write_samples': bool,
        parse.Required('target'): parse.Molecule_name,
        'pca_atoms': parse.In(['calpha', 'backbone','all']),
        'group_by': parse.In(['residues', 'mass', 'calpha', '']),
        'group_lambda': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'modes': [parse.All(parse.Coerce(int), parse.Range(min=0))],
        'n_samples': parse.All(parse.Coerce(int), parse.Range(min=1)),
        'rmsd': parse.All(parse.Coerce(float), parse.Range(min=0)),
        'minimize': bool,
        'minimization_tolerance': float,
        'minimization_iterations': int,
        'forcefields': [parse.Any(parse.ExpandUserPathExists, parse.In(_openmm_builtin_forcefields))],
        'auto_parametrize': [parse.Molecule_name],
        'parameters': [parse.All([parse.ExpandUserPathExists], 
                       parse.Length(min=2, max=2))],
        'frozen_atoms': parse.Any(None, parse.Coerce(str)),
        'system_options': parse.Any(dict, None)
    }

    def __init__(self, method='prody', target=None, modes=None, trajectory=None, n_samples=100, rmsd=2.0,
                 pca_atoms ='calpha',group_by='residues', group_lambda=None, samples_path=None,
                 path=None, write_modes=False, write_samples=False,
                 minimize=False, minimization_tolerance=10, minimization_iterations=1000,
                 forcefields=('amber99sbildn.xml',), auto_parametrize=None,
                 parameters=None, system_options=None, frozen_atoms=None, **kwargs):
        # Fire up!
        GeneProvider.__init__(self, **kwargs)
        self.method = method
        self.target = target
        self.trajectory = trajectory
        self.modes = modes if modes is not None else range(20)
        self.max_modes = max(self.modes) + 1
        self.n_samples = n_samples
        self.rmsd = rmsd
        self.group_by = None
        self.group_by_options = None
        self.path = None
        self.samples_path = samples_path
        self.write_modes = write_modes
        self.write_samples = write_samples
        self._need_express = False #Control of expression by GPathFinder
        self.minimize = minimize
        self.forcefields = forcefields
        self.auto_parametrize = auto_parametrize
        self.parameters = parameters
        self.minimization_tolerance = minimization_tolerance
        self.minimization_iterations = minimization_iterations
        self.system_options = system_options
        self.frozen_atoms = frozen_atoms
        self._discarded_samples = []
        if method == 'prody':
            if path is None:
                self.normal_modes_function = self.calculate_prody_normal_modes
                self.group_by = group_by
                self.group_by_options = {} if group_lambda is None else {'n': group_lambda}
            else:
                self.path = path
                self.normal_modes_function = self.read_prody_normal_modes
        
        elif method == 'pca':
            if path is None:
                self.normal_modes_function = self.calculate_pca_normal_modes 
                self.pca_atoms = pca_atoms
            else:
                self.path = path
                self.normal_modes_function = self.read_prody_normal_modes

        else:  # gaussian
            self.normal_modes_function = self.read_gaussian_normal_modes
            if path is None:
                raise ValueError('Path is required if method == gaussian')
            self.path = path

        if self.name not in self._cache:
            self._cache[self.name] = LRU(300)

    def __ready__(self):
        """
        Second stage of initialization
        It saves the parent coordinates, calculates the normal modes and initializes the allele
        """
        cached = self._CACHE.get('normal_modes')
        if not cached:
            normal_modes, normal_modes_samples, chimera2prody, prody_molecule = self.normal_modes_function()
            self._CACHE['normal_modes'] =  normal_modes
            self._CACHE['chimera2prody'] =  chimera2prody
            self._CACHE['original_coords'] =  chimeracoords2numpy(self.molecule)
            self._discarded_samples = []

            #MINIMIZATION
            if self.minimize == True:
                print('Minimizing samples...')
                self.energy_minimizer = energy.Energy(targets=[self.target],
                                            forcefields=self.forcefields, 
                                            parameters=self.parameters,
                                            platform='CPU',
                                            minimization_tolerance=self.minimization_tolerance,
                                            minimization_iterations=self.minimization_iterations,
                                            system_options=self.system_options,
                                            frozen_atoms=self.frozen_atoms)
                for i, sample in enumerate(normal_modes_samples):
                    self.allele = sample
                    self._need_express = True
                    self.express()                 
                    positions, kj = self.energy_minimizer.minimize_energy([self.molecule], self.parent)
                    self.unexpress()
                    normal_modes_samples[i] = numpy.array([list(p) for p in positions])
                    if abs(kj) > 1000000: #discard sample if it has very bad energy
                    	self._discarded_samples.append(i)
                        print('Sample ', str(i) , '/', str(len(normal_modes_samples)-1), ' : ', str(kj), 
                            ' -> DISCARDED due to high energy')
                    else:
                        print('Sample ', str(i) , '/', str(len(normal_modes_samples)-1), ' : ', str(kj))
                    self._need_express = False
                self.energy_minimizer = None
            self._CACHE['normal_modes_samples'] = normal_modes_samples
            if self.write_modes:
                title = os.path.join(self.parent.cfg.output.path,
                                     'modes.nmd')
                prody.writeNMD(title, normal_modes, prody_molecule)
            #Added Jose
            if self.write_samples:
                file_name = os.path.join(self.parent.cfg.output.path,
                                        'samples.samples')
                with open(file_name, 'wb') as f:
                    pickle.dump(self.NORMAL_MODES_SAMPLES, f)
                    f.close()
    	self.allele = []
    	if len(self.NORMAL_MODES_SAMPLES) == len(self._discarded_samples):
    		raise ValueError('All samples have very bad energies after minimization. \
    				Check the structure of the receptor and the parameters of \
    				path_normalmodes gene.')
    	while len(self.allele) == 0:
        	sample = random.randint(0, len(self.NORMAL_MODES_SAMPLES)-1)
        	if sample not in self._discarded_samples:
        		self.allele = self.NORMAL_MODES_SAMPLES[sample]

    def express(self):
        """
        Apply new coords as provided by current normal mode
        """
        if self._need_express: #Control of expression by GPathFinder
            c2p = self._chimera2prody
            for atom in self.molecule.atoms:
                index = c2p[atom.serialNumber]
                new_coords = self.allele[index]
                atom.setCoord(chimera.Point(*new_coords))

    def unexpress(self):
        """
        Undo coordinates change
        """
        if self._need_express: #Control of expression by GPathFinder
            for i, atom in enumerate(self.molecule.atoms):
                atom.setCoord(chimera.Point(*self._original_coords[i]))

    def mate(self, mate):
        """
        .. todo::
            Combine coords between two samples in NORMAL_MODES_SAMPLES?
            Or two samples between diferent NORMAL_MODES_SAMPLES?
            Or combine samples between two NORMAL_MODES_SAMPLES?
            For now : pass
        """
        pass

    def mutate(self, indpb):
        """
        Not needed. Controlled in ``path`` gene mutation.
        """
        pass

    def write(self, path, name, *args, **kwargs):
        return None

    #####
    @property
    def molecule(self):
        return self.parent.genes[self.target].compound.mol

    @property
    def _CACHE(self):
        return self._cache[self.name]

    @property
    def NORMAL_MODES(self):
        return self._CACHE.get('normal_modes')

    @property
    def NORMAL_MODES_SAMPLES(self):
        return self._CACHE.get('normal_modes_samples')

    @property
    def _chimera2prody(self):
        return self._CACHE.get('chimera2prody')

    @property
    def _original_coords(self):
        return self._CACHE.get('original_coords')

    def calculate_prody_normal_modes(self):
        """
        calculate normal modes, creates a diccionary between chimera and prody indices
        and calculate n_confs number of configurations using this modes
        """
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = prody_modes(prody_molecule, self.max_modes, self.group_by,
                            **self.group_by_options)
        samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                    n_confs=self.n_samples, rmsd=self.rmsd)
        samples.addCoordset(prody_molecule)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody, prody_molecule

    def calculate_pca_normal_modes(self):
        """
        First it creates a diccionary between chimera and prody indices, then calculates the pca normal modes 
        and calculate n_confs number of configurations using this modes
        """

        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule) 

        modes = pca_modes(prody_molecule, self.trajectory, self.max_modes, self.pca_atoms) 
        samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule, 
                                   n_confs=self.n_samples, rmsd=self.rmsd) 
        samples.addCoordset(prody_molecule)
        samples_coords = [sample.getCoords() for sample in samples]
        return modes, samples_coords, chimera2prody, prody_molecule

    def read_prody_normal_modes(self):
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = prody.parseNMD(self.path)[0]
        if not self.samples_path: 
            samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                        n_confs=self.n_samples, rmsd=self.rmsd)
            samples.addCoordset(prody_molecule)
            samples_coords = [sample.getCoords() for sample in samples]
        else: #Read the samples from .samples file
            with open(self.samples_path, 'rb') as f:
                samples_coords = pickle.load(f)
                f.close()
        return modes, samples_coords, chimera2prody, prody_molecule

    def read_gaussian_normal_modes(self):
        """
        read normal modes, creates a diccionary between chimera and prody indices
        and calculate n_confs number of configurations using this modes
        """
        prody_molecule, chimera2prody = convert_chimera_molecule_to_prody(self.molecule)
        modes = gaussian_modes(self.path)
        if not self.samples_path: 
            samples = prody.sampleModes(modes=modes[self.modes], atoms=prody_molecule,
                                        n_confs=self.n_samples, rmsd=self.rmsd)
            samples.addCoordset(prody_molecule)
            samples_coords = [sample.getCoords() for sample in samples]
        else: #Read the samples from .samples file
            with open(self.samples_path, 'rb') as f:
                samples_coords = pickle.load(f)
                f.close()
        return modes, samples_coords, chimera2prody, prody_molecule


####
def prody_modes(molecule, max_modes, algorithm=None, **options):
    """
    Parameters
    ----------
    molecule : prody.AtomGroup
    nax_modes : int
        number of modes to calculate
    algorithm : callable, optional, default=None
        coarseGrain(prm) wich make molecule.select().setBetas(i) where i
        is the index Coarse Grain group
        Where prm is prody AtomGroup
    options : dict, optional
        Parameters for algorithm callable
    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    modes = None
    if algorithm in ['residues', 'mass']:
        title = 'normal modes for {}'.format(molecule.getTitle())
        molecule = GROUPERS[algorithm](molecule, **options)
        modes = prody.RTB(title)
        modes.buildHessian(molecule.getCoords(), molecule.getBetas())
        modes.calcModes(n_modes=max_modes)
    elif algorithm == 'calpha':
        calphas_modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        calphas = molecule = molecule.select(algorithm)
        calphas_modes.buildHessian(calphas)
        calphas_modes.calcModes(n_modes=max_modes)
        modes = prody.extendModel(calphas_modes, calphas, molecule, norm=True)[0]
    else:
        return
        modes = prody.ANM('normal modes for {}'.format(molecule.getTitle()))
        modes.buildHessian(molecule)
        modes.calcModes(n_modes=max_modes)
    return modes

def pca_modes(molecule, trajectory, max_modes, pca_atoms):
    modes = None
    traj = prody.parseDCD(trajectory) #trajectory file (.dcd)
    if pca_atoms == "calpha":
        traj.setAtoms(molecule.calpha) 
        traj.setCoords(molecule)
        traj.superpose()
        pca = prody.EDA('normal modes for {}'.format(molecule.getTitle())) #EDA is better for MD trajectories
        pca.buildCovariance(traj)
        pca.calcModes(n_modes=max_modes)
        modes = prody.extendModel(pca, molecule.calpha, molecule, norm=True)[0]
    elif pca_atoms == "backbone":
        traj.setAtoms(molecule.backbone) 
        traj.setCoords(molecule)
        traj.superpose()
        pca = prody.EDA('normal modes for {}'.format(molecule.getTitle())) #EDA is better for MD trajectories
        pca.buildCovariance(traj)
        pca.calcModes(n_modes=max_modes)
        modes = prody.extendModel(pca, molecule.backbone, molecule, norm=True)[0]
    else:
        traj.setAtoms(molecule) 
        traj.setCoords(molecule)
        traj.superpose()
        pca = prody.EDA('normal modes for {}'.format(molecule.getTitle())) #EDA is better for MD trajectories
        pca.buildCovariance(traj)
        pca.calcModes(n_modes=max_modes)
    return modes
    
def gaussian_modes(path):
    """
    Read the modes
    Create a prody.modes instance
    Parameters
    ----------
    path : str
        gaussian frequencies output path
    Returns
    -------
    modes : ProDy modes ANM or RTB
    """
    gaussian_parser = Gaussian(path).parse()
    shape = gaussian_parser.vibdisps.shape
    modes_vectors = gaussian_parser.vibdisps.reshape(shape[0], shape[1]*shape[2]).T
    modes_frequencies = numpy.abs(gaussian_parser.vibfreqs)
    modes = prody.NMA()
    modes.setEigens(vectors=modes_vectors, values=modes_frequencies)
    return modes


def convert_chimera_molecule_to_prody(molecule):
    """
    Function that transforms a chimera molecule into a prody atom group
    Parameters
    ----------
    molecule : chimera.Molecule
    Returns
    -------
    prody_molecule : prody.AtomGroup()
    chimera2prody : dict
        dictionary: chimera2prody[chimera_atom.coordIndex] = i-thm element prody getCoords() array
    """
    prody_molecule = prody.AtomGroup()
    try:
        coords, elements, serials = [], [], []
        names, resnums, resnames = [], [], []
        chids, betas, masses = [], [], []
        chimera2prody = {}
        offset_chimera_residue = min(r.id.position for r in molecule.residues)

        for i, atm in enumerate(molecule.atoms):
            chimera2prody[atm.serialNumber] = i
            coords.append(tuple(atm.coord()))  # array documentation to improve
            elements.append(atm.element.name)
            serials.append(atm.serialNumber)
            names.append(atm.name)
            resnums.append(atm.residue.id.position - offset_chimera_residue)
            resnames.append(atm.residue.type)
            chids.append(atm.residue.id.chainId)
            masses.append(atm.element.mass)
            betas.append(atm.bfactor)

        prody_molecule.setCoords(coords)
        prody_molecule.setElements(elements)
        prody_molecule.setSerials(serials)
        prody_molecule.setNames(names)
        prody_molecule.setResnums(resnums)
        prody_molecule.setResnames(resnames)
        prody_molecule.setChids(chids)
        prody_molecule.setBetas(betas)
        prody_molecule.setMasses(masses)
        prody_molecule.setTitle(str(molecule.name))
        prody_molecule.setBonds([(chimera2prody[bond.atoms[0].serialNumber],
                                  chimera2prody[bond.atoms[1].serialNumber]) for bond in molecule.bonds])

    except AttributeError:
        raise TypeError('Attribute not found. Molecule must be a chimera.Molecule')

    return prody_molecule, chimera2prody


def group_by_residues(molecule, n=15):
    """
    Coarse Grain Algorithm 1: groups per residues
    Parameters
    ----------
    molecule : prody.AtomGroup
    n : int, optional, default=7
        number of residues per group
    Returns
    ------
    molecule : prody.AtomGroup
        New betas added
    """
    group = 1
    for chain in molecule.iterChains():
        residues_indices = sorted(list(set(chain.getResnums())))
        chain_name = chain.getChid()
        for a, b in chunker(len(residues_indices), n):
            try:
                start, end = residues_indices[a-1], residues_indices[b-1]
                selector = 'chain {} and resnum {} to {}'.format(chain_name, start, end)
                selection = molecule.select(selector)
                selection.setBetas(group)
                group += 1
            except AttributeError as e:
                logger.warning(str(e))
    return molecule


def group_by_mass(molecule, n=100):
    """
    Coarse Grain Algorithm 2: groups per mass percentage
    Parameters
    ----------
    molecule : prody.AtomGroup
    n: int, optional, default=100
        Intended number of groups. The mass of the system will be divided by this number,
        and each group will have the corresponding proportional mass. However, the final
        number of groups can be slightly different.
    Returns
    -------
    molecule: prody.AtomGroup
        New Betas added
    """
    group = 1

    total_mass = sum(molecule.getMasses())
    chunk_mass = total_mass/n

    for chain in molecule.iterChains():
        selection = molecule.select('chain {}'.format(chain.getChid()))
        mass_accumulator = 0.

        for atom in iter(selection):
            atom.setBeta(group)
            mass_accumulator += atom.getMass()
            if mass_accumulator > chunk_mass:
                mass_accumulator = 0.
                group += 1
        group += 1
    return molecule


def alg3(moldy, max_bonds=3, **kwargs):
    """
    TESTS PENDING!
    Coarse Grain Algorithm 3: Graph algorithm.
        New group when a vertice: have more than n,
                                  have 0 edges
                                  new chain
    Parameters
    ----------
    moldy : prody.AtomGroup
    n : int, optional, default=2
        maximum bonds number
    Returns
    -------
    moldy: prody.AtomGroup
        New Betas added
    """
    group = 1

    for chain in moldy.iterChains():
        selection = moldy.select('chain {}'.format(chain.getChid()))
        for atom in iter(selection):
            atom.setBeta(group)
            if atom.numBonds() >= max_bonds:
                group += 1
        group += 1
    return moldy


def chunker(end, n):
    """
    divide end integers in closed groups of n
    """
    for i in range(0, end-n+1, n):
        yield i+1, i+n
    if end % n:
        yield end-end % n+1, end


def chimeracoords2numpy(molecule):
    """
    Parameters
    ----------
    molecule : chimera.molecule
    Returns
    -------
    numpy.array with molecule.atoms coordinates
    """
    return get_atom_coordinates(molecule.atoms, transformed=False)


GROUPERS = {
    'residues': group_by_residues,
    'mass': group_by_mass,
    'calpha': 'calpha',
    '': None
}
