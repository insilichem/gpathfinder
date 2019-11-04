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
# Jaime Rodríguez-Guerra Pedregal and Jean-Didier Maréchal
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
This objective is a wrapper around the scoring functions provided by
`AutoDock Vina <http://vina.scripps.edu/>`_.
"""

import logging
from gpath.objectives import ObjectiveProvider
from gpath import parse
import math
import chimera

DIST_PROBE_ALPHA = {
    'ARG': (3.913, 10.308),
    'ASN': (4.001, 5.751),
    'ASP': (4.058, 5.993),
    'CYS': (3.532, 5.115),
    'GLN': (4.141, 7.649),
    'GLU': (4.428, 7.453),
    'HIS': (4.239, 7.680),
    'MET': (4.521, 6.547),
    'SER': (3.487, 5.058),
    'THR': (3.439, 4.976),
    'TYR': (6.699, 9.151),
    'LYS': (5.279, 9.267),
    'LEU': (3.748, 4.690),
    'TRP': (3.168, 8.189),
    'ILE': (4.016, 5.025),
    'ALA': (3.218, 5.440),
    'VAL': (3.717, 4.848),
    'PHE': (3.518, 4.885),
    'PRO': (3.629, 4.667),
    'ALL': (3.390, 7.218)
}

DIST_PROBE_BETA = {
    'ARG': (3.612, 8.966),
    'ASN': (3.911, 5.157),
    'ASP': (3.675, 4.991),
    'CYS': (2.794, 3.843),
    'GLN': (3.875, 6.139),
    'GLU': (4.252, 5.949),
    'HIS': (3.318, 7.014),
    'MET': (3.625, 5.101),
    'SER': (2.817, 4.146),
    'THR': (2.816, 4.126),
    'TYR': (6.344, 8.082),
    'LYS': (4.735, 8.052),
    'LEU': (4.317, 5.727),
    'TRP': (4.267, 6.920),
    'ILE': (4.871, 5.801),
    'ALA': (3.653, 6.463),
    'VAL': (4.345, 5.669),
    'PHE': (4.111, 5.769),
    'PRO': (4.758, 5.271),
    'ALL': (2.666, 6.351)
}

ANGLE_PAB = {
    'ARG': (-0.128, 1.824),
    'ASN': (0.668, 1.695),
    'ASP': (0.229, 1.629),
    'CYS': (0.021, 1.299),
    'GLN': (-0.217, 1.723),
    'GLU': (-0.164, 1.753),
    'HIS': (0.307, 1.455),
    'MET': (-0.132, 1.210),
    'SER': (0.278, 1.408),
    'THR': (0.358, 1.416),
    'TYR': (0.558, 1.422),
    'LYS': (0.020, 1.616),
    'LEU': (1.443, 2.560),
    'TRP': (0.362, 2.422),
    'ILE': (1.728, 2.238),
    'ALA': (1.394, 2.478),
    'VAL': (1.607, 2.213),
    'PHE': (1.593, 2.254),
    'PRO': (1.464, 2.661),
    'ALL': (0.092, 1.601),
}

DIST_PROBE_OXYGEN = (2.017, 2.946)

ANGLE_POC = (1.955, 3.163)

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Metal.validate(kwargs)
    return Metal(**kwargs)


class Metal(ObjectiveProvider):

    """
    Vina class
    Parameters
    ----------
    receptor : str
        Key of the gene containing the molecule acting as receptor (protein)
    ligand : str
        Key of the gene containing the molecule acting as ligand
    metal_atoms : list of Molecule/atom_serial_number
        List of atoms of the Ligand molecule that will be considered as metals in
        order to search possible coordinating sites in the receptor.
    residues : list of str, optional
        List of residues that can act as coordinating residues (named by 3-letter
        aminoacid code). Default to [ARG,ASN,ASP,CYS,GLN,GLU,HIS,MET,SER,THR,TYR]
    backbone : bool, optional
        Also consider possible coordination with backbone atoms. Default to False
    ####COORDINATING_LIGAND_ATOMS

    Returns
    -------
    float
        Metal interaction score.
    """
    _validate = {
        parse.Required('receptor'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        parse.Required('metal_atoms') : [parse.All(parse.Named_spec("molecule", "atom"))],
        'residues' : parse.All(parse.Coerce(str)),
        'backbone' : parse.Coerce(bool),
        }

    def __init__(self, receptor='Protein', ligand='Ligand', metal_atoms=None,
                    residues = ['ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'MET', 'SER', 'THR', 'TYR'],
                    backbone = False, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.receptor = receptor
        self.ligand = ligand
        self.metal_atoms = metal_atoms
        self.residues = residues
        self.backbone = backbone

    def evaluate(self, ind):
        protein = ind.find_molecule(self.receptor).compound.mol
        ligand = ind.find_molecule(self.ligand).compound.mol
        if self.metal_atoms is not None:
            self.probes = []
            for at in self.metal_atoms:
                mol, atom = at
                try:
                    molecule_gene = ind.find_molecule(mol)
                    probe = molecule_gene.find_atom(atom)
                except StopIteration:
                    pass
                else:
                    self.probes.append(probe)
        else:
            raise Exception('You should indicate the parameter ``metal_atoms`` containing the list of the Ligand Atoms that will be considered as metals')   

        atom_surr = self._surrounding_atoms(ind)
        possible_atoms = []
        for at in atom_surr:
            if at.name in ['CA', 'CB', 'O', 'C', 'N']:
                possible_atoms.append(at)

        #Set residues to consider as coordinating
        for r in DIST_PROBE_ALPHA.keys():
            if r not in self.residues:
                del(DIST_PROBE_ALPHA[r])
                del(DIST_PROBE_BETA[r])
                del(ANGLE_PAB[r])

        atoms_to_consider = self._create_dict_atoms(possible_atoms)
        
        good_probes = 0

        for probe in self.probes: # METALS TO CONSIDER
            probe_coords = probe.xformCoord()
            for res in atoms_to_consider.keys():
                # Check and discard possible clashes with backbone atoms
                try:
                    dist_oxygen = atoms_to_consider[res]['oxygen'].xformCoord().distance(probe_coords)
                    if dist_oxygen < 1.0:
                        good_probes = 0
                        break
                    dist_carbon = atoms_to_consider[res]['carbon'].xformCoord().distance(probe_coords)
                    if dist_carbon < 1.0:
                        good_probes = 0
                        break
                    dist_nitrogen = atoms_to_consider[res]['nitrogen'].xformCoord().distance(probe_coords)
                    if dist_nitrogen < 1.0:
                        good_probes = 0
                        break
                    dist_alpha = atoms_to_consider[res]['alpha'].xformCoord().distance(probe_coords)
                    if dist_alpha < 1.0:
                        good_probes = 0
                        break
                except:
                    pass

                # include sidechain coordinations
                # To be a potential sidechain coordination, needs to accomplish:
                # 1 - Range of distance probe - carbon alpha
                # 2 - Range of distance probe - carbon beta
                # 3 - Range of angle probe - carbon alpha - carbon beta

                res_name = str(res.type)
                # now res_name contain the name of the residue that is needed 
                # to check in the ranges tables (DIST_PROBE_ALPHA, 
                # DIST_PROBE_BETA, ANGLE_PAB)
                try:
                    if (res_name in DIST_PROBE_ALPHA.keys()):
                        if (DIST_PROBE_ALPHA[res_name][0] <= dist_alpha <= DIST_PROBE_ALPHA[res_name][1]):
                            dist_beta = atoms_to_consider[res]['beta'].xformCoord().distance(probe_coords)
                            if (DIST_PROBE_BETA[res_name][0] <= dist_beta <= DIST_PROBE_BETA[res_name][1]):                        
                                dist_alpha_beta = atoms_to_consider[res]['alpha'].xformCoord().distance(atoms_to_consider[res]['beta'].xformCoord())
                                angle_PAB = math.acos((dist_alpha ** 2 + dist_alpha_beta ** 2 - dist_beta ** 2)/(2*dist_alpha*dist_alpha_beta))
                                if (ANGLE_PAB[res_name][0] <= angle_PAB <= ANGLE_PAB[res_name][1]):  
                                    good_probes += 1
                except:
                    continue

                # include backbone coordinations
                if self.backbone:
                    try:
                        if DIST_PROBE_OXYGEN[0] <= dist_oxygen <= DIST_PROBE_OXYGEN[1]:
                            dist_oxygen_carbon = atoms_to_consider[res]['oxygen'].xformCoord().distance(atoms_to_consider[res]['carbon'].xformCoord())
                            angle_POC = math.acos((dist_oxygen ** 2 + dist_oxygen_carbon ** 2 - dist_carbon ** 2)/(2*dist_oxygen*dist_oxygen_carbon))
                            if ANGLE_POC[0] <= angle_POC <= ANGLE_POC[1]:
                                good_probes += 1
                    except:
                        continue
        return good_probes

    def _surrounding_atoms(self, ind):
        """
        Get atoms in the search zone, based on the molecule, (possible) rotamer
        genes and the radius
        """
        self.zone.clear()
        #Add all atoms of probes molecules
        self.zone.add(self.probes)

        #Surrounding zone from probes
        self.zone.merge(chimera.selection.REPLACE,
                        chimera.specifier.zone(self.zone, 'atom', None,
                                               10.5, [ind.find_molecule(self.receptor).compound.mol]))
        return self.zone.atoms()
    
    def _create_dict_atoms(self, atoms):
        atoms_to_consider = {}
        for at in atoms:
            if at.name == 'CA':
                if at.residue in atoms_to_consider.keys():
                    atoms_to_consider[at.residue]['alpha'] = at
                else:
                    atoms_to_consider[at.residue] = {'alpha': at}
            elif at.name == 'CB':
                if at.residue in atoms_to_consider.keys():
                    atoms_to_consider[at.residue]['beta'] = at
                else:
                    atoms_to_consider[at.residue] = {'beta': at}
            elif at.name == 'O':
                if at.residue in atoms_to_consider.keys():
                    atoms_to_consider[at.residue]['oxygen'] = at
                else:
                    atoms_to_consider[at.residue] = {'oxygen': at}
            elif at.name == 'C':
                if at.residue in atoms_to_consider.keys():
                    atoms_to_consider[at.residue]['carbon'] = at
                else:
                    atoms_to_consider[at.residue] = {'carbon': at}
            elif at.name == 'N':
                if at.residue in atoms_to_consider.keys():
                    atoms_to_consider[at.residue]['nitrogen'] = at
                else:
                    atoms_to_consider[at.residue] = {'nitrogen': at}
        return atoms_to_consider