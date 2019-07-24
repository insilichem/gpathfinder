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
This modules allows to explore side chains flexibility
in proteins, as well as mutation.

It needs that at least a :class:`gpath.genes.rotamers.molecule.Molecule` has been
requested in the input file. Residues of those are referenced in the `residues` argument.

It also needs a GPathFinder `path` gene.

"""

# Python
import random
from collections import OrderedDict
import logging
# Chimera
import chimera
import Midas
from chimera import BondRot, dihedral
from chimera.phipsi import chiAtoms, AtomsMissingError
from Rotamers import getRotamerParams, getRotamers, NoResidueRotamersError
# External dependencies
import deap.tools
# GPATH
from gpath import parse, box
from gpath.genes import GeneProvider

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Rotamers.validate(kwargs)
    return Rotamers(**kwargs)


class Rotamers(GeneProvider):

    """
    Rotamers class

    Parameters
    ----------
    residues : list of str
        Residues that should be analyzed. This has to be in the form:

            [ Protein/233, Protein/109 ]

        where the first element (before slash) is the gpath.genes.molecule name
        and the second element (after slash) is the residue position number in that
        molecule.

        This list of str is later parsed to the proper chimera.Residue objects

    library : {'Dunbrack', 'Dynameomics'}
        The rotamer library to use.

    with_original: bool, defaults to True
        Whether to include the original set of chi angles as part of the
        rotamer library.

    Attributes
    ----------
    allele : list of float
        For i residues, it contains i floats within [0, 1), that will point
        to the selected rotamer torsions for each residue.
    """
    _validate = {
        'residues': [parse.Named_spec("molecule", "residue")],
        'library': parse.Any('Dunbrack', 'dunbrack', 'Dynameomics', 'dynameomics'),
        'with_original': parse.Boolean,
       }

    # Avoid unnecesary calls to expensive get_rotamers if residue is known
    # to not have any rotamers
    _residues_without_rotamers = set(('ALA', 'GLY'))

    def __init__(self, residues=[], library='Dunbrack', with_original=True, **kwargs):
        GeneProvider.__init__(self, **kwargs)
        self._kwargs = kwargs
        self._residues = residues
        self.library = library
        self.with_original = with_original
        self._need_express = False #Control of expression by Gpath
        self.allele = []
        # set caches
        try:
            self.residues = self._cache[self.name + '_residues']
        except KeyError:
            self.residues = self._cache[self.name + '_residues'] = OrderedDict()
        try:
            self.rotamers = self._cache[self.name + '_rotamers']
        except KeyError:
            self.rotamers = self._cache[self.name + '_rotamers'] = OrderedDict()


    def __ready__(self):
        """
        Second stage of initialization.

        It parses the requested residues strings to actual residues.
        """
        self.allele = []
        self.residues = {}
        for molname, pos in self._residues:
            residues = self.parent.find_molecule(molname).find_residues(int(pos))
            if len(residues) > 1 and pos != '*':
                logger.warn('Found one more than residue for %s/%s', molname, pos)
            for r in residues:
                self.patch_residue(r)
                self.residues[(molname, r.id.position)] = r
                self.allele.append(random.random())

    def __deepcopy__(self, memo):
        new = self.__class__(residues=self._residues, library=self.library, **self._kwargs)
        new.residues = self.residues
        new.allele = self.allele[:]
        return new

    def express(self):
        if self._need_express: #Control of expression by Gpath
            for ((molname, pos), residue), i in zip(self.residues.items(), self.allele):
                if residue.type not in self._residues_without_rotamers:
                    try:
                        rotamers = self.retrieve_rotamers(molname, pos, residue,
                                                          library=self.library.title())
                    #####
                    except:  # ALA, GLY...
                        logger.warn('%s/%s (%s) has no rotamers', molname, pos, residue.type)
                        self._residues_without_rotamers.add(residue.type)
                    else:
                        rotamer_chis = rotamers[int(i * len(rotamers))]
                        self.update_rotamer(residue, rotamer_chis)

    def unexpress(self):
        if self._need_express: #Control of expression by Gpath
            for res in self.residues.values():
                for torsion in res._rotamer_torsions:
                    torsion.adjustAngle(-torsion.angle, torsion.rotanchor)

    def mate(self, mate):
        """
        Not needed in GPathFinder
        """
        pass

    def mutate(self, indpb):
        """
        Not needed in GPathFinder
        """
        pass

    def write(self, path, name, *args, **kwargs):
        """
        Ìt is not necessary to write anything by ``rotamers`` gene. All 
        the output files are generated by ``path`` gene.
        """
        return None

    def retrieve_rotamers(self, molecule, position, residue, library='Dunbrack'):
        try:
            rotamers = self.rotamers[(molecule, position)]
        except KeyError:
            def sort_by_rmsd(ref, query):
                ref_atoms = [ref.atomsMap[a.name][0] for a in query.atoms]
                return Midas.rmsd(ref_atoms, query.atoms)
            chis = getRotamerParams(residue, lib=library)[2]
            rotamers_mols = getRotamers(residue, lib=library)[1]
            reference = residue if self.with_original else rotamers_mols[0].residues[0]
            rotamers_and_chis = zip(rotamers_mols, [c.chis for c in chis])
            rotamers_and_chis.sort(key=lambda rc: sort_by_rmsd(reference, rc[0]))
            rotamers = zip(*rotamers_and_chis)[1]
            if self.with_original:
                rotamers = (self.all_chis(residue),) + rotamers
            self.rotamers[(molecule, position)] = rotamers
            for rot in rotamers_mols:
                rot.destroy()
        return rotamers

    @staticmethod
    def update_rotamer(residue, chis):
        for bondrot, chi in zip(residue._rotamer_torsions, chis):
            bondrot.adjustAngle(chi - bondrot.chi, bondrot.rotanchor)

    @staticmethod
    def patch_residue(residue):
        if getattr(residue, '_rotamer_torsions', None):
            return
        residue._rotamer_torsions = [] # BondRot objects cache
        alpha_carbon = next((a for a in residue.atoms if a.name == 'CA'), residue.atoms[0])
        for chi in range(1, 5):
            try:
                atoms = chiAtoms(residue, chi)
                bond = atoms[1].bondsMap[atoms[2]]
                br = BondRot(bond)
                br.rotanchor = box.find_nearest(alpha_carbon, bond.atoms)
                br.chi = dihedral(*[a.coord() for a in atoms])
                residue._rotamer_torsions.append(br)
            except AtomsMissingError:
                break
            except (chimera.error, ValueError) as v:
                if "cycle" in str(v) or "already used" in str(v):
                    continue  # discard bonds in cycles and used!
                break

    @staticmethod
    def all_chis(residue):
        chis = []
        for i in range(1, 5):
            try:
                chis.append(getattr(residue, 'chi{}'.format(i)))
            except AttributeError:
                break
        return chis
