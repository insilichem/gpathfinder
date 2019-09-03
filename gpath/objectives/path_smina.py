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
`Smina <https://sourceforge.net/projects/smina/>`_.
"""

import logging
import os
import sys
from subprocess import call, check_output, CalledProcessError
from tempfile import _get_default_tempdir, _get_candidate_names
from gpath.objectives import ObjectiveProvider
from gpath import parse

logger = logging.getLogger(__name__)


def enable(**kwargs):
    kwargs = Smina.validate(kwargs)
    return Smina(**kwargs)


class Smina(ObjectiveProvider):

    """
    Smiina class
    Parameters
    ----------
    receptor : str
        Key of the gene containing the molecule acting as receptor (protein)
    ligand : str
        Key of the gene containing the molecule acting as ligand
    scoring : str, optional
        Specify alternative builtin scoring function (e.g. vinardo). 
        Defaults to None
    custom_scoring : str, optional
        Specify a custom scoring function file. Defaults to None
    custom_atoms : str, optional
        Specify a custom atom type parameters file. Defaults to None
    Returns
    -------
    float
        Interaction energy in kcal/mol, as reported by Smina--score-only.
    """
    _validate = {
        parse.Required('receptor'): parse.Molecule_name,
        parse.Required('ligand'): parse.Molecule_name,
        'scoring': str,
        'custom_scoring': parse.RelPathToInputFile(),
        'custom_atoms': parse.RelPathToInputFile(),
        }

    def __init__(self, receptor='Protein', ligand='Ligand', scoring=None,
                custom_scoring=None, custom_atoms=None, *args, **kwargs):
        ObjectiveProvider.__init__(self, **kwargs)
        self.receptor = receptor
        self.ligand = ligand
        self.scoring = scoring
        self.custom_scoring = custom_scoring
        self.custom_atoms = custom_atoms
        self._paths = {}
        self._tmpfile = None
        if os.name == 'posix' and os.path.exists('/dev/shm'):
            self.tmpdir = '/dev/shm'
        else:
            self.tmpdir = _get_default_tempdir()

    def evaluate(self, ind):
        """
        Run a subprocess calling Smina binary with provided options,
        and parse the results. Clean tmp files at exit.
        """
        receptor = ind.find_molecule(self.receptor)
        ligand = ind.find_molecule(self.ligand)

        receptor_path = self.prepare_receptor(receptor)
        ligand_path = self.prepare_ligand(ligand)
        command = ['smina', '--score_only', '--cpu', '1',
                   '-r', receptor_path, '-l', ligand_path]
        if self.scoring:
            command.extend(['--scoring', self.scoring])
        if self.custom_scoring:
            command.extend(['--custom_scoring', self.custom_scoring])
        if self.custom_atoms:
            command.extend(['--custom_atoms', self.custom_atoms])

        try:
            stream = check_output(command, universal_newlines=True)
        except CalledProcessError:
            logger.warning("Could not run Smina with command %s", command)
            return -100000 * self.weight
        else:
            return self.parse_output(stream)
        finally:
            self.clean()

    def prepare_receptor(self, molecule):
        receptor_path = '{}_receptor.pdb'.format(self.tmpfile)
        pdb = molecule.write(absolute=receptor_path, filetype='pdb')
        self._paths['receptor'] = receptor_path
        return receptor_path

    def prepare_ligand(self, molecule):
        ligand_path = '{}_ligand.pdb'.format(self.tmpfile)
        pdb = molecule.write(absolute=ligand_path, filetype='pdb')
        self._paths['ligand'] = ligand_path
        return ligand_path


    def parse_output(self, stream):
        for line in stream.splitlines():
            if line[:9] == "Affinity:":
                return float(line.split()[1])
        return -1000 * self.weight

    def clean(self):
        for p in self._paths.values():
            os.remove(p)
        self._paths = {}

    @property
    def tmpfile(self):
        if self._tmpfile is None:
            self._tmpfile = os.path.join(self.tmpdir, next(_get_candidate_names()))
        return self._tmpfile