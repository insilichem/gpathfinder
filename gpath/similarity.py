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
This module contains the similarity functions that are used
to discard individuals that are not different enough.

"""
from __future__ import print_function, division
import logging

logger = logging.getLogger(__name__)

def rmsd(ind1, ind2, subjects, threshold, *args, **kwargs):
    """
    Returns the RMSD between two individuals

    Parameters
    ----------
    ind1 : gpath.base.Individual
    ind2 : gpath.base.Individual
    subjects : list of str
        Name of gpath.genes.molecule instances to measure
    threshold : float
        Maximum RMSD value to consider two individuals as similar.
        If ``rmsd > threshold``, they are considered different.

    Returns
    -------
    bool
        True if ``rmsd`` is within threshold, False otherwise.
        It will always return False if number of atoms is not
        equal in the two Individuals.

    """
    molecules1 = [ind1.find_molecule(s) for s in subjects]
    molecules2 = [ind2.find_molecule(s) for s in subjects]

    logger.debug("Comparing RMSD between #%s and #%s", id(ind1), id(ind2))
    for m1, m2 in zip(molecules1, molecules2):
        coords1 = m1._expressed_coordinates
        coords2 = m2._expressed_coordinates
        if coords1.shape[0] != coords2.shape[0]:
            return False
        rmsd_squared = _rmsd_squared(coords1, coords2)
        if rmsd_squared > threshold*threshold:
            return False
    return True


def _rmsd_squared(coords1, coords2):
    diff = coords1 - coords2
    return (diff * diff).sum() / coords1.shape[0]
