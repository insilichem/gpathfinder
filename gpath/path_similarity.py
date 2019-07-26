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
This module contains the similarity functions that are used to discard 
individuals that are not different enough. 
This criterion of similarity will be applied in the case of two 
``pathways`` individuals with the same score. Then, if they are similar
enough according to this module, one of them will be discarded.
"""
from __future__ import print_function, division
import logging
import numpy as np

logger = logging.getLogger(__name__)

def pathways_rmsd(ind1, ind2, subject, threshold, *args, **kwargs):
    """
    Calculates the RMSD between the positions of the ``pathways`` genes
    belonging two the two individuals object of study. If the squared 
    RMSD is less or equal than the squared threshold, we consider that
    the two pathways are identical and one of them will be discarded.

    Parameters
    ----------
    ind1 : gpath.base.Individual
    ind2 : gpath.base.Individual
    subject: str
        Name of Gpath ``pathway`` gene instance to measure.
    threshold : float
        Maximum RMSD value in Angstroms to consider two individuals as 
        similar.
        If ``rmsd > threshold``, they are considered different.
    
    Returns
    -------
    bool
        True if ``rmsd`` is within threshold, False otherwise.
        It will always return False if number of points of the pathway 
        is not equal in the two Individuals.
    """
    coords1 =  np.array([elem[:] for elem in \
                                    ind1.genes[subject].allele['positions']])
    coords2 =  np.array([elem[:] for elem in \
                                    ind2.genes[subject].allele['positions']])

    if coords1.shape[0] != coords2.shape[0]:
        return False
    rmsd_squared = _rmsd_squared(coords1, coords2)
    if rmsd_squared > threshold*threshold:
        return False
    return True

def _rmsd_squared(coords1, coords2):
    diff = coords1 - coords2
    return (diff * diff).sum() / coords1.shape[0]