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
:mod:`gpath.cli.gpath_view` is a wrapper around several GUI programs
that can help visualize GPathFinder results.

As of now, it implements:

- GPathView: An extension for UCSF Chimera.

"""

from __future__ import print_function
import os
import sys
from subprocess import call
from pychimera.pychimera import guess_chimera_path 


def launch(filename, viewer=None):
    """
    if viewer in (None, 'gpathview'):
        visualize_with_gpathview(filename)
    else:
        sys.exit("Viewer {} not supported".format(viewer))
    """
    sys.exit("GPathView is still not implemented")

def visualize_with_gpathview(filename):
    chimera_path = 'chimera'
    chimera_paths = guess_chimera_path(common_locations=True)
    for path in chimera_paths:
        if 'headless' not in path.lower():
            chimera_path = os.path.join(path, 'bin', 'chimera')
    call([chimera_path, filename])
