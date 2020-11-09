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
The GPathFinder package is comprised of several core modules that establish the base
architecture to build an extensible platform of molecular design.

The main module is :mod:`gpath.base`, which defines the :class:`gpath.base.Individual`,
whose instances represent the potential solutions to the proposed problem. Two plugin
packages allow easy customization of how individuals are defined (:mod:`gpath.genes`) and
how they are evaluated (:mod:`gpath.objectives`). Additionally:

- :mod:`gpath.algorithms` is the place to look for the actual GA implementation
- :mod:`gpath.box` is a placeholder for several small functions that are used across GPathFinder.
- :mod:`gpath.exceptions` defines custom exceptions.
- :mod:`gpath.parse` contains parsing utilities to retrieve the configuration files.
- :mod:`gpath.parallel` contains helpers to deal with parallel GPathFinder jobs.
- :mod:`gpath.plugin` holds some magic to make the plugin system work.
- :mod:`gpath.similarity` defines the diversity enhancers.
- :mod:`gpath.path_similarity` defines the diversity enhancers for pathways in GPathFinder.
"""

# Logging
import logging
try:  # Python 2.7+
    from logging import NullHandler
except ImportError:
    class NullHandler(logging.Handler):

        def emit(self, record):
            pass

logging.getLogger(__name__).addHandler(NullHandler())

__author__ = 'José-Emilio Sánchez Aparicio, Giuseppe Sciortino, \
              Daniel Villadrich Herrmannsdoerfer, Pablo Orenes Chueca, \
              Jaime Rodríguez-Guerra Pedregal and Jean-Didier Maréchal'
__copyright__ = '2019, InsiliChem'
__url__ = 'https://github.com/insilichem/gpathfinder'
__title__ = 'GPathFinder'
__description__ = ('Identification of ligand binding pathways \
                    by a multi-objective genetic algorithm')
__version__ = '1.3.0'
