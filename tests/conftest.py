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

from pychimera import patch_environ, enable_chimera
patch_environ()
enable_chimera()

import os
import pytest
from gpath.base import expressed, MolecularIndividual, Environment
from gpath.box import suppress_ksdssp
import chimera

TESTPATH = os.path.dirname(os.path.abspath(__file__))
chimera.triggers.addHandler("Model", suppress_ksdssp, None)


def datapath(path):
    return os.path.join(TESTPATH, 'data', path)


@pytest.fixture
def individual():
    individual = MolecularIndividual(dummy=True)
    yield individual
    individual.clear_cache()
    chimera.closeSession()

individual2 = individual

@pytest.fixture
def environment():
    environment = Environment()
    yield environment
    environment.clear_cache()
