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
Helper functions to deal with parallel execution of GPathFinder jobs.abs

Useful for benchmarks.
"""

from multiprocessing import Pool, cpu_count


def run_parallel(fn, args=(), processes=None,
                 initializer=None, initargs=(), maxtasksperchild=1,
                 map_timeout=9999999, map_chunksize=1, map_callback=None):
    """
    Create a pool instance with built-in exception handling.
    """

    pool = Pool(processes=processes, maxtasksperchild=maxtasksperchild)
    try:
        pool.map_async(fn, args, chunksize=map_chunksize,
                       callback=map_callback).get(map_timeout)
    except KeyboardInterrupt:
        print('Exiting...')
        pool.terminate()
    except Exception as e:
        print('An error ocurred:', type(e).__name__, e.message)
        pool.terminate()
    finally:
        pool.close()
        pool.join()