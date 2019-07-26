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
`gpath.cli.gpath_run` is the main hub for launching GPATH jobs.

It sets up the configuration environment needed by DEAP (responsible for the GA)
and ties it up to the GPATH custom classes that shape up the invididuals and
objectives. All in a loosely-coupled approach based on Python modules called on-demand.

**Usage**. Simply, type:

.. code-block :: console

    gpath run /path/to/job.gpath-input
"""

# Python
from __future__ import print_function
from importlib import import_module
import logging
import os
import sys
# External dependencies
try:
    import chimera
except ImportError:
    sys.exit('Chimera not importable from this environment. Please, install '
             'pychimera with `pip install pychimera`.')
import numpy
import deap.creator
import deap.tools
import deap.base
import deap.algorithms
# from multiprocess import Pool
# GPATH
import gpath.algorithms
import gpath.base
import gpath.box
import gpath.genes
import gpath.objectives
import gpath.parse
import gpath.plugin
import gpath.similarity

def launch(cfg):
    """
    Runs a GPATH job

    Parameters
    ----------
    cfg : gpath.parse.Settings
        Parsed YAML dict with attribute-like access
    """
    gpath.plugin.import_plugins(*cfg.genes)
    gpath.plugin.import_plugins(*cfg.objectives)
    import_module(cfg.similarity.module.rsplit('.', 1)[0])

    # DEAP setup: Fitness, Individuals, Population
    toolbox = deap.base.Toolbox()
    toolbox.register("call", (lambda fn, *args, **kwargs: fn(*args, **kwargs)))
    toolbox.register("individual", toolbox.call, gpath.base.MolecularIndividual, cfg)
    toolbox.register("population", deap.tools.initRepeat, list, toolbox.individual)
    population = toolbox.population(n=cfg.ga.population)

    environment = gpath.base.Environment(cfg)

    toolbox.register("evaluate", lambda ind: environment.evaluate(ind))
    toolbox.register("mate", (lambda ind1, ind2: ind1.mate(ind2)))
    toolbox.register("mutate", (lambda ind, indpb: ind.mutate(indpb)), indpb=cfg.ga.mut_indpb)
    toolbox.register("similarity", (lambda ind1, ind2: ind1.similar(ind2)))
    toolbox.register("select", deap.tools.selNSGA2)
    
    # Multiprocessing pool (Resolve pickle!)
    # pool = Pool(initializer=lambda: setattr(sys, 'stdout', open(str(os.getpid()) + ".out", "w")))
    # toolbox.register("map", pool.map)

    if cfg.output.history:
        history = deap.tools.History()
        # Decorate the variation operators
        toolbox.decorate("mate", history.decorator)
        toolbox.decorate("mutate", history.decorator)
        history.update(population)

    if cfg.output.pareto:
        elite = deap.tools.ParetoFront(lambda ind1, ind2: False)
    else:
        elite_size = int(cfg.ga.population * cfg.ga.mu)
        elite = MyHallOfFame(elite_size)
    
    if cfg.output.verbose:
        stats = deap.tools.Statistics(lambda ind: ind.fitness.values)
        numpy.set_printoptions(precision=cfg.output.precision)
        stats.register("avg", numpy.mean, axis=0)
        stats.register("std", numpy.std, axis=0)
        stats.register("min", numpy.min, axis=0)
        stats.register("max", numpy.max, axis=0)
    else:
        stats = None

    # Begin evolution
    mu = int(cfg.ga.mu * cfg.ga.population)
    lambda_ = int(cfg.ga.lambda_ * cfg.ga.population)
    population, log = gpath.algorithms.ea_mu_plus_lambda(
        population, toolbox, cfg=cfg, mu=mu, lambda_=lambda_,
        cxpb=cfg.ga.cx_pb, mutpb=cfg.ga.mut_pb, 
        ngen=cfg.ga.generations, halloffame=elite,
        verbose=cfg.output.verbose, stats=stats,
        prompt_on_exception=cfg.output.prompt_on_exception)

    return population, log, elite


def enable_logging(path=None, name=None, debug=False):
    """
    Register loggers and handlers for both stdout and file
    """

    class CustomFormatter(logging.Formatter):

        CUSTOM_FORMATS = {
            logging.DEBUG: "DEBUG: %(module)s: %(lineno)d: %(message)s",
            logging.INFO: "%(message)s",
            logging.WARNING: "Warning: %(message)s",
            logging.ERROR: "[!] %(message)s",
            logging.CRITICAL: "CRITICAL: %(message)s",
            100: "%(message)s"
        }

        def format(self, record):
            format_orig = self._fmt
            self._fmt = self.CUSTOM_FORMATS.get(record.levelno, format_orig)
            result = logging.Formatter.format(self, record)
            self._fmt = format_orig
            return result

    logger = logging.getLogger('gpath')
    logger.setLevel(logging.DEBUG)

    # create CONSOLE handler and set level to error
    handler = logging.StreamHandler()
    handler.setLevel(logging.ERROR)
    formatter = CustomFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    # create debug file handler and set level to debug
    if path and name:
        handler = logging.FileHandler(os.path.join(path, name + ".gpath-log"), 'w')
        if debug:
            handler.setLevel(logging.DEBUG)
        else:
            handler.setLevel(15)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "%Y.%m.%d %H:%M:%S")
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


def unbuffer_stdout():
    class Unbuffered(object):

        def __init__(self, stream):
            self.stream = stream

        def write(self, data):
            self.stream.write(data)
            self.stream.flush()

        def __getattr__(self, attr):
            return getattr(self.stream, attr)

    sys.stdout = Unbuffered(sys.stdout)


#@gpath.box.do_cprofile
def main(cfg, debug=False):
    """
    Starts a GPATH job

    Parameters
    ----------
    cfg : str or gpath.parse.Settings
        Path to YAML input file or an already parsed YAML file
        via gpath.parse.Settings class
    debug : bool, optional, default=False
        Whether to enable verbose logging or not.
    """
    # Parse input file
    if isinstance(cfg, basestring) and os.path.isfile(cfg):
        cfg = gpath.parse.Settings(cfg)

    # Enable logging to stdout and file
    unbuffer_stdout()
    logger = enable_logging(cfg.output.path, cfg.output.name, debug=debug)
    logger.log(100, 'Loaded input %s', cfg._path)

    # Check names of the modules (should be gpath)
    for m in cfg.genes + cfg.objectives + [cfg.similarity]:
        if 'gpath' not in m.module:
            logger.error("From version 1.0.1 all genes, objectives and similarity " \
                        "modules should be imported from gpath. Please, check your " \
                        ".yaml file (sections: genes, objectives and similarity) and " \
                        "change old references to gaudi modules by gpath. For example: " \
                        "gaudi.genes.molecule should be gpath.genes.molecule")
            sys.exit()
    
    # Round population to multiple of 4 if necessary
    if cfg.ga.population % 4 > 0:
        cfg.ga.population = int(4 * round(float(cfg.ga.population)/4))
        logger.log(100, 'Population rounded to: %s', str(cfg.ga.population))

    # Place a copy of input inside output path
    with open(os.path.join(cfg.output.path, os.path.basename(cfg._path)), 'w') as f:
        f.write(cfg.toYAML())

    # Disable Chimera's auto ksdssp
    chimera.triggers.addHandler("Model", gpath.box.suppress_ksdssp, None)

    # Run simulation
    try:
        logger.log(100, 'Launching job with...')
        logger.log(100, '  Genes: %s', ', '.join([g.name for g in cfg.genes]))
        logger.log(100, '  Objectives: %s', ', '.join([o.name for o in cfg.objectives]))
        pop, log, best = launch(cfg)
    except Exception as e:
        log_path = os.path.join(cfg.output.path, cfg.output.name + ".gpath-log")
        logger.error('An exception occurred: %s\n    '
                     'Check traceback in logfile %s', e, log_path)
        logger.log(15, "An exception occurred", exc_info=True)
        sys.exit(1)

    gpath.algorithms.dump_population(best, cfg)
    logger.handlers = []

#Class HallOfFame with a bug fixed in the update method
class MyHallOfFame(deap.tools.HallOfFame):
   def update(self, population):
        """Update the hall of fame with the *population* by replacing the
        worst individuals in it by the best individuals present in
        *population* (if they are better). The size of the hall of fame is
        kept constant.
        
        :param population: A list of individual with a fitness attribute to
                           update the hall of fame with.
        """     
        for ind in population:
            if len(self) == 0 and self.maxsize !=0:
                # Working on an empty hall of fame is problematic for the
                # "for else"
                self.insert(population[0])
                continue
            if ind.fitness > self[-1].fitness or len(self) < self.maxsize:
                for hofer in self:
                    # Loop through the hall of fame to check for any
                    # similar individual
                    if self.similar(ind, hofer):
                        break
                else:
                    # The individual is unique and strictly better than
                    # the worst
                    if len(self) >= self.maxsize:
                        self.remove(-1)
                    self.insert(ind)