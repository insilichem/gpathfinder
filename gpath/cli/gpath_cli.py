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
``gpath.cli.gpath_cli`` is the CLI entry point for all GPathFinder scripts.

Available commands:

    - run
    - view
    - prepare

"""
# Python
import os
import sys
import time
from datetime import timedelta
from importlib import import_module
from textwrap import dedent
# 3rd party
import click
# insilichem
import pychimera
import gpath


# Helpers
def test_import(command, module):
    try:
        imported = import_module('gpath.cli.' + module)
    except ImportError as e:
        sys.exit("ERROR: {} (and its dependencies) must be installed"
                 " to use <{}> command.\n{}".format(module, command, e))
    else:
        return imported


def timeit(func, *args, **kwargs):
    def wrapped():
        ts = time.time()
        result = func(*args, **kwargs)
        te = time.time()
        click.echo(
            'Finished after {:0>8}'.format(timedelta(seconds=int(te-ts))))
        return result
    return wrapped


def load_chimera(nogui=True):
    pychimera.patch_environ(nogui=nogui)
    pychimera.enable_chimera()


def echo_banner():
    logo = dedent('''
    .g8"""""""bgd `7MM"""Mq.         mm   `7MM        
  .dP'         `M   MM   `MM.        MM     MM        
  dM'           `   MM   ,M9 ,6"Yb.mmMMmm   MMpMMMb.  
  MM                MMmmdM9 8)   MM  MM     MM    MM  
 MM                 MM       ,pm9MM  MM     MM    MM  
 MM                 MM      8M   MM  MM     MM    MM  
 MM                JMML.    `Moo9^Yo.`Mbmo.JMML  JMML.
 MM          `7MMF _____ ____ _   _ ____  _____ ____  
  MM           MM |  ___|_  _| \ | |  _ \|  ___|  _ \ 
  MM.          MM | |_    || |  \| | | | | |_  | |_) |
  `Mb.         MM |  _|  _||_|     | |_| |  _| |  _ / 
    `"bmmmmmmdPY  |_|   |____|_|\__|____/|____\|_| \_\ ''')
    banner = ["{}\n{}".format(logo, '-'*len(logo.splitlines()[-1])),
              'GPathFinder: Identification of ligand pathways by a\n'
              'multi-objective genetic algorithm',
              '{} · v{}\n'.format(gpath.__copyright__, gpath.__version__)]
    return '\n'.join(banner)


@click.group(context_settings=dict(help_option_names=['-h', '--help']))
@click.version_option(version=gpath.__version__)
def cli(prog_name='gpath'):
    """
    GPathFinder: Identification of ligand pathways by 
    a multi-objective genetic algorithm

    \b
    (C) 2019, InsiliChem
    https://github.com/insilichem/gpathfinder
    """
    if not os.environ.get('CHIMERA'):
        click.echo(echo_banner())


@cli.command()
@click.option('--debug', help='Dump debug info to logfile', is_flag=True)
@click.argument('filename', required=True, type=click.Path(exists=True))
def run(filename, debug):
    """
    Launch a GPATH input file.
    """
    load_chimera()
    gpath_run = test_import('run', 'gpath_run')
    ts = time.time()
    gpath_run.main(filename, debug)
    te = time.time()
    click.echo('Finished after {:0>8}'.format(timedelta(seconds=int(te-ts))))


@cli.command()
@click.option('--viewer', required=False, type=click.Choice(['gpathview']))
@click.argument('filename', required=True, type=click.Path(exists=True))
def view(filename, viewer):
    """
    Analyze the results in a GPATH output file.
    """
    load_chimera(nogui=False)
    gpath_view = test_import('view', 'gpath_view')
    gpath_view.launch(filename, viewer=viewer)


@cli.command()
@click.argument('filename', required=False, type=click.Path())
def prepare(filename):
    """
    Create or edit a GPATH input file.
    """
    # gpathinspect = test_import('prepare', 'gpathinspect')
    # gpathinspect.cli.prepare_input(filename)
    click.echo("This argument is still unimplemented.")


if "__main__" == __name__:
    cli(prog_name='gpath')
