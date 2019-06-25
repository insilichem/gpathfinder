.. GPathFinder: Identification of ligand binding pathways 
.. by a multi-objective genetic algorithm

   https://github.com/insilichem/gaudi/tree/gpathfinder

   Copyright 2019 José-Emilio Sánchez Aparicio, Giuseppe Sciortino,
   Daniel Villadrich Herrmannsdoerfer, Pablo Orenes Chueca, 
   Jaime Rodríguez-Guerra Pedregal and Jean-Didier Maréchal
   
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.


GPathFinder
===========

.. image:: https://travis-ci.org/insilichem/gaudi.svg?branch=master
    :target: https://travis-ci.org/insilichem/gaudi

.. image:: https://readthedocs.org/projects/gpathfinder/badge/?version=latest
    :target: https://gpathfinder.readthedocs.io/en/latest/

.. image:: https://anaconda.org/insilichem/gaudi/badges/installer/conda.svg
    :target: https://anaconda.org/josan_bcn/gpathfinder

.. image:: https://img.shields.io/badge/python-2.7.16-blue.svg
   :target: https://www.python.org/downloads/release/python-2716

.. image:: https://img.shields.io/github/license/josan82/gpathfinder.svg?color=orange
   :target: http://www.apache.org/licenses/LICENSE-2.0

.. image:: https://img.shields.io/static/v1.svg?label=platform&message=linux%20|%20macOS&color=lightgrey

.. .. image:: https://img.shields.io/badge/doi-10.1002%2Fjcc.24847-blue.svg
..   :target: http://onlinelibrary.wiley.com/doi/10.1002/jcc.24847/full

GPathFinder is an extension built over GaudiMM core to allow the identification 
of ligand binding pathways at atomistic level.

.. image:: docs/data/gpathfinderlogo-whitebg.jpg
    :alt: GPathFinder logo

Features
--------

**Different options for generate pathways**

- Unbinding routes from a known binding site
- Binding routes to a known binding site
- Channel analysis (given starting and final points)

**Flexibility for the ligand**

- Dihedral angles

**Different levels of flexibility for the receptor**

- Side-chain flexibility using rotamer libraries
- Global movements by Normal Mode Analysis sampling

**Different options for evalute and optimize the solutions**

- Steric clashes
- Vina scoring function
- Smoothness of the ligand movements

Documentation and support
-------------------------

Documentation source is available in ``docs/`` subdirectory, and also compiled as HTML at `this webpage <https://gpathfinder.readthedocs.io/en/latest/>`_.

If you need help with GPathFinder, please use the `issues page <https://github.com/insilichem/gaudi/issues>`_ of our `GitHub repo <https://github.com/insilichem/gaudi>`_. You can drop me a message at `joseemilio.sanchez@uab.cat <mailto:joseemilio.sanchez@uab.cat>`_ too.

**Developer friendly**

If the provided genes and objectives are not enough, you can always code your own ones. Check out the `developer docs <https://gpathfinder.readthedocs.io/en/latest/developers.html>`_!

License
-------

GPathFinder and GaudiMM are licensed under the Apache License, Version 2.0. Check the details in the `LICENSE <https://raw.githubusercontent.com/insilichem/gaudi/master/LICENSE>`_ file.

History of versions
-------------------

- **v1.0.0:** Release version. Used in the benchmark and cases study of the article.

OS Compatibility
----------------

GPathFinder is compatible with Linux and macOS.

If you find some dificulties when installing it in a concrete distribution, please use the `issues page <https://github.com/insilichem/gaudi/issues>`_ to report them.
