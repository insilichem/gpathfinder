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
This script complements the Gpath gene by generating a refined pathway from the 
frames it provides. The new interpolated frames avoid high clashscore conformations
through a RRT path planning algorithm.
"""

import time
import random
import math
import copy
import heapq
import shutil
import yaml
import zipfile
from argparse import ArgumentParser

import pychimera
pychimera.patch_environ()
pychimera.enable_chimera()

from Combine import combine
import chimera
import mdtraj
import numpy as np
import Matrix as M
import itertools, operator
import os

from FitMap.search import random_rotation
from Molecule import atom_positions

import pprint
pp = pprint.PrettyPrinter(4)

import gpath
from gpath import base
from gpath.genes import path_normalmodes
from gpath.genes import molecule
from gpath.genes import path_torsion as torsion
from gpath.genes import path_rotamers as rotamers
from gpath.objectives import contacts


def parse_rotamers(chi_groups,node):
	"""
	Slices rotamer angles into a list of lists for all chi angles in a given rotamer.

	Parameters
	----------
	chi_groups : 
	    unparsed list of rotmaer angles
	node : TYPE
	    Description
	
	Returns
	-------
	TYPE
	    Description
	"""
	parsed_chi_groups=[]
	chi_groups = sorted(chi_groups, key=operator.itemgetter(2))
	for residue, group in itertools.groupby(chi_groups, operator.itemgetter(2)):
		g_list=list(group)
		chis=[item[3] for item in g_list]
		parsed_chi_groups.append([residue,chis])

	chi_num=[len(chi_item[1]) for chi_item in parsed_chi_groups]
	it = iter(node['rotamers'])
	sliced =[list(itertools.islice(it, 0, i)) for i in chi_num]

	parsed_chi_groups=[ [residue,node_chis] for (residue,__), node_chis in zip(parsed_chi_groups,sliced) ]
	return parsed_chi_groups

def merge_two_dicts(residue_dict_start, residue_dict_end):
	"""Merge two dictionaries into one.
	
	Parameters
	----------
	residue_dict_start : dict
	residue_dict_end : dict
	
	Returns
	-------
	residue_dict_merged: dict
	"""
	residue_dict_merged = residue_dict_start.copy()
	residue_dict_merged.update(residue_dict_end)
	return residue_dict_merged

def slerp(p0, p1, t):
	"""
	Computes slerp interpolation of the quaternions associated with the rotation matrices of the ligand.
	
	Parameters
	----------
	p0 : list
		Initial quaternion
	p1 : list
		End quaternion
	t : float
		interpolation fraction
	
	Returns
	-------
	list
		interpolated quaternion as list of floats
	"""
	dot = np.dot(p0 / np.linalg.norm(p0), p1 / np.linalg.norm(p1))

	dot_threshold = 0.9995
	if dot > dot_threshold:
		result = p0 + t * (p1 - p0)
		result = result / np.linalg.norm(result)
		result = result / np.linalg.norm(result)
		return result

	omega = np.arccos(dot)
	so = np.sin(omega)
	result = np.sin((1.0 - t) * omega) / so * p0 + np.sin(t * omega) / so * p1
	result = result / np.linalg.norm(result)
	return result

def correct_geometry(v):
	"""

	If needed, identifies expanded angles with a corresponding one in the range (-180,180)
	
	Parameters
	----------
	v : list
		torsion or rotamer angle list that needs to be corrected to the range (-180,180)
	
	Returns
	-------
	list
		corrected list of rotamer angles
	"""

	for v_comp in v:
		if v_comp > 180.0:
			v_comp = v_comp - 360.0
		if v_comp < -180.0:
			v_comp = v_comp + 360.0
	return v


def correct_quaternion_geometry(node1, node2):
	"""
	Ensure that the dot product of two quaterions is always positive as required by the slerp method
	
	Parameters
	----------
	node1 : dict
	node2 : dict
	"""
	if sum([node1_comp * node2_comp for node1_comp, node2_comp in zip(node1['rotation'], node2['rotation'])]) < 0:
		node2['rotation'] = [-node2_comp for node2_comp in node2['rotation']]


def diff_node(node1, node2):
	"""
	Compute the difference of two nodes. Used to obtain the direction of expansion towards the nearest node
	and to compute the distance between two nodes. Also applies the required geometric corrections
	correct_geometry and correct_quaternion_geometry
	
	Parameters
	----------
	node1 : dict
	node2 : dict
	
	Returns
	-------
	diff_dict : dict
		dictionary of node attributes, with values equal to the difference of the attributes of the two nodes
	"""
	diff_dict = {}

	for key in node1.keys():
		if key == 'mode':
			diff_dict[key] = [node2_comp - node1_comp for node2_comp, node1_comp in
										  zip(node2[key],node1[key])]
		if key == 'torsions':
			diff_dict[key] = [node2_comp - node1_comp for node2_comp,
																	  node1_comp in
										  zip(node2[key],node1[key])]
			correct_geometry(diff_dict[key])
		if key == 'rotamers':
			diff_dict[key] = [node2_comp - node1_comp for node2_comp,
																	  node1_comp in
										  zip(node2[key],node1[key])]
			correct_geometry(diff_dict[key])
		if key == 'rotation':
			correct_quaternion_geometry(node1, node2)
			diff_dict[key] = [node2_comp - node1_comp for node2_comp,
																	  node1_comp in
										  zip(node2[key],node1[key])]
		if key == 'translation':
			diff_dict[key] = [node2_comp - node1_comp for node2_comp,
																	  node1_comp in
										  zip(node2[key],node1[key])]
	return diff_dict

def make_norm_dict(node1, node2):
	"""
	Creates a dictionary of the distance of each attribute of two nodes
	Used in order to normalize the expansion vector created in expand()
	
	Parameters
	----------
	node1 : dict
	node2 : dict
	
	Returns
	-------
	norm_dict : dict
		dictionary of distances
	"""
	diff_dict = diff_node(node1, node2)
	norm_dict = {}
	for key in diff_dict.keys():
		if key == 'mode':
			norm_dict[key] = math.sqrt(
				sum([dx_comp * dx_comp for dx_comp in diff_dict[key]]))
		if key == 'torsions':
			norm_dict[key] = math.sqrt(
				sum([dx_comp * dx_comp for dx_comp in diff_dict[key]]))
		if key == 'rotamers':
			norm_dict[key] = math.sqrt(
				sum([dx_comp * dx_comp for dx_comp in diff_dict[key]]))
		if key == 'rotation':
			norm_dict[key] = math.sqrt(
				sum([dx_comp * dx_comp for dx_comp in diff_dict[key]]))
		if key == 'translation':
			norm_dict[key] = math.sqrt(
				sum([dx_comp * dx_comp for dx_comp in diff_dict[key]]))
	return norm_dict


def weighted_norm(node1, node2):
	"""
	Apply weights to the norm dictionary from norm_dict to give equal
	weight to all attributes of the node
	
	Parameters
	----------
	node1 : dict
	node2 : dict
	
	Returns
	-------
	float
		weighted norm between two nodes
	"""

	norm_dict = make_norm_dict(node1, node2)
	weight_norm_dict={}
	for key in norm_dict.keys():
		if key == 'mode':
			weight_norm_dict[key] = norm_dict[key] / 1.0
		if key == 'torsions':
			weight_norm_dict[key] = norm_dict[key] / 360.0
		if key == 'rotamers':
			weight_norm_dict[key] = norm_dict[key] / 360.0
		if key == 'rotation':
			weight_norm_dict[key] = norm_dict[key] / 1.0
		if key == 'translation':
			weight_norm_dict[key] = norm_dict[key] / 1.0

	norm=sum(norm_dict.values())
	return norm


def get_nearest_list_index(node_list, guide_node):
	"""
	Finds nearest nodes among node_list, using the metric given by weighted_norm
	and chooses one of them at random.
	
	Parameters
	----------
	node_list : list
		list of nodes corresponding to one of the two search trees growing towards each other.
	guide_node : dict
		node that has been randomly chosen to expand towards
	
	Returns
	-------
	min_ind : int
		index of the chosen node
	min_dist_choice : float
		distance between the chosen node and the guide_node
	"""
	k_nearest = int(len(node_list) / 100) + 1
	dlist = [weighted_norm(node, guide_node) for node in node_list]
	k_min_dist_list = heapq.nsmallest(k_nearest, dlist)
	min_dist_choice = random.choice(k_min_dist_list)
	min_ind = dlist.index(min_dist_choice)

	return min_ind, min_dist_choice


class RRT:
	"""
	Contains the express algorithm to generate the path and all the objects needed for it
	(the self.protein, the ligand, etc.)
	_start configurations correspond to the state at the frame frame_num, _end at frame frame_num+1
	
	"""
	def __init__(self, input_path, output_path, path_num, frame_num, greed=10, max_iter=10000,
				 relax_clashscore=True, max_clashscore=100.0,step_num=10):#74.96
		"""Summary
		
		Parameters
		----------
		input_path : str
			directory path from which to gather the input files
		output_path : str
			path to the directory in which the result files will be placed
		path_num : int
		    path number of the .zip file in the directory
		frame_num : int
			frame number used as the starting frame. the end frame is taken as frame_num+1
		greed : int, optional
			greed parameter of the RRT-algorithm
		max_iter : int, optional
			maximum number of iterations
		relax_clashscore : bool, optional
			Optional argument to allow steadily increasing clashscore values or use the maximum permitted value from the start
		max_clashscore : float, optional
			value of the maximal permitted clashscore
		step_num : int, optional
			expected number of interpolation steps between the start and end configuration.
		"""
		self.input_path = input_path
		self.output_path = output_path
		self.frame_num = frame_num
		self.greed = greed
		self.max_iter = max_iter
		self.relax_clashscore = relax_clashscore
		self.max_clashscore = max_clashscore
		self.step_num = step_num
		self.path_num=path_num

	def ready(self):
		"""
		
		Attributes
		----------
		clash_diff : float
		    Difference between the clashscore of the start node as measured by the refinement script and the one saved from the Gpath results
		    refinement clashscore can be higher since it may include more rotamers from frame_num and frame_num+1 as part of the interpolation 
		config_end : dict
		    start node configuration corresponding to the state at frame frame_num+1
		config_start : dict
		    start node configuration corresponding to the state at frame frame_num
		contacts_objective : gpath objective
		    
		frame_vdw_volume : float
		    total Van der Wals volume of the ligand and active rotamers 
		goal: bool
			records whether the algorithm has been able to connect the initial and final configurations
		ind : gpath individual
		  
		iteration_counter : int
		    counter that will keep track of the total number of iterations of the RRT-algorithm
		ligand : chimera molecule

		ligand_gene : gpath gene

		next_normal_mode_sample : numpy array
		    atom coordiantes corresponding to the normal mode deformation of the frame framenum+1
		next_sample_number : int
		    normal modes sample number for the frame frame_number+1
		nm_gene : gpath gene
		    Description
		normal_mode_sample : numpy array
		    atom coordinates corresponding to the normal mode deformation of the frame framenum
		path : list
		    list of nodes connecting the start and end configuration
		protein : chimera molecule
		    
		protein_gene : gpath gene
		   
		rotamer_gene : gpath gene
		    
		rotamer_start : list
		    list of rotamer angles with values corresponding to their start configuration
		sample_number : int
		    normal modes sample number for the frame frame_number
		to_zero : tuple
		    matrix that centers the ligand to the origin before applying the second transformation to set it in the right location
		torsion_gene : gpath gene



		"""
		self.config_end = {}
		self.config_start = {}

		modes_path = os.path.join(self.input_path, 'modes.nmd')
		samples_path = os.path.join(self.input_path, 'samples.samples')

		for file in os.listdir(self.input_path):
			if str(self.path_num).zfill(3)+'.zip' in file:
				gpathzip = file

		zip_path = os.path.join(self.input_path, gpathzip)

		with zipfile.ZipFile(zip_path) as myzip:
			for file_name in myzip.namelist():
				if 'Protein' in file_name:
					protein_path = myzip.extract(
						file_name, path=self.output_path)

				if 'Ligand' in file_name:
					ligand_path = myzip.extract(
						file_name, path=self.output_path)

				if '.zip' in file_name:
					pathway_zip_path = myzip.extract(
						file_name, path=self.output_path)

		with zipfile.ZipFile(pathway_zip_path) as myzip:
			allele_path = myzip.extract('allele.txt', self.output_path)
			scores_path = myzip.extract('scores.txt', self.output_path)
			yaml_path=os.path.join(self.input_path,os.path.basename(self.input_path)+'.yaml')
			frame_path = myzip.extract('frame_{:03d}.pdb'.format(self.frame_num), self.output_path)
			framep1_path = myzip.extract('frame_{:03d}.pdb'.format(self.frame_num + 1), self.output_path)

		with open(allele_path, 'r') as allele_file:
			try:
				points=eval(allele_file.read())
			except:
				print('allele_file load error')

		with open(scores_path, 'r') as scores_file:
			try:
				scores=eval(scores_file.read())
			except:
				print('scores_file load error')

		with open(yaml_path, 'r') as yaml_file:
			try:
				settings_dict=yaml.load(yaml_file)
				gene_settings={}
				for gene in settings_dict['genes']:
					if 'torsion' in gene['module']:
						gene_settings['torsions']=gene
						gene_settings['torsions'].pop('module')
						gene_settings['torsions'].pop('name')
						gene_settings['torsions'].pop('target')
						gene_settings['torsions'].pop('anchor')
						continue
					if 'rotamers' in gene['module']:
						gene_settings['rotamers']=gene
						gene_settings['rotamers'].pop('module')
						gene_settings['rotamers'].pop('name')
						continue
					if 'normalmodes' in gene['module']:
						gene_settings['normalmodes']=gene
						gene_settings['normalmodes'].pop('module')
						gene_settings['normalmodes'].pop('name')
						gene_settings['normalmodes'].pop('target')
						continue
			except yaml.YAMLError as exc:
				print(exc)

		if 'torsions' in points:
			torsion_list = points['torsions']
			torsion_anchor = points['torsion_anchor']

		xf_list = zip(points['positions'], points['rotations'])

		self.ind = base.MolecularIndividual()

		self.ligand_gene = molecule.Molecule(path=ligand_path, parent=self.ind)
		self.ind.genes['Ligand'] = self.ligand_gene
		self.ligand_gene.express()
		self.ligand = self.ligand_gene.compound.mol

		or_x, or_y, or_z = np.average(atom_positions(self.ligand.atoms), axis=0)
		self.to_zero = ((1.0, 0.0, 0.0, -or_x),
						(0.0, 1.0, 0.0, -or_y),
						(0.0, 0.0, 1.0, -or_z))

		self.protein_gene = molecule.Molecule(path=protein_path, parent=self.ind)
		self.ind.genes['Protein'] = self.protein_gene
		self.protein_gene.express()
		self.protein = self.protein_gene.compound.mol

		if 'torsions' in points:
			self.torsion_gene = torsion.Torsion('Ligand', parent=self.ligand_gene.parent, anchor=torsion_anchor,**gene_settings['torsions'])

		if 'normal_modes' in points:
			self.nm_gene = path_normalmodes.NormalModes(target='Protein', parent=self.ind, path=modes_path, samples_path=samples_path,**gene_settings['torsions'])
			self.ind.genes['NM'] = self.nm_gene
		else:
			print('No normal modes in allele.txt')

		# parse rotamer information from allele.txt and frame_00x.pdb, if allele.txt used rotamers
		if 'rotamers' in points:
			if points['rotamers'][0] is not []:
				frame = chimera.openModels.open(frame_path)
				frame_protein = frame[0]
				framep1 = chimera.openModels.open(framep1_path)
				framep1_protein = framep1[0]
				residues = points['rotamers']
				chi_init, residue_init = residues[self.frame_num][0], residues[self.frame_num][1]
				chi_end, residue_end = residues[self.frame_num + 1][0], residues[self.frame_num + 1][1]
				residue_dict_start = dict(zip(chi_init, residue_init))
				residue_dict_end = dict(zip(chi_end, residue_end))
				residue_dict_zero = merge_two_dicts(residue_dict_start, residue_dict_end)
				residue_dict_zero = residue_dict_zero.fromkeys(residue_dict_zero, 0.0)

				residue_dict_start = merge_two_dicts(
					residue_dict_zero, residue_dict_start)
				residue_dict_end = merge_two_dicts(
					residue_dict_zero, residue_dict_end)

				residues=[residue.split('/') for residue in residue_dict_start.keys()]
				residues=[[residue[0],int(residue[1])] for residue in residues]
				self.rotamer_gene= rotamers.Rotamers(residues=residues,parent=self.protein_gene.parent,**gene_settings['rotamers'])
				self.ind.genes['Rotamers'] = self.rotamer_gene
				self.ind.__ready__()
				self.rotamer_gene.__ready__()
				self.rotamer_start = []
				rotamer_end = []
				for rot in residue_dict_start.keys():
					res_num = int(rot.split('/')[1])

					residue_dict_start[res_num] = residue_dict_start[rot]
					residue_dict_end[res_num] = residue_dict_end[rot]
					del residue_dict_end[rot]
					del residue_dict_start[rot]

					for rotamer,protein in zip([self.rotamer_start,rotamer_end],[frame_protein,framep1_protein]):		
						residue = protein.findResidue(res_num - 1)
						for chi_num in range(1, 5):
							chi_value = getattr(residue, 'chi' + str(chi_num))
							if chi_value is not None:
								rotamer.append([res_num, chi_num, self.protein.findResidue(res_num - 1), chi_value])

				chimera.openModels.close(frame)
				chimera.openModels.close(framep1)
				self.rotamer_start.sort()
				rotamer_end.sort()

				rotamer_start_values = [item[3] for item in self.rotamer_start]
				rotamer_end_values = [item[3] for item in rotamer_end]
				self.config_start['rotamers'] = rotamer_start_values
				self.config_end['rotamers'] = rotamer_end_values
			else:
				print('No rotamers in frame %s' % self.frame_num)
		else:
			print('No rotamers in allele.txt')

		# initialize self.protein normal modes
		if 'normal_modes' in points:
			self.sample_number = points['normal_modes'][self.frame_num]
			self.next_sample_number = points['normal_modes'][self.frame_num + 1]
			if self.frame_num == 0:
				self.normal_mode_sample = np.asarray(
					self.nm_gene._original_coords)
			else:
				self.normal_mode_sample = np.asarray(
					self.nm_gene.NORMAL_MODES_SAMPLES[self.sample_number])

			self.next_normal_mode_sample = np.asarray(
				self.nm_gene.NORMAL_MODES_SAMPLES[self.next_sample_number])
			
			self.config_start['mode'] = [0.0]
			self.config_end['mode'] = [1.0]

		rotamer_residues = [item[2] for item in self.rotamer_start]
		self.contacts_objective = contacts.Contacts(['Ligand'], which='clashes', rotamer_residues=rotamer_residues)

		self.config_start['rotation'] = list(M.chimera_xform(xf_list[self.frame_num][1]).getQuaternion())
		self.config_end['rotation'] = list(M.chimera_xform(xf_list[self.frame_num + 1][1]).getQuaternion())
		self.config_start['translation'] = [xf_list[self.frame_num][0][i][3] for i in range(3)]
		self.config_end['translation'] = [xf_list[self.frame_num+1][0][i][3] for i in range(3)]
		
		if 'torsions' in points:
			torsion_list[0] = [0.0 for i in range(
				len(self.torsion_gene.rotatable_bonds))]
			self.config_start['torsions'] = torsion_list[self.frame_num]
			self.config_end['torsions'] = torsion_list[self.frame_num + 1]


		self.config_start['clashscore'] = None
		self.config_end['clashscore'] = None

		self.config_start['iteration_counter'] = None
		self.config_end['iteration_counter'] = None

		self.config_start['parent'] = None
		self.config_end['parent'] = None

		self.iteration_counter = -2
		self.collision_check(self.config_start)

		self.iteration_counter = -3
		self.collision_check(self.config_end)

		self.max_clashscore=max(self.max_clashscore,self.config_start['clashscore'],self.config_end['clashscore'])

		self.path=[]

		lig_vdw=[4/3*np.pi*lig_atom.defaultRadius**3 for lig_atom in self.ligand.atoms]
		res_atoms=[]
		for item in self.rotamer_start:
			res_atoms.extend(item[2].atoms)
		res_vdw=[4/3*np.pi*res_atom.defaultRadius**3 for res_atom in res_atoms]
		self.frame_vdw_volume=sum(lig_vdw+res_vdw)

	def save_connect_nodes(self,iteration, node_list, guide_node):
		"""

		Saves the nodes from each RRT-tree nearest to each other
		
		Parameters
		----------
		iteration : int
			current iteration of the RRT-method
		node_list : list
			list of nodes to be expanded on
		guide_node : dict
			
		Returns
		-------
		connect_node_start : dict
			node of the start tree
		connect_node_end : dict
			node of the end tree
		    Description
		"""
		if iteration % 2 == 1:
			connect_node_start = node_list[-1]
			connect_node_end = guide_node
		else:
			connect_node_start = guide_node
			connect_node_end = node_list[-1]

		return 	connect_node_start,connect_node_end

	def expand(self, new_node, guide_node,step_size_dict):
		"""
		expand towards guide_node from nearest_node
		
		Parameters
		----------
		new_node : dict
			new node to be added to the current tree being expanded
		guide_node : dict
			node from the opposite tree that gives the guiding direction
		step_size_dict : dict
		    dictionary of step size values vor each component of the coordinate vector
		"""
		norm_dict = make_norm_dict(new_node, guide_node)
		diff_dict=diff_node(new_node,guide_node)
		new_node_dict=new_node

		for key in new_node_dict.keys():
			if key == 'mode':
				if norm_dict[key] != 0:
					new_node_dict[key] = [new_comp + step_size_dict[key] * diff_comp / norm_dict[key]
						for new_comp, diff_comp in zip(new_node_dict[key], diff_dict[key])]
			if key == 'torsions':
				if norm_dict[key] != 0:
					new_node_dict[key] = [new_comp + step_size_dict[key] * diff_comp / norm_dict[key]
						for new_comp, diff_comp in zip(new_node_dict[key], diff_dict[key])]
					correct_geometry(new_node_dict[key])
			if key == 'rotamers':
				if norm_dict[key] != 0:
					new_node_dict[key] = [new_comp + step_size_dict[key] * diff_comp / norm_dict[key]
						for new_comp, diff_comp in zip(new_node_dict[key], diff_dict[key])]
					correct_geometry(new_node_dict[key])
			if key == 'rotation':
				if norm_dict[key] != 0:
					new_node_dict[key] = list(slerp(np.asarray(new_node_dict[key]), np.asarray(
						guide_node['rotation']), step_size_dict[key] / (norm_dict[key])))
			if key == 'translation':
				if norm_dict[key] != 0:
					new_node_dict[key] = [new_comp + step_size_dict[key] * diff_comp / norm_dict[key]
						for new_comp, diff_comp in zip(new_node_dict[key], diff_dict[key])]

		new_node=new_node_dict

	def apply_node_config(self, node):
		"""
		Apply the rotamers, torsions, normal modes and rotations of the given node to the ligand
		and the protein
		
		Parameters
		----------
		node : dict
		
		"""
		quaternion_vector = node['rotation']
		translation = node['translation']

		if 'torsions' in self.config_start:
			self.torsion_gene.allele = node['torsions']
			self.torsion_gene.express()

		if 'mode' in self.config_start:
			if self.frame_num == 0:
				self.normal_mode_sample = self.nm_gene._original_coords
			else:
				self.normal_mode_sample = np.asarray(
					self.nm_gene.NORMAL_MODES_SAMPLES[self.sample_number])
			self.next_normal_mode_sample = np.asarray(
				self.nm_gene.NORMAL_MODES_SAMPLES[self.next_sample_number])
			interpol_mode = (1.0 - node['mode'][0]) * self.normal_mode_sample + node['mode'][0] * self.next_normal_mode_sample
			self.nm_gene.allele = interpol_mode
			self.nm_gene._need_express = True
			self.nm_gene.express()

		if 'rotamers' in self.config_start:
			parsed_rotamers=parse_rotamers(self.rotamer_start,node)
			for residue,chis in parsed_rotamers:
				if residue.type not in self.rotamer_gene._residues_without_rotamers:
					try:
						self.rotamer_gene.update_rotamer(residue,chis)
					except NoResidueRotamersError:  # ALA, GLY...
						print('no rotamers')
						self.rotamer_gene._residues_without_rotamers.add(residue.type)

		#apply rotation and translation matrix onto the ligand
		rotation_matrix = M.xform_matrix(chimera.Xform().quaternion(*quaternion_vector))
		translation_matrix = M.translation_matrix(translation)
		matrices = (translation_matrix,) + (rotation_matrix,) + (self.to_zero,)
		matrices = M.multiply_matrices(*matrices)
		self.ligand.openState.xform = M.chimera_xform(matrices)

	def collision_check(self, node):
		"""
		Checks if the node surpasses the maximum clashscore value,
		if it doesn't it adds the clashscore value as an attribute of the node
		
		Parameters
		----------
		node : dict
			Node to be checked for collision state
		
		Returns
		-------
		Boolean
			Returns True in case the clashscore value exceeds the the maximum clashscore
			value, False otherwise.
		"""
		self.iteration_counter += 1
		self.apply_node_config(node)
		clashscore = self.contacts_objective.evaluate_clashes(self.ind)
		if self.iteration_counter>0 and clashscore > self.max_clashscore:
			
			return True  # collision
		else:
			node['clashscore'] = clashscore
			node['iteration_counter'] = self.iteration_counter

			return False  # no collision, safe

	def write_refinement(self):
		"""
		Applies the transformation specified for each node in the generated path and writes
		the protein and ligand molecule objects as a trajectoryxxx.pdb file, that will be later
		converted to an MD Traj .pdb file
		"""
		os.chdir(self.output_path)
		for node_num, node in enumerate(self.path):
			self.apply_node_config(node)
			clashscore_written = self.contacts_objective.evaluate_clashes(self.ind)


			totalnum = str(self.frame_num).zfill(3) + str(node_num).zfill(3)
			if node_num==0:
				chimera.pdbWrite([self.ligand, self.protein], chimera.Xform(),'check%s.pdb' % totalnum)
				self.clash_diff=clashscore_written-node['clashscore']

			print(clashscore_written, 'Clashscore of interpol%s.pdb' % totalnum, 'vs stored:',
				  node['clashscore'], node['iteration_counter'])

			combination = combine([self.ligand, self.protein], self.protein)
			chimera.pdbWrite([combination], chimera.Xform(),
							 'trajectory%s.pdb' % totalnum)

			combination.destroy()
		if 'normal_modes' in self.config_start:
			self.nm_gene.unexpress()
		if 'torsions' in self.config_start:
			self.torsion_gene.unexpress()

		self.ligand_gene.unexpress()
		self.protein_gene.unexpress()

	def nodeadder(self):
		"""
		Alternatively expands one of the two RRT trees For a number of max_iter iterations
		
		Returns
		-------
		goal : bool
			whether the algorithm has been able to connect the initial and final node or not
		"""
		minimal_distance = 1000
		self.iteration_counter=0
		connect_node_start = self.config_start
		connect_node_end = self.config_end
		node_list_start = [self.config_start]
		node_list_end = [self.config_end]
		norm_start_end_dict = make_norm_dict(self.config_start, self.config_end)
		norm_start_end_nearest = weighted_norm(self.config_start, self.config_end)
		halt_step = norm_start_end_nearest / self.step_num
		min_clashscore = max(self.config_start['clashscore'],self.config_end['clashscore'])
		relax_stepsize = (self.max_clashscore - min_clashscore) / self.max_iter
		step_size_dict = {}
		for key in norm_start_end_dict.keys():
			step_size_dict[key] = norm_start_end_dict[key] / self.step_num

		center_trans = [np.mean((config_start_comp,config_end_comp)) for config_start_comp, config_end_comp
						in zip(self.config_start['translation'], self.config_end['translation'])]

		node_list = node_list_start
		other_node_list = node_list_end

		for iteration in range(self.max_iter):

			node_list,other_node_list = other_node_list,node_list

			if self.relax_clashscore:
				self.max_clashscore = relax_stepsize * iteration + min_clashscore

			guide_dict = {}
			for key in self.config_start.keys():

				if key == 'mode':
					guide_dict[key] = [random.uniform(-1, 2)]
				if key == 'torsions':
					guide_dict[key] = [
						random.uniform(-180, 180) for i in range(0, len(self.config_end['torsions']))]
				if key == 'rotamers':
					guide_dict[key] = [
						random.uniform(-180, 180) for i in range(0, len(self.config_end['rotamers']))]
				if key == 'rotation':
					guide_dict[key] = list(M.chimera_xform(
						random_rotation()).getQuaternion())
				if key == 'translation':

					guide_dict[key] = [random.uniform(trans_comp - 2 * norm_start_end_dict['translation'],
													  trans_comp + 2 * norm_start_end_dict['translation']) for
									   trans_comp in center_trans]
			guide_node = guide_dict

			[min_ind, node_distance] = get_nearest_list_index(node_list, guide_node)
			nearest_node = node_list[min_ind]

			new_node = copy.deepcopy(nearest_node)

			self.expand(new_node, guide_node,step_size_dict)
			if not self.collision_check(new_node):
				new_node['parent'] = min_ind
				node_list.append(copy.deepcopy(new_node))
				print(node_distance,minimal_distance,halt_step,'node distance')
				if random.randint(0, 100) < self.greed:
					min_ind, min_dist = get_nearest_list_index(other_node_list, new_node)
					guide_node = other_node_list[min_ind]

					self.expand(new_node, guide_node,step_size_dict)
					while not self.collision_check(new_node):

						new_node['parent'] = len(node_list) - 1
						node_list.append(copy.deepcopy(new_node))

						node_distance = weighted_norm(new_node, guide_node)
						print(node_distance,minimal_distance,halt_step,'node distance')
						if weighted_norm(new_node, guide_node) < minimal_distance:
							minimal_distance = node_distance
							print(minimal_distance,halt_step, 'min norm VS reached')
							connect_node_start,connect_node_end=self.save_connect_nodes(iteration, node_list, guide_node)

						if weighted_norm(new_node, guide_node) < halt_step:
							print('goal reached!')
							connect_node_start,connect_node_end=self.save_connect_nodes(iteration, node_list, guide_node)

							return True,connect_node_start,connect_node_end,node_list_start,node_list_end
						self.expand(new_node, guide_node,step_size_dict)

		print('not reached... ')
		return False,connect_node_start,connect_node_end,node_list_start,node_list_end

	def express(self):
		"""
		Function that will construct the random tree and when finished return the found path between two frames
		"""
		self.ready()
		self.goal,connect_node_start,connect_node_end,node_list_start,node_list_end = self.nodeadder()

		path_start = [connect_node_start]
		path_end = [connect_node_end]

		for path,node_list in zip([path_start,path_end],[node_list_start,node_list_end]):
			last_index = path[0]['parent']
			while last_index is not None:
				node = node_list[last_index]
				path.append(node)
				last_index = node['parent']

		path_start.reverse()
		self.path = path_start + path_end
		self.write_refinement()

def Refinement(directory, path_num, frame_start, frame_end):
	"""
	Creates the output path, runs the RRT-algorithm and saves the resulting path .pdb
	files as a MDtraj file
	
	Parameters
	----------
	directory : str
		directory from which the files are 
	path_num : TYPE
	    Description
	frame_start : int

	frame_end : int
	"""

	init_path = os.getcwd()
	input_path = os.path.join(os.getcwd(), 'refinement_input_files', directory)
	output_path = os.path.join(os.getcwd(), 'refinement_output_files', directory,
							   'refinement' + '_' + str(frame_start).zfill(3) + '_' + str(frame_end).zfill(3) + '_' + str(time.time()))

	if os.path.exists(output_path):
		shutil.rmtree(output_path)
		os.makedirs(output_path)
	else:
		os.makedirs(output_path)

	framerange = range(frame_start, frame_end)
	logfile = open(os.path.join(output_path, 'logfile' + '.txt'), 'w+')
	log_dict_list = []
	refined_frame_num=1
	for frame_num in framerange:
		start_time = time.time()
		rrt = RRT(input_path, output_path,path_num, frame_num)
		rrt.express()
		end_time = time.time()

		log_dict = {}
		log_dict['1_frame_num'] = rrt.frame_num
		log_dict['2_goal_reached'] = rrt.goal
		log_dict['3_iterations'] = rrt.iteration_counter
		log_dict['4_time'] = end_time - start_time
		log_dict['5_clash_diff'] = rrt.clash_diff
		log_dict['6_clashscore_iteration'] = [(node,refined_frame_num+i) for i,node in enumerate(rrt.path)]
		log_dict['7_frame_vdw_volume'] = rrt.frame_vdw_volume
		log_dict_list.append(log_dict)
		refined_frame_num+=len(rrt.path)
	logfile.write(pp.pformat(log_dict_list))

	os.chdir(output_path)
	file_list = []
	trajectory_list = [output_file for output_file in os.listdir(
		output_path) if 'trajectory' in output_file]
	trajectory_list.sort(key=lambda f: int(filter(str.isdigit, f)))

	for the_file in trajectory_list:
		file_path = os.path.join(output_path, the_file)
		file_list.append(file_path)
	md_trajectory = mdtraj.load(file_list)
	md_trajectory.save(os.path.join(output_path, 'refinement.pdb'))

	for the_file in trajectory_list:
		file_path = os.path.join(output_path, the_file)
		try:
			if os.path.isfile(file_path):
				os.unlink(file_path)
		except Exception as e:
			print(e)
	os.chdir(init_path)
	chimera.closeSession()

def parse_cli():
	"""
	Cli input parser
	"""
	p = ArgumentParser()
	p.add_argument('directory', type=str,
				   help='directory containing the results from the Gpath inside refinement_test_input_files')
	p.add_argument('path_num', type=int,
				   help='path number of the .zip file in the directory')
	p.add_argument('frame_start', type=int,
				   help='frame number from which the interpolation has to be done')
	p.add_argument('frame_end', type=int,
				   help='frame number until which the interpolation has to be done')
	return p.parse_args()

if __name__ == '__main__':

	args = parse_cli()
	directory=args.directory
	frame_start=args.frame_start
	frame_end=args.frame_end
	path_num=args.path_num

	output_path=os.path.join(os.getcwd(),'examples','refinement_output_files',directory)
	Refinement(directory,path_num,frame_start,frame_end)
