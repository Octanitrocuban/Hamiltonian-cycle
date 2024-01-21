#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contain the functions to create random and determinit hamiltonian
cycle for a regular (square or rectangular) mesh.
"""
import numpy as np
from scipy.spatial.distance import cdist
#==============================================================================
def kurskal(dots):
	"""
	Function to construct the minimum spanning tree through kurskal algorithm.

	Parameters
	----------
	dots : numpy.ndarray
		Dots we want to connect in the minimal tree.

	Returns
	-------
	tree : numpy.ndarray
		Vectorized dictionnary of the dots-connections.

	Example
	-------
	In [0] : kurskal(half_of_centers(centers)
					+ rd.uniform(-.1, .1, half_cents.shape))

	Out [0]: [np.array([0, 1]), np.array([1, 4, 0]), np.array([2, 5]),
			  np.array([3, 4, 6]), np.array([4, 1, 5, 3]),
			  np.array([5, 4, 8, 2]), np.array([6, 7, 3]), np.array([7, 6]),
			  np.array([8, 5])]

	"""
	# calculates the distance matrix
	m_dist = cdist(dots, dots, metric='euclidean').T
	length = len(dots)
	# list of array
	tree = list(np.arange(length)[:, np.newaxis])
	mask = (np.arange(length)-np.arange(length)[:, np.newaxis]) > 0
	# lists of index matrices
	indices = list(np.meshgrid(range(length), range(length)))
	# vector 1d to track connections in the tree and avoid loop formation
	state = np.arange(length)
	# We flatten the 2d matrix by keeping less than half of the distance
	# matrix not to re-evaluate relationships between pairs of points.
	sort_d = m_dist[mask]
	# The same is done for index matrices
	p_j = indices[0][mask]
	p_i = indices[1][mask]
	# Indices triés par ordre croissant selon les valeurs des distances
	rank = np.argsort(sort_d)
	# Sorting indices and distance values
	p_i = p_i[rank]
	p_j = p_j[rank]
	sort_d = sort_d[rank]
	for i in range(len(sort_d)):
		# To have no recontection with loops in the tree
		if state[p_i[i]] != state[p_j[i]]:
			tree[p_i[i]] = np.append(tree[p_i[i]], p_j[i])
			tree[p_j[i]] = np.append(tree[p_j[i]], p_i[i])
			# Update of the 'state' vector
			minima = np.min([state[p_i[i]], state[p_j[i]]])
			state[state == state[p_i[i]]] = minima
			state[state == state[p_j[i]]] = minima
			# early stoping to avoid useless loop
			if len(state[state != minima]) == 0:
				break

	return tree

def make_square_mesh(n, m, step=1):
	"""
	Function to construct the nodes mesh position.

	Parameters
	----------
	n : int
		Number of lines in the mesh.
	m : int
		Number of columns in the mesh.

	Returns
	-------
	nodes : list
		List of x and y positions of the nodes.

	Example
	-------
	In [0] : make_square_mesh(6, 6)
	Out [0]: [array([[0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5],
					 [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5],
					 [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5]]),
			  array([[0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1],
					 [2, 2, 2, 2, 2, 2], [3, 3, 3, 3, 3, 3],
					 [4, 4, 4, 4, 4, 4], [5, 5, 5, 5, 5, 5]])]

	"""
	nodes = np.meshgrid(range(0, n, step), range(0, m, step))
	return nodes

def tree_square_mesh(n, m):
	"""
	Function to construct the tree of connections between nodes mesh
	position.

	Parameters
	----------
	n : int
		Number of lines in the mesh.
	m : int
		Number of columns in the mesh.

	Returns
	-------
	nodes : numpy.ndarray
		A object type numpy array. It list the connections between nodes.
		The i-th nodes is connected to the nodes numbers in the array at the
		i-th position.

	Example
	-------
	In [0] : tree_square_mesh(6, 6)
	Out [0]: np.array([np.array([6, 1]), np.array([0, 7, 2]),
					   np.array([1, 8, 3]), np.array([2, 9, 4]),
					   np.array([3, 10, 5]), np.array([4, 11]),
					   np.array([0, 12, 7]), np.array([6, 1, 13, 8]),
					   np.array([7, 2, 14, 9]), np.array([8, 3, 15, 10]),
					   np.array([9, 4, 16, 11]), np.array([10, 5, 17]),
					   np.array([6, 18, 13]), np.array([12, 7, 19, 14]),
					   np.array([13, 8, 20, 15]), np.array([14, 9, 21, 16]),
					   np.array([15, 10, 22, 17]), np.array([16, 11, 23]),
					   np.array([12, 24, 19]), np.array([18, 13, 25, 20]),
					   np.array([19, 14, 26, 21]), np.array([20, 15, 27, 22]),
					   np.array([21, 16, 28, 23]), np.array([22, 17, 29]),
					   np.array([18, 30, 25]), np.array([24, 19, 31, 26]),
					   np.array([25, 20, 32, 27]), np.array([26, 21, 33, 28]),
					   np.array([27, 22, 34, 29]), np.array([28, 23, 35]),
					   np.array([24, 31]), np.array([30, 25, 32]),
					   np.array([31, 26, 33]), np.array([32, 27, 34]),
					   np.array([33, 28, 35]), np.array([34, 29])])

	"""
	# create 2d array from 0 to n*m
	square = np.arange(n*m).reshape((m, n))
	# cropping
	square_rec = np.zeros((square.shape[0]+2, square.shape[1]+2), dtype=int)-1
	square_rec[1:-1, 1:-1] = square
	nod_sq = make_square_mesh(n, m)
	nod_sq = np.array([np.ravel(nod_sq[0]), np.ravel(nod_sq[1])]).T+1
	around = nod_sq+np.array([[[-1, 0]], [[0, -1]], [[0, 1]], [[1, 0]]])
	# making the tree
	connexions = square_rec[around[:, :, 1], around[:, :, 0]].T
	arbre = []
	for i in range(len(connexions)):
		arbre.append(connexions[i][connexions[i] != -1])

	arbre = np.array(arbre, dtype=object)
	return arbre

def tree_mesh_unit(dots):
	"""
	Function to construct the nodes mesh tree.

	Parameters
	----------
	dots : numpy.ndarray
		Array listing the position of the dots.

	Returns
	-------
	arbre : numpy.ndarray
		A 2d array listing the indices of the connected dots. The shape is
		(N, 2) where N is the number of connections. For a connexion between
		two dots it is indicate by : ([indice of the first dot, indice of
		the second dot]).

	Example
	-------
	In [0] : tree_mesh_unit(grid_dots(6, 6))
	Out [0]: np.array([[0, 1], [0, 6], [1, 2], [1, 7], [2, 3], [2, 8],
					   [3, 4], [3, 9], [4, 5], [4, 10], [5, 11], [6, 7],
					   [6, 12], [7, 8], [7, 13], [8, 9], [8, 14], [9, 10],
					   [9, 15], [10, 11], [10, 16], [11, 17], [12, 13],
					   [12, 18], [13, 14], [13, 19], [14, 15], [14, 20],
					   [15, 16], [15, 21], [16, 17], [16, 22], [17, 23],
					   [18, 19], [18, 24], [19, 20], [19, 25], [20, 21],
					   [20, 26], [21, 22], [21, 27], [22, 23], [22, 28],
					   [23, 29], [24, 25], [24, 30], [25, 26], [25, 31],
					   [26, 27], [26, 32], [27, 28], [27, 33], [28, 29],
					   [28, 34], [29, 35], [30, 31], [31, 32], [32, 33],
					   [33, 34], [34, 35]])

	"""
	# compute the distance matrix
	distances = cdist(dots, dots)
	# where the value is ~1 but not 0
	masque = (distances <= 1)&(distances > 0)
	arbre = np.argwhere(masque)
	arbre = np.unique(np.sort(arbre, axis=1), axis=0)
	return arbre

def derand_centres(minimal_tree, dots):
	"""
	Function to construct the nodes mesh position.

	Parameters
	----------
	minimal_tree : numpy.ndarray
		Vectorized dictionnary of the dots-connections.
	dots : numpy.ndarray
		Dots we want to use for interpolation between the nods connexion.

	Returns
	-------
	uses_centers : numpy.ndarray
		Dots list with the intrepolation between the connected dots.

	Example
	-------
	In [0] : derand_centres(kurskal(
					make_square_mesh(n_lines-1, n_columns-1, 2)
					+ rd.uniform(-.1, .1, half_cents.shape)),
				make_square_mesh(n_lines-1, n_columns-1, 2))

	Out [0]: np.array([[0.5, 0.5], [0.5, 2.5], [0.5, 4.5], [1.5, 0.5],
					   [1.5, 2.5], [1.5, 4.5], [2.5, 0.5], [2.5, 1.5],
					   [2.5, 2.5], [2.5, 3.5], [2.5, 4.5], [3.5, 0.5],
					   [4.5, 0.5], [4.5, 1.5], [4.5, 2.5], [4.5, 3.5],
					   [4.5, 4.5]])

	"""
	uses_centers = []
	for i in range(len(minimal_tree)):
		for j in range(1, len(minimal_tree[i])):
			# first dot
			uses_centers.append(
				[dots[minimal_tree[i][0], 0],
				 dots[minimal_tree[i][0], 1]])

			# interpolated dot
			uses_centers.append(
				[(dots[minimal_tree[i][0], 0]+
				  dots[minimal_tree[i][j], 0])/2,
				 (dots[minimal_tree[i][0], 1]+
				  dots[minimal_tree[i][j], 1])/2])

			# final dot
			uses_centers.append(
				[dots[minimal_tree[i][j], 0],
				 dots[minimal_tree[i][j], 1]])

	# to not have many time the same dot
	uses_centers = np.unique(np.array(uses_centers), axis=0)
	return uses_centers

def tree_nodes_centers(nodes, centers):
	"""
	Function to construct the nodes mesh position.

	Parameters
	----------
	nodes : numpy.ndarray
		Nodes of a mesh that have gride shape.
	centers : numpy.ndarray
		Nodes of a mesh that have gride shape where nodes are at the centers
		of the sub square of 'nodes'.

	Returns
	-------
	tree : numpy.ndarray
		A 2d array listing the indices of the connected dots. The shape is
		(N, 2) where N is the number of connections. For a connexion between
		two dots it is indicate by : ([indice of the nodes, indice of the
		centers]).

	Example
	-------
	In [0] : nodes, centers
	Out [0] : (np.array([[0 0], [1 0], [2 0], [3 0], [4 0], [5 0], [0 1],
						 [1 1], [2 1], [3 1], [4 1], [5 1], [0 2], [1 2],
						 [2 2], [3 2], [4 2], [5 2], [0 3], [1 3], [2 3],
						 [3 3], [4 3], [5 3], [0 4], [1 4], [2 4], [3 4],
						 [4 4], [5 4], [0 5], [1 5], [2 5], [3 5], [4 5],
						 [5 5]]),
			   np.array([[0.5, 0.5], [0.5, 1.5], [0.5, 2.5], [0.5, 3.5],
						 [0.5, 4.5], [1.5, 2.5], [1.5, 4.5], [2.5, 0.5],
						 [2.5, 1.5], [2.5, 2.5], [2.5, 4.5], [3.5, 0.5],
						 [3.5, 4.5], [4.5, 0.5], [4.5, 1.5], [4.5, 2.5],
						 [4.5, 4.5]]))

	In [2] : tree_nodes_centers(nodes, centers)
	Out [2]: np.array([[0, 1], [1, 2], [2, 4], [3, 5], [4, 7], [5, 9],
					   [6, 10], [7, 8], [7, 11], [8, 9], [9, 12], [10, 13],
					   [11, 15], [12, 16], [13, 14], [14, 15]])

	"""
	# compute the distance matrix
	distances = cdist(nodes, centers)
	# where the value is <= 2**0.5
	masque = (distances <= (2**0.5))
	tree = np.argwhere(masque)
	return tree

def make_path(n_lines, n_columns, centres, raveled_nodes):
	"""
	Function to construct the nodes mesh position.

	Parameters
	----------
	n_lines : int
		Number of lines.
	n_columns : int
		Number of n_columns.
	centres : numpy.ndarray
		Centers of squares drawn by nodes and use for construction of
		kurskal tree. Shape = (n, 2).
	raveled_nodes : numpy.ndarray
		Nodes of a mesh that have gride shape.

	Returns
	-------
	chemin : numpy.ndarray
		Vector 1d which contains the order in which the nodes are seen.

	Example
	-------
	# with 4 by 4 nodes
	In [0] : n_lines, n_columns, centres, raveled_nodes
	Out [0]: (4, 4, np.array([[0.5, 0.5], [0.5, 2.5], [1.5, 0.5], [1.5, 2.5],
						      [2.5, 0.5], [2.5, 1.5], [2.5, 2.5]]),
			  np.array([[0, 0], [1, 0], [2, 0], [3, 0], [0, 1], [1, 1],
						[2, 1], [3, 1], [0, 2], [1, 2], [2, 2], [3, 2],
						[0, 3], [1, 3], [2, 3], [3, 3]]))

	In [1] : make_path(n_lines, n_columns, centres, raveled_nodes)
	Out [1]: np.array([0, 1, 2, 3, 7, 6, 5, 9, 10, 11, 15, 14, 13, 12, 8, 4])

	"""
	# map where the algo will 'walk'
	mapp = np.zeros((2*n_columns-1, 2*n_lines-1), dtype='int8')
	posi_centers_t2 = (centres*2).astype(int)
	mapp[posi_centers_t2[:, 1], posi_centers_t2[:, 0]] = 1
	# cropping the map for debug
	recadre = np.zeros((2*n_columns+1, 2*n_lines+1), dtype='int8')
	recadre[1:-1, 1:-1] = mapp
	# starting position
	self_x = 1
	self_y = 1
	start = [self_x, self_y]
	# starting direction
	direction = 'e'
	path = []
	path.append([self_x, self_y])
	while True:
		# the four corners
		corn_1 = recadre[self_x-1, self_y-1]
		corn_2 = recadre[self_x-1, self_y+1]
		corn_3 = recadre[self_x+1, self_y+1]
		corn_4 = recadre[self_x+1, self_y-1]

		# if direction is east => can go {'south', 'north', 'east'}
		if direction == 'e':
			if corn_3 == 0:
				self_x += 2
				direction = 's'
			elif (corn_2 == 1)&(corn_3 == 1):
				self_x -= 2
				direction = 'n'
			elif corn_3 == 1:
				self_y += 2
			else:
				raise

		# if direction is south => can go {'west', 'east', 'south'}
		elif direction == 's':
			if corn_4 == 0:
				self_y -= 2
				direction = 'o'
			elif (corn_4 == 1)&(corn_3 == 1):
				self_y += 2
				direction = 'e'
			elif corn_4 == 1:
				self_x += 2
			else:
				raise

		# if direction is west => can go {'north', 'south', 'west'}
		elif direction == 'o':
			if corn_1 == 0:
				self_x -= 2
				direction = 'n'
			elif (corn_1 == 1)&(corn_4 == 1):
				self_x += 2
				direction = 's'
			elif corn_1 == 1:
				self_y -= 2
			else:
				raise

		# if direction is north => can go {'east', 'west', 'north'}
		elif direction == 'n':
			if corn_2 == 0:
				self_y += 2
				direction = 'e'
			elif (corn_1 == 1)&(corn_2 == 1):
				self_y -= 2
				direction = 'o'
			elif corn_2 == 1:
				self_x -= 2
			else:
				raise

		else:
			raise

		# stop when reach the start
		if [self_x, self_y] == start:
			break
		else:
			path.append([self_x, self_y])

	# recompute for the different shape between a matrix and plot method
	path = np.roll((np.array(path)-1)/2, 1, axis=1)
	# compute the distance matrix
	distances = cdist(path, raveled_nodes.astype(float))
	path = np.argwhere(distances < 1)[:, 1]
	return path

def random_hamiltonian_cycle(n_lines, n_columns, return_min_tree=False):
	"""
	Function to create -if it is possible- a random hamiltonian cycle path.

	Parameters
	----------
	n_lines : int
		Number of rows. It must be >= 2 and an even value.
	n_columns : int
		Number of columns. It must be >= 2 and an even value.
	return_min_tree : bool, optional
		If we want to get the minimum spaning tree. The default is False.

	Raises
	------
	ValueError
		The number of rows and columns must be even.

	Returns
	-------
	cycle : numpy.ndarray
		Vector 1d which contains the order in which the nodes are seen.
	nodes_ravel : numpy.ndarray
		Array with the nodes used.
	minimum_tree : numpy.ndarray, optional
		Array with the dots indices for connections.
	use_centers : numpy.ndarray, optional
		 Dots used to create the minimum spanning tree used as a "wall".

	Example
	-------
	In [0] : hlt.random_hamiltonian_cycle(6, 6, False)
	Out [0]: (np.array([ 0,  1,  7, 13, 14,  8,  2,  3,  4,  5, 11, 10,  9,
						15, 21, 27, 28, 22, 16, 17, 23, 29, 35, 34, 33, 32,
						31, 30, 24, 25, 26, 20, 19, 18, 12,  6], dtype=int64),
			 (np.array([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0],
				        [0, 1], [1, 1], [2, 1], [3, 1], [4, 1], [5, 1],
						[0, 2], [1, 2], [2, 2], [3, 2], [4, 2], [5, 2],
						[0, 3], [1, 3], [2, 3], [3, 3], [4, 3], [5, 3],
						[0, 4], [1, 4], [2, 4], [3, 4], [4, 4], [5, 4],
						[0, 5], [1, 5], [2, 5], [3, 5], [4, 5], [5, 5]]))

	Note
	----
	If return_min_tree is used, it will have to recompute the minimum spaning
	tree and thus take longer.
	The part taking the more time are make_path and kurskal.

	"""
	if ((n_lines%2) == 1)|((n_columns%2) == 1):
		raise ValueError(
			"The number of rows and columns must be even.")

	# create the nodes mesh position to travel
	nodes = make_square_mesh(n_lines, n_columns)
	nodes_ravel = np.array([np.ravel(nodes[0]),
							np.ravel(nodes[1])]).T

	# create the tree of connections between nodes mesh position
	conect_nods = tree_square_mesh(n_lines, n_columns)
	# create the nodes mesh position to create the walls
	half_centers = make_square_mesh(n_lines-1, n_columns-1, 2)
	half_centers = np.array([np.ravel(half_centers[0]),
						     np.ravel(half_centers[1])]).T +.5

	half_centers_random = half_centers + np.random.uniform(-.1, .1,
														   half_centers.shape)

	# compute the walls
	minimum_tree = kurskal(half_centers_random)
	# interpolation to have continuous walls
	use_centers = derand_centres(minimum_tree, half_centers)
	# nodes connection betw
	cycle = make_path(n_lines, n_columns, use_centers, nodes_ravel)
	if return_min_tree:
		minimum_tree = kurskal(use_centers)
		return cycle, nodes_ravel, minimum_tree, use_centers

	else:
		return cycle, nodes_ravel

def deterministic_cycle(num_lines, num_cols):
	"""
	Function to create a hamiltionian cycle through a determinisc method.

	Parameters
	----------
	num_lines : int
		Number of lines in the graph.
	num_cols : int
		Number of columns in the graph.

	Returns
	-------
	path : numpy.ndarray
		Array listing the position followed by the hamiltonian cycle path.
	nodes : numpy.ndarray
		Array listing the position of nodes.

	Raises
	-------
	ValueError
		The number of rows and columns cannot be odd at the same time.

	Example
	-------
	In [0] : deterministic_cycle(5, 6)
	Out [0]: (np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
					    13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
						26, 27, 28, 29, 30, 31, 32, 33, 34, 35]),
		      np.array([[0, 0], [1, 0], [2, 0], [3, 0], [4, 0], [5, 0],
						[5, 1], [5, 2], [5, 3], [5, 4], [5, 5], [4, 5],
						[4, 4], [4, 3], [4, 2], [4, 1], [3, 1], [3, 2],
						[3, 3], [3, 4], [3, 5], [2, 5], [2, 4], [2, 3],
						[2, 2], [2, 1], [1, 1], [1, 2], [1, 3], [1, 4],
						[1, 5], [0, 5], [0, 4], [0, 3], [0, 2], [0, 1]]))

	In [1] : deterministic_cycle(4, 5)
	Out [1]: (np.array([[0, 0], [0, 1], [0, 2], [0, 3], [1, 3], [2, 3],
					    [3, 3], [4, 3], [5, 3], [5, 2], [4, 2], [3, 2],
						[2, 2], [1, 2], [1, 1], [2, 1], [3, 1], [4, 1],
						[5, 1], [5, 0], [4, 0], [3, 0], [2, 0], [1, 0]]),
		      np. array([[0, 0], [0, 1], [0, 2], [0, 3], [1, 0], [1, 1],
				         [1, 2], [1, 3], [2, 0], [2, 1], [2, 2], [2, 3],
						 [3, 0], [3, 1], [3, 2], [3, 3], [4, 0], [4, 1],
						 [4, 2], [4, 3]]))

	"""
	if ((num_lines%2) == 1)&((num_cols%2) == 1):
		raise ValueError(
			"The number of rows and columns cannot be odd at the same time")

	if ((num_lines%2) == 0):
		nodes = np.arange(num_lines)
		nodes = np.array([np.zeros(num_lines, dtype=int), nodes]).T
		for i in range(0, num_lines, 1):
			if (i%2) == 0:
				longueur = np.arange(1, num_cols+1)
			else:
				longueur = np.arange(num_cols, 0, -1)

			hauteur = np.zeros(num_cols, dtype=int)+(num_lines-i-1)
			sub_path = np.array([longueur, hauteur]).T
			nodes = np.concatenate((nodes, sub_path))

	else:
		nodes = np.arange(num_cols)
		nodes = np.array([nodes, np.zeros(num_cols, dtype=int)]).T
		for i in range(0, num_cols, 1):
			if (i%2) == 0:
				hauteur = np.arange(1, num_lines+1)
			else:
				hauteur = np.arange(num_lines, 0, -1)

			longueur = np.zeros(num_lines, dtype=int)+(num_cols-i-1)
			sub_path = np.array([longueur, hauteur]).T
			nodes = np.concatenate((nodes, sub_path))

	path = np.arange(len(nodes))
	return path, nodes
