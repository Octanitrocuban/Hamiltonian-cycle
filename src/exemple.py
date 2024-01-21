# -*- coding: utf-8 -*-
"""
This script contain some example about the use of functions from
hamiltonian.py
"""
import numpy as np
import matplotlib.pyplot as plt
import hamiltonian as hlt
from time import time
from tqdm import tqdm
from scipy.spatial.distance import cdist
#=============================================================================
# Algorithm idea from Code Bullet on snake ia:
# https://www.youtube.com/watch?v=tjQIO1rqTBE itself inspired by JOHN TAPSELL:
# https://johnflux.com/2015/05/02/nokia-6110-part-3-algorithms/

# fill to_do with 'random_cycle' and/or 'deterministic_cycle'.
to_do = ['']

if 'random_cycle' in to_do:
	n_lines, n_columns = 100, 100
	return_min_tree = False

	if return_min_tree:
		outs = hlt.random_hamiltonian_cycle(n_lines, n_columns,
										    return_min_tree)
		
		cycle, nodes_r, min_tree, cents = outs

	else:
		cycle, nodes_r = hlt.random_hamiltonian_cycle(n_lines, n_columns,
												      return_min_tree)
	
	plt.figure(figsize=(15, 15))
	plt.plot(nodes_r[:, 0], nodes_r[:, 1], '.', label='nodes')
	if return_min_tree:
		for i in range(len(min_tree)):
			for j in range(1, len(min_tree[i])):
				if (i == 0)&(j == 1):
					plt.plot([cents[min_tree[i][0], 0],
					          cents[min_tree[i][j], 0]],
							 [cents[min_tree[i][0], 1],
							  cents[min_tree[i][j], 1]],
							  'r-', label='minimum\nspanning\ntree wall')
				else:
					plt.plot([cents[min_tree[i][0], 0],
					          cents[min_tree[i][j], 0]],
							 [cents[min_tree[i][0], 1],
							  cents[min_tree[i][j], 1]],
							  'r-')
	
	plt.plot(nodes_r[cycle, 0], nodes_r[cycle, 1], 'm', label='path')
	plt.xlim(-1, n_columns)
	plt.ylim(-1, n_lines)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.legend(loc=[1.01, 0.5], fontsize=16, markerscale=2)
	plt.show()

elif 'deterministic_cycle' in to_do:
	n_lines, n_columns = 100, 101
	cycle_d, nodes = hlt.deterministic_cycle(n_lines, n_columns)

	plt.figure(figsize=(15, 15))
	plt.plot(nodes[:, 0], nodes[:, 1], '.', label='nodes')
	plt.plot(nodes[cycle_d, 0], nodes[cycle_d, 1], 'm', label='path')
	plt.xlim(-1, n_columns+1)
	plt.ylim(-1, n_lines+1)
	plt.xticks(fontsize=16)
	plt.yticks(fontsize=16)
	plt.legend(loc=[1.01, 0.5], fontsize=16, markerscale=2)
	plt.show()
