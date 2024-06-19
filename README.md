# Hamiltonian-cycle
Module with functions to create random and deterministic hamiltonian cycle.


### Scripts

**hamiltonian.py**: is the module with the functions to create random and determinit hamiltonian cycle for a regular (square or rectangular) mesh.

  - kurskal: function to construct the minimum spanning tree through kurskal algorithm.
  - make_square_mesh: function to construct the nodes mesh position.
  - tree_square_mesh: function to construct the tree of connections between nodes mesh position.
  - derand_centres: function to interpolate the nodes positions that will be used as wall.
  - tree_nodes_centers: function to construct the connection tree between the nodes that are from the cycle and the nodes that are from the walls.
  - make_path: function to construct the random hamiltonian cycle from the constructed walls. It use the 'right hand method' which usually serves yo solve maze. 
  - random_hamiltonian_cycle: function to create -if it is possible- a random hamiltonian cycle path.
  - deterministic_cycle: function to create a hamiltionian cycle through a determinisc method.

**exemple.py**: is the script that contain some example about the use of functions from hamiltonian.py

  - if 'random_cycle' is in the list to_do: create and plot a random hamiltonian cycle with the given parameters.
  - if 'deterministic_cycle' is in the list to_do: create and plot a deterministic hamiltonian cycle with the given parameters.

### Examples

Here are some examples of what we can have for random and deterministic method.

A 20 by 20 random cycle without the plot of the walls.
![Exemple picture](img/random_20_20_no_wall.png)


A 20 by 20 random cycle with the plot of the walls.
![Exemple picture](img/random_20_20.png)


A 100 by 100 random cycle without the plot of the walls.
![Exemple picture](img/random_100_100_no_wall.png)


A 100 by 100 random cycle with the plot of the walls.
![Exemple picture](img/random_100_100.png)


A 20 by 21 deterministic cycle.
![Exemple picture](img/determinist_20_21.png)


A 101 by 100 deterministic cycle.
![Exemple picture](img/determinist_101_100.png)

### Statistics
As we can see on the two following histograms, the deterministic algorithm is much faster than the random one. This come from the fact that the first algorithm have only one step to create the cycle and use only one loop. On the contrary the second algorithm have at least nine step with four heavy loop. Also as we can see on the third plot the time consumption of some steps are really heavy. The two heaviest steps are: first the construction of the cycle with the 'right hand' method and secondly the construction of the minimum spaning tree through the kurskal algorithm.

![Exemple picture](img/time_conso_pdf_determ_101_100_runs_100000.png)
![Exemple picture](img/time_conso_pdf_random_100_100_runs_1000.png)
![Exemple picture](img/time_conso_elem_100_100.png)
